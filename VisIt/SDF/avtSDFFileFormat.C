/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtSDFFileFormat.C                           //
// ************************************************************************* //

#include <avtSDFFileFormat.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <avtParallel.h>

#include <string>

#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>

#include <avtDatabaseMetaData.h>
#include <avtIntervalTree.h>
#include <avtVariableCache.h>
#include <avtStructuredDomainBoundaries.h>
#include <avtMaterial.h>
#include <avtSpecies.h>
#include <DebugStream.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidFilesException.h>
#include <InvalidVariableException.h>

#include <dlfcn.h>

using     std::string;
int avtSDFFileFormat::extension_not_found = 0;


// ****************************************************************************
//  Memory management stack
// ****************************************************************************

struct stack;

struct stack {
    sdf_block_t *block;
    struct stack *next;
};

static struct stack *stack_head = NULL;
static struct stack *stack_tail = NULL;

#define MAX_MEMORY 2147483648 // 2GB
static uint64_t memory_size = 0;


static inline void stack_alloc(sdf_block_t *b)
{
    struct stack *tail;
    b->data = calloc(b->nlocal, b->type_size_out);
    memory_size += b->nlocal * b->type_size_out;
    stack_tail->next = tail = (struct stack*)malloc(sizeof(struct stack));
    tail->block = b;
    tail->next = NULL;
    stack_tail = tail;
}


static inline void stack_free(void)
{
    sdf_block_t *b;
    struct stack *head;

    if (memory_size < MAX_MEMORY) return;

    while (stack_head->next) {
        head = stack_head;
        stack_head = stack_head->next;
        free(head);
        b = stack_head->block;
        stack_head->block = NULL;
        free(b->data);
        b->data = NULL;
        b->done_data = 0;
        memory_size -= b->nlocal * b->type_size_out;
        if (memory_size < MAX_MEMORY) break;
    }
}


static inline void stack_init(void)
{
    if (!stack_head) stack_head =
        stack_tail = (struct stack*)calloc(1, sizeof(struct stack));
}


sdf_extension_t *avtSDFFileFormat::sdf_extension_load(sdf_file_t *h)
{
    if (avtSDFFileFormat::extension_not_found) return NULL;

    sdf_extension_handle = dlopen("sdf_extension.so", RTLD_LAZY);

    if (!sdf_extension_handle) {
        avtSDFFileFormat::extension_not_found = 1;
        cerr << dlerror() << endl;
        return NULL;
    }

    sdf_extension_create_t *sdf_extension_create =
        (sdf_extension_create_t *)dlsym(sdf_extension_handle,
        "sdf_extension_create");

    sdf_extension_t *ext = sdf_extension_create(h);

    if (!ext) avtSDFFileFormat::extension_not_found = 1;

    return ext;
}


void avtSDFFileFormat::sdf_extension_unload(void)
{
    if (!sdf_extension_handle) return;

    sdf_extension_destroy_t *sdf_extension_destroy =
        (sdf_extension_destroy_t *)dlsym(sdf_extension_handle,
        "sdf_extension_destroy");

    sdf_extension_destroy(ext);

    dlclose(sdf_extension_handle);

    sdf_extension_destroy = NULL;
    ext = NULL;

    return;
}


// ****************************************************************************
//  Helper routine for opening SDF files.
// ****************************************************************************

void
avtSDFFileFormat::OpenFile(int open_only)
{
    if (!h) h = sdf_open(filename, rank, comm, 0);
    if (!h) EXCEPTION1(InvalidFilesException, filename);
    h->use_float = use_float;
    step = h->step;
    time = h->time;
    debug1 << "avtSDFFileFormat:: " << __LINE__ << " h:" << h << endl;

    if (open_only) {
        // Retrieve the extended interface library from the plugin manager
        ext = sdf_extension_load(h);

        sdf_close(h);
        h = NULL;
        return;
    }

    if (h->blocklist) {
        if (ext) ext->timestate_update(ext, h);
        return;
    }

    sdf_read_blocklist(h);
    // Append derived data to the blocklist using built-in library.
    sdf_add_derived_blocks(h);

    if (ext) {
        char **preload;

        preload = ext->preload(ext, h);
        // For each entry in the preload array, try to find the block
        // and populate its data.
        if (preload) {
            sdf_block_t *var, *cur;
            int n = 0;
            cur = h->current_block;
            while(preload[n]) {
                var = sdf_find_block_by_id(h, preload[n]);
                if (var && !var->data) {
                    h->current_block = var;
                    sdf_read_data(h);
                }
                free(preload[n]);
                n++;
            }
            free(preload);
            h->current_block = cur;
        }

        // Append derived data to the blocklist using the extension library.
        ext->read_blocklist(ext, h);
    }
}


// ****************************************************************************
//  Method: avtSDFFileFormat constructor
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

avtSDFFileFormat::avtSDFFileFormat(const char *filename,
    DBOptionsAttributes *readOpts) : avtSTMDFileFormat(&filename, 1)
{
    debug1 << "avtSDFFileFormat::avtSDFFileFormat(filename:" << filename
           << ") " << this << endl;
    // INITIALIZE DATA MEMBERS
#ifdef PARALLEL
    comm = VISIT_MPI_COMM;
    debug1 << "avtSDFFileFormat:: parallel" << endl;
#else
    debug1 << "avtSDFFileFormat:: serial" << endl;
#endif
    rank = PAR_Rank();
    ncpus = PAR_Size();
    ndomains = ncpus;
    debug1 << rank << " " << ncpus << " " << comm << " f:" << filename << endl;
    this->filename = new char[strlen(filename)+1];
    memcpy(this->filename, filename, strlen(filename)+1);
    gotMetadata = false;
    h = NULL;
    sdf_extension_handle = NULL;

    use_float = 0;
    if (readOpts) {
        bool opt =
            readOpts->GetBool("Read double variables as floats to save memory");
        if (opt)
            use_float = 1;
        else
            use_float = 0;
    }

    stack_init();

    debug1 << "avtSDFFileFormat::OpenFile() " << __LINE__ << endl;
    OpenFile(1);
}


// ****************************************************************************
//  Method: avtSDFFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

void
avtSDFFileFormat::FreeUpResources(void)
{
    debug1 << "avtSDFFileFormat::FreeUpResources(void) " << this << endl;
    sdf_free_blocklist_data(h);
}


avtSDFFileFormat::~avtSDFFileFormat(void)
{
    debug1 << "avtSDFFileFormat::~avtSDFFileFormat(void) " << this << endl;
    FreeUpResources();
    if (filename) delete [] filename;
    filename = NULL;
    if (h) sdf_close(h);
    h = NULL;
    sdf_extension_unload();
}


// ****************************************************************************
//  Method: avtSDFFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

void
avtSDFFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    debug1 << "avtSDFFileFormat::PopulateDatabaseMetaData(md:" << md
           << ") " << this << endl;

    if (md == NULL) return;

    //
    // CODE TO ADD A MESH
    //

    debug1 << "avtSDFFileFormat::OpenFile() " << __LINE__ << endl;
    OpenFile(0);

    gotMetadata = true;

    char buf[1024];
    snprintf(buf, 1024, " Job ID: %d.%d\n Code name: %s\n "
        "Code I/O version: %d\n File revision: %d\n Restart flag: %d\n "
        "Other domains: %d", h->jobid1, h->jobid2, h->code_name,
        h->code_io_version, h->file_revision, h->restart_flag,
        h->other_domains);

    md->SetDatabaseComment(buf);

    md->SetMustAlphabetizeVariables(false);

    sdf_block_t *b, *c, *next = h->blocklist;
    for (int i = 0; i < h->nblocks; i++) {
        b = h->current_block = next;
        next = b->next;
        if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH ||
                b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
            debug1 << "Found mesh: " << b->id << " " << b->name << endl;
            avtMeshType meshtype;
            int topol;
            if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
                meshtype = AVT_RECTILINEAR_MESH;
                topol = b->ndims;
            } else {
                meshtype = AVT_POINT_MESH;
                topol = 0;
            }
            avtMeshMetaData *mmd = new avtMeshMetaData(b->name, 1, 0, 0, 0,
                b->ndims, topol, meshtype);
            mmd->xUnits = b->dim_units[0];
            if (b->ndims > 1) mmd->yUnits = b->dim_units[1];
            if (b->ndims > 2) mmd->zUnits = b->dim_units[2];
            mmd->xLabel = b->dim_labels[0];
            if (b->ndims > 1) mmd->yLabel = b->dim_labels[1];
            if (b->ndims > 2) mmd->zLabel = b->dim_labels[2];
            mmd->hasSpatialExtents = true;
            for (int n=0; n < b->ndims; n++) {
                mmd->minSpatialExtents[n] = b->extents[n];
                mmd->maxSpatialExtents[n] = b->extents[b->ndims+n];
            }
            md->Add(mmd);
        } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE ||
                b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE ||
                b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED ||
                b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {

            // Now fill the metadata for a 1d or nd scalar variable
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            if (!mesh) continue;

            if (mesh->ndims == 1) {
                avtCurveMetaData *cmd = new avtCurveMetaData(b->name);
                cmd->originalName = b->name;
                cmd->validVariable = true;
                cmd->yUnits = b->units;
                cmd->xUnits = mesh->dim_units[0];
                cmd->xLabel = mesh->dim_labels[0];
                cmd->hasDataExtents = false;
                //cmd->meshName = mesh->name;
                md->Add(cmd);
                //continue;
            } else {
                avtScalarMetaData *smd = new avtScalarMetaData();
                avtCentering cent;
                // For the time being, most data is plotted as zon-centred
                // This will probably change in the future.
                if (b->stagger == SDF_STAGGER_VERTEX)
                    cent = AVT_NODECENT;
	        else
                    cent = AVT_ZONECENT;
                smd->name = b->name;
                smd->meshName = mesh->name;
                smd->centering = cent;
                smd->hasDataExtents = false;
                smd->treatAsASCII = false;
                smd->hasUnits = true;
                smd->units = b->units;
                md->Add(smd);
            }
        } else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_TENSOR
                || b->blocktype == SDF_BLOCKTYPE_MULTI_TENSOR) {
            std::string definition;
            definition.append("{");
            sdf_block_t *matvar;
            bool done = false;
            for (int i = 0; i < b->ndims; i++) {
                matvar = sdf_find_block_by_id(h, b->variable_ids[i]);
                if (!matvar) continue;
                if (done) definition.append(",");
                definition.append("<");
                definition.append(matvar->name);
                definition.append(">");
                done = true;
            }
            definition.append("}");

            Expression expr;
            expr.SetName(b->name);
            expr.SetDefinition(definition);
            expr.SetType(Expression::VectorMeshVar);
            md->AddExpression(&expr);
        } else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL
                || b->blocktype == SDF_BLOCKTYPE_MULTI_MATERIAL) {
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            if (!mesh) continue;

            avtCentering cent;
            if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                cent = AVT_ZONECENT;
            else
                cent = AVT_NODECENT;

            vector<string> mnames;
            char **matptr = b->material_names;
            for (int n = 0 ; n < b->ndims ; n++) {
                mnames.push_back(*matptr);
                matptr++;
            }

            AddMaterialToMetaData(md, b->name, mesh->name, b->ndims, mnames);
        } else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_SPECIES
                || b->blocktype == SDF_BLOCKTYPE_MULTI_SPECIES) {
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            if (!mesh) continue;
            sdf_block_t *mat = sdf_find_block_by_id(h, b->material_id);
            if (!mat) continue;

            avtCentering cent;
            if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                cent = AVT_ZONECENT;
            else
                cent = AVT_NODECENT;

            vector<string> mnames;
            char **matptr = b->material_names;
            for (int n = 0 ; n < b->ndims ; n++) {
                mnames.push_back(*matptr);
                matptr++;
            }

            vector<int> nspec;
            vector<vector<string> > specnames;
            for (int n = 0 ; n < mat->ndims ; n++) {
                if (strcmp(mat->material_names[n],b->material_name) == 0) {
                    specnames.push_back(mnames);
                    nspec.push_back(mnames.size());
                } else {
                    vector<string> tmp;
                    specnames.push_back(tmp);
                    nspec.push_back(0);
                }
            }

            AddSpeciesToMetaData(md, b->name, mesh->name, mat->name,
                mat->ndims, nspec, specnames);
        }
    }

    md->SetFormatCanDoDomainDecomposition(true);

    //
    // CODE TO ADD A VECTOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int vector_dim = 2;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddVectorVarToMetaData(md, varname, mesh_for_this_var, cent,vector_dim);
    //

    //
    // CODE TO ADD A TENSOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int tensor_dim = 9;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddTensorVarToMetaData(md, varname, mesh_for_this_var, cent,tensor_dim);
    //

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      domain      The index of the domain.  If there are NDomains, this
//                  value is guaranteed to be between 0 and NDomains-1,
//                  regardless of block origin.
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

vtkDataSet *
avtSDFFileFormat::GetMesh(int domain, const char *meshname)
{
    debug1 << "avtSDFFileFormat::GetMesh(domain:" << domain << ", meshname:"
           << meshname << ") " << this << endl;

    sdf_block_t *b = sdf_find_block_by_name(h, meshname);
    if (!b) EXCEPTION1(InvalidVariableException, meshname);
    h->current_block = b;

    debug1 << "found block:" << b->id << " for mesh:" << meshname << endl;

    ncpus = PAR_Size();
    sdf_set_ncpus(h, ncpus);

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE
            || b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        return GetCurve(domain, b);

    sdf_read_data(h);

    if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
        vtkPoints *points  = vtkPoints::New();
        vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
        points->SetNumberOfPoints(b->nlocal);
        ugrid->Allocate(b->nlocal);
        ugrid->SetPoints(points);

        if (b->datatype_out == SDF_DATATYPE_REAL4) {
            float *x = (float *)b->grids[0];
            float *y = NULL;
            float *z = NULL;
            if (b->ndims > 1) {
                y = (float *)b->grids[1];
                if (b->ndims > 2) z = (float *)b->grids[2];
            }

            vtkIdType vertex;
            float yy = 0, zz = 0;
            for (int i=0; i < b->nlocal; i++) {
                if (y) {
                    yy = y[i];
                    if (z) zz = z[i];
                }
                vertex = i;
                points->SetPoint(i, x[i], yy, zz);
                ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
            }
        } else {
            double *x = (double *)b->grids[0];
            double *y = NULL;
            double *z = NULL;
            if (b->ndims > 1) {
                y = (double *)b->grids[1];
                if (b->ndims > 2) z = (double *)b->grids[2];
            }

            vtkIdType vertex;
            double yy = 0, zz = 0;
            for (int i=0; i < b->nlocal; i++) {
                if (y) {
                    yy = y[i];
                    if (z) zz = z[i];
                }
                vertex = i;
                points->SetPoint(i, x[i], yy, zz);
                ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
            }
        }

        points->Delete();

#ifdef SDF_DEBUG
        debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
        return ugrid;
    }

    vtkDataArray *xx, *yy, *zz;

    if (b->datatype_out == SDF_DATATYPE_REAL4) {
        xx = vtkFloatArray::New();
        yy = vtkFloatArray::New();
        zz = vtkFloatArray::New();
    } else if (b->datatype_out == SDF_DATATYPE_REAL8) {
        xx = vtkDoubleArray::New();
        yy = vtkDoubleArray::New();
        zz = vtkDoubleArray::New();
    }

    xx->SetVoidArray(b->grids[0], b->local_dims[0], 1);
    yy->SetVoidArray(b->grids[1], b->local_dims[1], 1);
    zz->SetVoidArray(b->grids[2], b->local_dims[2], 1);

    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    rgrid->SetDimensions(b->local_dims);
    rgrid->SetXCoordinates(xx);
    rgrid->SetYCoordinates(yy);
    rgrid->SetZCoordinates(zz);

    xx->Delete();
    yy->Delete();
    zz->Delete();

    SetUpDomainConnectivity();

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return rgrid;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetCurve
//
//  Purpose:
//      Helper function for GetMesh. Constructs and returns a curve using
//      the 1D data read from an SDF dump.
//
//  Arguments:
//      domain      The index of the domain.  If there are NDomains, this
//                  value is guaranteed to be between 0 and NDomains-1,
//                  regardless of block origin.
//      varname     The name of the variable requested.
//
//  Programmer: Keith Bennett
//  Creation:   Tue Dec 14 2010
//
// ****************************************************************************

vtkDataSet *
avtSDFFileFormat::GetCurve(int domain, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);

    int nlocal;

    h->current_block = mesh;
    sdf_read_data(h);

    h->current_block = b;
    sdf_read_data(h);

    if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE) {
        nlocal = b->nlocal;
    } else {
        nlocal = b->dims[0];
    }

    //
    // Add all of the points to an array.
    //
    vtkPoints *pts = vtkPoints::New();
    pts->SetNumberOfPoints(nlocal);

    if (b->datatype_out == SDF_DATATYPE_REAL4) {
        float *x = (float*)mesh->grids[0];
        float *y = (float*)b->data;

        for (int i = 0 ; i < nlocal; i++)
            pts->SetPoint(i, x[i], y[i], 0.0);
    } else {
        double *x = (double*)mesh->grids[0];
        double *y = (double*)b->data;

        for (int i = 0 ; i < nlocal; i++)
            pts->SetPoint(i, x[i], y[i], 0.0);
    }

    //
    // Connect the points up with line segments.
    //
    vtkCellArray *line = vtkCellArray::New();
    for (int i = 1 ; i < nlocal; i++) {
        line->InsertNextCell(2);
        line->InsertCellPoint(i - 1);
        line->InsertCellPoint(i);
    }

    vtkPolyData *pd  = vtkPolyData::New();
    pd->SetPoints(pts);
    pd->SetLines(line);

    pts->Delete();
    line->Delete();

    return pd;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetArray
//
//  Purpose:
//      Helper function for GetVar. Fills in the blocklist data for the
//      requested variable.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

sdf_block_t *
avtSDFFileFormat::GetArray(int domain, const char *varname)
{
    debug1 << "avtSDFFileFormat::GetArray(domain:" << domain << ", varname:"
           << varname << ") " << this << endl;

    ncpus = PAR_Size();
    sdf_set_ncpus(h, ncpus);

    sdf_block_t *b = sdf_find_block_by_name(h, varname);
    if (!b) return NULL;

    debug1 << "found block:" << b->id << " for var:" << varname <<
              " type " << b->blocktype << endl;

    if (b->data) return b;
    h->current_block = b;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE ||
            b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE) {
        sdf_read_data(h);

    } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED ||
               b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {
        sdf_block_t *var;
        for (int i = 0; i < b->n_ids; i++) {
            // Fill in derived components which are not already cached
            if (b->must_read[i]) {
                var = sdf_find_block_by_id(h, b->variable_ids[i]);
                if (var && !var->data) GetArray(domain, var->name);
            }
        }

        // Allocate derived variable data if required
        if (!b->data && !b->dont_allocate) {
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            b->ndims = mesh->ndims;
            memcpy(b->local_dims, mesh->local_dims, b->ndims*sizeof(int));

            if (b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {
                b->nlocal = mesh->npoints;
            } else {
                b->nlocal = 1;
                for (int i=0; i < b->ndims; i++) {
                    if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                        b->local_dims[i]--;
                    b->nlocal *= b->local_dims[i];
                }
            }

            if (!b->datatype_out) {
                b->type_size_out = mesh->type_size_out;
                b->datatype_out = mesh->datatype_out;
            } else
                b->type_size_out = SDF_TYPE_SIZES[b->datatype_out];

            stack_alloc(b);
        }

        // Execute callback to fill in the derived variable
        if (b->populate_data) b->populate_data(h, b);
    } else
        sdf_read_data(h);

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return b;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

vtkDataArray *
avtSDFFileFormat::GetVar(int domain, const char *varname)
{
    debug1 << "avtSDFFileFormat::GetVar(domain:" << domain << ", varname:"
           << varname << ") " << this << endl;
    stack_free();

    sdf_block_t *b = GetArray(domain, varname);
    if (!b) EXCEPTION1(InvalidVariableException, varname);

    vtkDataArray *rv;
    if (b->datatype_out == SDF_DATATYPE_INTEGER4)
        rv = vtkIntArray::New();
    else if (b->datatype_out == SDF_DATATYPE_REAL4)
        rv = vtkFloatArray::New();
    else if (b->datatype_out == SDF_DATATYPE_REAL8)
        rv = vtkDoubleArray::New();

    rv->SetVoidArray(b->data, b->nlocal, 1);

    if (b->blocktype != SDF_BLOCKTYPE_POINT_MESH &&
        b->blocktype != SDF_BLOCKTYPE_POINT_VARIABLE)
            SetUpDomainConnectivity();

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return rv;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

vtkDataArray *
avtSDFFileFormat::GetVectorVar(int domain, const char *varname)
{
    debug1 << "avtSDFFileFormat::GetVectorVar(domain:" << domain << ", varname:"
           << varname << ") " << this << endl;
    // YOU MUST IMPLEMENT THIS
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a vector variable, here is some code that may be helpful.
    //
    // int ncomps = YYY;  // This is the rank of the vector - typically 2 or 3.
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // int ucomps = (ncomps == 2 ? 3 : ncomps);
    // rv->SetNumberOfComponents(ucomps);
    // rv->SetNumberOfTuples(ntuples);
    // float *one_entry = new float[ucomps];
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      int j;
    //      for (j = 0 ; j < ncomps ; j++)
    //           one_entry[j] = ...
    //      for (j = ncomps ; j < ucomps ; j++)
    //           one_entry[j] = 0.;
    //      rv->SetTuple(i, one_entry); 
    // }
    //
    // delete [] one_entry;
    // return rv;
    //
    return NULL;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetVectorVar
//
//  Purpose:
//      This method is called each time we change to a new time state.
//      Guaranteed to be called on every processor.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//   
// ****************************************************************************

void
avtSDFFileFormat::ActivateTimestep(void)
{
    debug1 << "avtSDFFileFormat::ActivateTimestep(void) " << this << endl;

#ifdef PARALLEL
    comm = VISIT_MPI_COMM;
    debug1 << "avtSDFFileFormat::A parallel" << endl;
    ncpus = PAR_Size();
    debug1 << "cpu1: " << ncpus << endl;
    MPI_Comm_size(VISIT_MPI_COMM, &ncpus);
#else
    debug1 << "avtSDFFileFormat::A serial" << endl;
#endif
    rank = PAR_Rank();
    ncpus = PAR_Size();
    debug1 << rank << " " << ncpus << " " << comm << " f:" << filename << endl;

    debug1 << "avtSDFFileFormat::OpenFile() " << __LINE__ << endl;
    OpenFile(0);
}


void
avtSDFFileFormat::SetUpDomainConnectivity(void)
{
    debug1 << "SetUpDomainConnectivity()" << endl;
    if (ncpus < 2) return;

    avtRectilinearDomainBoundaries *rdb =
        new avtRectilinearDomainBoundaries(true);

    rdb->SetNumDomains(ndomains);

    int start[3], local[3];

    for(int n=0; n < ndomains; n++) {
        int extents[6];

        sdf_get_domain_extents(h, n, start, local);

        debug1 << "SetUpDomainConnectivity0(" << n << ": ";
        for (int i = 0; i < 3; i++) {
            extents[2*i] = start[i];
            extents[2*i+1] = start[i] + local[i] - 1;
        debug1 << start[i] << "," << local[i] << " ";
        }
        debug1 << ")" << endl;

        debug1 << "SetUpDomainConnectivity(" << n << ": ";
        for (int i = 0; i < 6; i++) debug1 << extents[i] << " ";
        debug1 << endl;

        rdb->SetIndicesForRectGrid(n, extents);
    }
    rdb->CalculateBoundaries();

    void_ref_ptr vr =
        void_ref_ptr(rdb, avtStructuredDomainBoundaries::Destruct);

    cache->CacheVoidRef("any_mesh",
                       AUXILIARY_DATA_DOMAIN_BOUNDARY_INFORMATION, -1, -1, vr);
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetMaterialType
//
//  Purpose:
//      Gets an avtMaterial object for the specified domain
//
//  Arguments:
//      var        The variable of interest.
//      domain     The domain of interest.
//
//  Programmer: Keith Bennett
//  Creation:   Nov 8, 2010
//
// ****************************************************************************

template <typename Real>
void *
avtSDFFileFormat::GetMaterialType(sdf_block_t *sblock, int domain)
{
    // Read volume fraction blocks for each material
    char *var_id;
    int nm = sblock->ndims;
    sdf_block_t **vfm_blocks = new sdf_block_t *[nm];
    for (int i = 0; i < nm; i++) {
        var_id = sblock->variable_ids[i];
        vfm_blocks[i] = sdf_find_block_by_id(h, var_id);
        if (!vfm_blocks[i]) EXCEPTION1(InvalidVariableException, var_id);
        h->current_block = vfm_blocks[i];
        sdf_read_data(h);
    }
#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif

    sdf_block_t *v = vfm_blocks[0];
    int ndims = v->ndims;
    int nlocal = v->nlocal;

    int *material_list = new int[nlocal];
    int mixed_size = 0;
    int *mat_numbers = new int[nm];
    for (int n = 0; n < nm; n++) mat_numbers[n] = n + 1;

    // Fill in the pure cell array and find the size of the mixed cell arrays
    Real *vfm;
    for (int i = 0; i < nlocal; i++) {
        int material_number = 0, nmats = 0;
        Real vf;
        // Find number of materials for this cell
        for (int n = 0; n < nm; n++) {
            vfm = (Real *)vfm_blocks[n]->data;
            vf = vfm[i];
            if (vf > 0) {
                nmats++;
                material_number = mat_numbers[n];
            }
        }

        if (nmats > 1) {
            mixed_size += nmats;
            material_list[i] = -1;
        } else {
            material_list[i] = material_number;
        }
    }

    int *mix_zone = new int[mixed_size];
    int *mix_mat = new int[mixed_size];
    int *mix_next = new int[mixed_size];
    float *mix_vf = new float[mixed_size];
    int mix_index = 1;

    // Fill in the mixed cell arrays
    for (int i = 0; i < nlocal; i++) {
        // Skip pure cells
        if (material_list[i] >= 0) continue;

        material_list[i] = -mix_index;

        int material_number = 0, nmats = 0, idx = 0;
        Real vf;
        for (int n = 0; n < nm; n++) {
            vfm = (Real *)vfm_blocks[n]->data;
            vf = vfm[i];
            if (vf > 0) {
                idx = mix_index - 1;
                mix_zone[idx] = i;
                mix_mat[idx] = mat_numbers[n];
                mix_vf[idx] = vf;
                mix_next[idx] = ++mix_index;
            }
        }
        mix_next[mix_index-2] = 0;
    }

    char dom_string[128];
    sprintf(dom_string, "Domain %d", domain);

    avtMaterial *mat = new avtMaterial(nm, mat_numbers, sblock->material_names,
            ndims, v->local_dims, 0, material_list, mixed_size, mix_mat,
            mix_next, mix_zone, mix_vf, dom_string, 0); 

    delete [] mix_vf;
    delete [] mix_next;
    delete [] mix_mat;
    delete [] mix_zone;
    delete [] mat_numbers;
    delete [] material_list;
    delete [] vfm_blocks;

    debug1 << "GetMaterial() done" << endl;

    SetUpDomainConnectivity();

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return (void *)mat;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetMaterial
//
//  Purpose:
//      Gets an avtMaterial object for the specified domain
//
//  Arguments:
//      var        The variable of interest.
//      domain     The domain of interest.
//
//  Programmer: Keith Bennett
//  Creation:   Nov 8, 2010
//
// ****************************************************************************

void *
avtSDFFileFormat::GetMaterial(const char *var, int domain)
{
    debug1 << "avtSDFFileFormat::GetMaterial(var:" << var << ", domain:"
           << domain << ")" << endl;

    sdf_block_t *sblock = sdf_find_block_by_name(h, var);
    if (!sblock) EXCEPTION1(InvalidVariableException, var);
    h->current_block = sblock;

    debug1 << "found block:" << sblock->id << " for material:" << var << endl;

    if (sblock->datatype_out == SDF_DATATYPE_REAL4)
        return GetMaterialType<float>(sblock, domain);
    else
        return GetMaterialType<double>(sblock, domain);
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetSpeciesType
//
//  Purpose:
//      Gets an avtSpecies object for the specified domain
//
//  Arguments:
//      var        The variable of interest.
//      domain     The domain of interest.
//
//  Programmer: Keith Bennett
//  Creation:   Dec 7, 2010
//
// ****************************************************************************

template <typename Real>
void *
avtSDFFileFormat::GetSpeciesType(sdf_block_t *sblock, int domain)
{
    sdf_block_t *mblock = sdf_find_block_by_id(h, sblock->material_id);
    if (!mblock) EXCEPTION1(InvalidVariableException, sblock->material_id);

    //
    // we need the mixed material format computed during a GetMaterial call
    //
    void_ref_ptr vrTmp = cache->GetVoidRef(mblock->name,
            AUXILIARY_DATA_MATERIAL, -1, domain);

    if (*vrTmp == 0) {
        void *p = (void*) GetMaterial(mblock->name, domain);

        vrTmp = void_ref_ptr(p, avtMaterial::Destruct);

        cache->CacheVoidRef(mblock->name, AUXILIARY_DATA_MATERIAL, -1, domain,
            vrTmp);
    }

    avtMaterial *mat = (avtMaterial*) *vrTmp;
    if (mat == 0) EXCEPTION1(InvalidVariableException, mblock->name);

    const int* matlist  = mat->GetMatlist();
    const int* mixmat   = mat->GetMixMat();
    const int* mixnext  = mat->GetMixNext();
    int mixlen          = mat->GetMixlen();
    int *mixlist = new int[mixlen];
    for (int i=0; i<mixlen; i++) mixlist[i] = 0;

    int nmat = mat->GetNMaterials();
    std::vector<std::string> materials = mat->GetMaterials();

    int nspec = sblock->ndims;
    int specmat = 0;
    int *nspec_mats = new int[nmat];
    for (int n = 0 ; n < nmat ; n++) {
        if (strcmp(materials[n].c_str(), sblock->material_name) == 0) {
            specmat = n;
            nspec_mats[n] = nspec;
        } else {
            nspec_mats[n] = 0;
        }
    }

    // Read mass fraction blocks for each species
    char *var_id;
    sdf_block_t *vfm_block;
    Real **vfm_ptrs = new Real *[nspec];
    for (int i = 0; i < nspec; i++) {
        var_id = sblock->variable_ids[i];
        vfm_block = sdf_find_block_by_id(h, var_id);
        if (!vfm_block) EXCEPTION1(InvalidVariableException, var_id);
        h->current_block = vfm_block;
        sdf_read_data(h);
        vfm_ptrs[i] = (Real*)vfm_block->data;
    }

    int nlocal = vfm_block->nlocal;
    int ndims = vfm_block->ndims;
    int dims[ndims];
    for (int i = 0; i < ndims; i++) dims[i] = vfm_block->dims[i];

    int nmf = 0;
    int mix = 0;
    for (int i = 0; i < nlocal; i++) {
        if (matlist[i] == specmat) {
            nmf += nspec;
        } else if (matlist[i] < 0) {
            while (true) {
                if (mixmat[mix] == specmat)
                    nmf += nspec;
                if (mixnext[mix] == 0) break;
                mix = mixnext[mix] - 1;
            }
        }
    }

    float *specmf = new float[nmf * nlocal];
    int *speclist = new int[nlocal];

    nmf = 0;
    mix = 0;
    float *specptr = specmf;
    for (int i = 0; i < nlocal; i++) {
        if (matlist[i] == specmat) {
            // Chemistry material
            speclist[i] = nmf + 1;
            nmf += nspec;
            for (int m = 0; m < nspec; m++) {
                *specptr = vfm_ptrs[m][i];
                specptr++;
            }
        } else if (matlist[i] < 0) {
            // Mixed material
            speclist[i] = -(mix + 1);
            while (true) {
                if (mixmat[mix] == specmat) {
                    mixlist[mix] = nmf;
                    nmf += nspec;
                    for (int m = 0; m < nspec; m++) {
                        *specptr = vfm_ptrs[m][i];
                        specptr++;
                    }
                } else
                    mixlist[mix] = 0;
                if (mixnext[mix] == 0) break;
                mix = mixnext[mix] - 1;
            }
        } else {
            // Non-chemistry material
            speclist[i] = 0;
        }
    }

    avtSpecies *spec = new avtSpecies(nmat, nspec_mats, ndims, dims, speclist,
                                      mixlen, mixlist, nmf, specmf);

    delete [] mixlist;
    delete [] nspec_mats;
    delete [] vfm_ptrs;
    delete [] speclist;
    delete [] specmf;

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return spec;
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetSpecies
//
//  Purpose:
//      Gets an avtSpecies object for the specified domain
//
//  Arguments:
//      var        The variable of interest.
//      domain     The domain of interest.
//
//  Programmer: Keith Bennett
//  Creation:   Dec 7, 2010
//
// ****************************************************************************

void *
avtSDFFileFormat::GetSpecies(const char *var, int domain)
{
    debug1 << "avtSDFFileFormat::GetSpecies(var:" << var << ", domain:"
           << domain << ")" << endl;

    sdf_block_t *sblock = sdf_find_block_by_name(h, var);
    if (!sblock) EXCEPTION1(InvalidVariableException, var);
    h->current_block = sblock;

    debug1 << "found block:" << sblock->id << " for material:" << var << endl;

    if (sblock->datatype_out == SDF_DATATYPE_REAL4)
        return GetSpeciesType<float>(sblock, domain);
    else
        return GetSpeciesType<double>(sblock, domain);
}


// ****************************************************************************
//  Method: avtSDFFileFormat::GetAuxiliaryData
//
//  Purpose:
//      Gets the auxiliary data specified.
//
//  Arguments:
//      var        The variable of interest.
//      domain     The domain of interest.
//      type       The type of auxiliary data.
//      args       The arguments for that type -- not used.
//      df         Destructor function.
//
//  Returns:    The auxiliary data.
//
//  Programmer: Keith Bennett
//  Creation:   Nov 8, 2010
//
// ****************************************************************************

void *
avtSDFFileFormat::GetAuxiliaryData(const char *var, int domain,
         const char *type, void *args, DestructorFunction &df)
{
    debug1 << "avtSDFFileFormat::GetAuxiliaryData(var:" << var << ", domain:"
           << domain << ", type:" << type << ", ...)" << endl;
    void *rv = NULL;

    if (strcmp(type, AUXILIARY_DATA_MATERIAL) == 0) {
        rv = GetMaterial(var, domain);
        df = avtMaterial::Destruct;
    } else if (strcmp(type, AUXILIARY_DATA_SPECIES) == 0) {
        rv = GetSpecies(var, domain);
        df = avtSpecies::Destruct;
    }

#ifdef SDF_DEBUG
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return rv;
}
