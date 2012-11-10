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
#include <avtSDFOptions.h>
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
#include "commit_info.h"
#include "build_date.h"

using std::string;
using namespace SDFDBOptions;
int avtSDFFileFormat::extension_not_found = 0;

#define IJK2(i,j,k) ((i)+ng + (nx+2*ng) * ((j)+ng + (ny+2*ng) * ((k)+ng)))

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
    if (b->done_data || b->dont_own_data) return;
    b->data = calloc(b->nlocal, SDF_TYPE_SIZES[b->datatype_out]);
    memory_size += b->nlocal * SDF_TYPE_SIZES[b->datatype_out];
    stack_tail->next = tail = (struct stack*)malloc(sizeof(struct stack));
    tail->block = b;
    tail->next = NULL;
    stack_tail = tail;
}


static inline void stack_free_block(sdf_block_t *b)
{
    struct stack *old_stack_entry = stack_head;
    struct stack *stack_entry = stack_head;

    while (stack_entry) {
        if (stack_entry->block == b) {
            free(b->data);
            b->data = NULL;
            b->done_data = 0;
            memory_size -= b->nlocal * SDF_TYPE_SIZES[b->datatype_out];
            old_stack_entry->next = stack_entry->next;
            if (stack_entry == stack_tail) stack_tail = old_stack_entry;
            free(stack_entry);
            return;
        }
        old_stack_entry = stack_entry;
        stack_entry = stack_entry->next;
    }
}


static inline void stack_push_to_bottom(sdf_block_t *b)
{
    struct stack *old_stack_entry = stack_head;
    struct stack *stack_entry = stack_head;

    while (stack_entry) {
        if (stack_entry->block == b) {
            old_stack_entry->next = stack_entry->next;
            stack_tail->next = stack_entry;
            stack_tail = stack_entry;
            stack_tail->next = NULL;
            return;
        }
        old_stack_entry = stack_entry;
        stack_entry = stack_entry->next;
    }
}


static inline void stack_freeup_memory(void)
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
        memory_size -= b->nlocal * SDF_TYPE_SIZES[b->datatype_out];
        if (memory_size < MAX_MEMORY) break;
    }
}


static inline void stack_free(void)
{
    sdf_block_t *b;
    struct stack *head;

    while (stack_head->next) {
        head = stack_head;
        stack_head = stack_head->next;
        free(head);
        b = stack_head->block;
        stack_head->block = NULL;
        free(b->data);
        b->data = NULL;
        b->done_data = 0;
    }
    memory_size = 0;
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
        if (!h->rank) {
            debug1 << dlerror() << endl;
        }
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
    if (!h) h = sdf_open(filename, comm, SDF_READ, 0);
    if (!h) EXCEPTION1(InvalidFilesException, filename);
    step = h->step;
    time = h->time;
    debug1 << "avtSDFFileFormat::OpenFile h:" << h << endl;

    if (open_only) {
        sdf_close(h);
        h = NULL;
        return;
    }

    h->use_float = use_float;
    h->use_random = use_random;
    h->sdf_extension_version  = SDF_EXTENSION_VERSION;
    h->sdf_extension_revision = SDF_EXTENSION_REVISION;

    // If nblocks is negative then the file is corrupt
    if (h->nblocks <= 0) {
        int block = (-h->nblocks) / 64;
        int err = -h->nblocks - 64 * block;
        cerr << "Error code " << sdf_error_codes_c[err] << " at block "
             << block << " found for SDF file:" << endl;
        cerr << "\"" << filename << "\"" << endl;
        EXCEPTION1(InvalidFilesException, filename);
    }

    // Retrieve the extended interface library from the plugin manager
    if (!ext && !extension_not_found) ext = sdf_extension_load(h);

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
                    stack_alloc(h->current_block);
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

    // Append additional derived data for blocks added by the extension.
    sdf_add_derived_blocks_final(h);
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
    MPI_Comm_dup(VISIT_MPI_COMM, &comm);
    debug1 << "avtSDFFileFormat:: parallel" << endl;
#else
    debug1 << "avtSDFFileFormat:: serial" << endl;
#endif
    rank = PAR_Rank();
    ncpus = PAR_Size();
    ndomains = ncpus;
    debug1 << "avtSDFFileFormat:: rank:" << rank << ", ncpus:" << ncpus
           << ", comm:" << comm << ", filename:" << filename << endl;
    this->filename = new char[strlen(filename)+1];
    memcpy(this->filename, filename, strlen(filename)+1);
    gotMetadata = false;
    h = NULL;
    sdf_extension_handle = NULL;

    use_float = 0;
    use_random = 0;
    for (int i = 0; readOpts && i < readOpts->GetNumberOfOptions(); i++) {
        if (readOpts->GetName(i) == SDF_RDOPT_CONVERT_FLOAT)
            use_float = readOpts->GetBool(SDF_RDOPT_CONVERT_FLOAT) ? 1 : 0;
        else if (readOpts->GetName(i) == SDF_RDOPT_RANDOMISE)
            use_random = readOpts->GetBool(SDF_RDOPT_RANDOMISE) ? 1 : 0;
        else
            debug1 << "Ignoring unknown option \"" << readOpts->GetName(i)
                   << "\"" << endl;
    }

    stack_init();
    ext = NULL;

#ifdef MDSERVER
    debug1 << "avtSDFFileFormat::OpenFile(1) call " << __LINE__ << endl;
    OpenFile(1);
#endif
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
    stack_free();
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

    debug1 << "avtSDFFileFormat::OpenFile(0) call " << __LINE__ << endl;
    OpenFile(0);

    gotMetadata = true;

    char buf[1024];
    snprintf(buf, 1024, "\n SDF reader commit ID: %s\n "
        "SDF reader commit date: %s\n SDF reader build date: %s\n "
        "Job ID: %d.%d\n Code name: %s\n Code I/O version: %d\n "
        "File revision: %d\n Restart flag: %d\n Other domains: %d",
        SDF_READER_COMMIT_ID, SDF_READER_COMMIT_DATE, SDF_READER_BUILD_DATE,
        h->jobid1, h->jobid2, h->code_name, h->code_io_version,
        h->file_revision, h->restart_flag, h->other_domains);

    md->SetDatabaseComment(buf);

    md->SetMustAlphabetizeVariables(false);

    sdf_block_t *b, *next = h->blocklist;
    for (int i = 0; i < h->nblocks; i++) {
        b = h->current_block = next;
        next = b->next;
        if (b->dont_display) continue;

        if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH ||
                b->blocktype == SDF_BLOCKTYPE_POINT_MESH ||
                b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH) {
            debug1 << "avtSDFFileFormat:: Found mesh: id:" << b->id
                   << ", name:" << b->name << endl;
            avtMeshType meshtype;
            int topol, ndims = b->ndims;
            if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
                meshtype = AVT_RECTILINEAR_MESH;
                topol = b->ndims;
            } else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
                meshtype = AVT_POINT_MESH;
                topol = 0;
                // VisIt does not seem to handle scatter plots of 1D variables
                // properly. This hack makes changes 1D particle data to 2D
                // which works around the issue.
                if (b->ndims == 1) ndims = 2;
            } else {
                meshtype = AVT_CURVILINEAR_MESH;
                topol = b->ndims;
            }
            avtMeshMetaData *mmd = new avtMeshMetaData(b->name, 1, 0, 0, 0,
                ndims, topol, meshtype);
            mmd->xUnits = b->dim_units[0];
            if (b->ndims > 1) mmd->yUnits = b->dim_units[1];
            if (b->ndims > 2) mmd->zUnits = b->dim_units[2];
            mmd->xLabel = b->dim_labels[0];
            if (b->ndims > 1) mmd->yLabel = b->dim_labels[1];
            if (b->ndims > 2) mmd->zLabel = b->dim_labels[2];
            mmd->hasSpatialExtents = true;
            for (unsigned int n=0; n < b->ndims; n++) {
                mmd->minSpatialExtents[n] = b->extents[n];
                mmd->maxSpatialExtents[n] = b->extents[b->ndims+n];
            }
            if (b->ndims == 1 && ndims == 2) {
                mmd->yUnits = b->dim_units[0];
                mmd->yLabel = b->dim_labels[0];
                mmd->minSpatialExtents[1] = 0;
                mmd->maxSpatialExtents[1] = 0;
            }
            md->Add(mmd);
        } else if (b->blocktype == SDF_BLOCKTYPE_UNSTRUCTURED_MESH) {
            debug1 << "avtSDFFileFormat:: Found mesh: id:" << b->id
                   << ", name:" << b->name << endl;
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->subblock->mesh_id);
            avtMeshType meshtype = AVT_UNSTRUCTURED_MESH;
            int topol = b->ndims-1, ndims = b->ndims;
            avtMeshMetaData *mmd = new avtMeshMetaData(b->name, 1, 0, 0, 0,
                ndims, topol, meshtype);
            mmd->xUnits = mesh->dim_units[0];
            if (b->ndims > 1) mmd->yUnits = mesh->dim_units[1];
            if (b->ndims > 2) mmd->zUnits = mesh->dim_units[2];
            mmd->xLabel = mesh->dim_labels[0];
            if (b->ndims > 1) mmd->yLabel = mesh->dim_labels[1];
            if (b->ndims > 2) mmd->zLabel = mesh->dim_labels[2];
            mmd->hasSpatialExtents = false;
            md->Add(mmd);
        } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE ||
                b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE ||
                b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED ||
                b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {

            // Now fill the metadata for a 1d or nd scalar variable
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            if (!mesh) continue;

            if (mesh->ndims == 1 &&
                    (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE ||
                    b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED)) {
                avtCurveMetaData *cmd = new avtCurveMetaData(b->name);
                cmd->originalName = b->name;
                cmd->validVariable = true;
                cmd->yUnits = b->units;
                cmd->xUnits = mesh->dim_units[0];
                cmd->xLabel = mesh->dim_labels[0];
                cmd->hasDataExtents = false;
                // Older versions of VisIt don't support this
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
                || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_TENSOR) {
            std::string definition;
            definition.append("{");
            sdf_block_t *matvar;
            bool done = false;
            for (unsigned int i = 0; i < b->ndims; i++) {
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
                || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL) {
            sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
            if (!mesh) continue;

            avtCentering cent;
            if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                cent = AVT_ZONECENT;
            else
                cent = AVT_NODECENT;

            vector<string> mnames;
            char **matptr = b->material_names;
            for (unsigned int n = 0 ; n < b->ndims ; n++) {
                mnames.push_back(*matptr);
                matptr++;
            }

            // UGLY HACK
            // Look for an obstacle block which links with this material.
            // If found, add to the list of material names.
            if (b->subblock) {
                sdf_block_t *ob = b->subblock;
                char **matptr = ob->material_names;
                for (unsigned int n = 0 ; n < ob->ndims ; n++) {
                    mnames.push_back(*matptr);
                    matptr++;
                }
            }

            AddMaterialToMetaData(md, b->name, mesh->name, mnames.size(),
                mnames);
        } else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_SPECIES
                || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_SPECIES) {
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
            for (unsigned int n = 0 ; n < b->ndims ; n++) {
                mnames.push_back(*matptr);
                matptr++;
            }

            vector<int> nspec;
            vector<vector<string> > specnames;
            for (unsigned int n = 0 ; n < mat->ndims ; n++) {
                if (strcmp(mat->material_names[n],b->material_name) == 0
                        || strcmp(mat->id,b->material_name) == 0) {
                    specnames.push_back(mnames);
                    nspec.push_back(mnames.size());
                } else {
                    vector<string> tmp;
                    specnames.push_back(tmp);
                    nspec.push_back(0);
                }
            }

            int matdims = mat->ndims;
            // UGLY HACK
            // Look for an obstacle block which links with this material.
            // If found, add to the list of material names.
            if (mat->subblock) {
                sdf_block_t *ob = mat->subblock;
                for (unsigned int n = 0 ; n < ob->ndims ; n++) {
                    vector<string> tmp;
                    specnames.push_back(tmp);
                    nspec.push_back(0);
                }
                matdims += ob->ndims;
            }

            AddSpeciesToMetaData(md, b->name, mesh->name, mat->name,
                matdims, nspec, specnames);
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
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
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

    stack_freeup_memory();

    sdf_block_t *b = sdf_find_block_by_name(h, meshname);
    if (!b) EXCEPTION1(InvalidVariableException, meshname);
    h->current_block = b;

    debug1 << "avtSDFFileFormat:: Found block: id:" << b->id << " for mesh:"
           << meshname << endl;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE
            || b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        return GetCurve(domain, b);

    if (b->populate_data)
        b->populate_data(h, b);
    else
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
            float yy = 0, zz = 0;
            vtkIdType vertex;

            if (b->ndims > 2) {
                y = (float *)b->grids[1];
                z = (float *)b->grids[2];
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], y[i], z[i]);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            } else if (b->ndims > 1) {
                y = (float *)b->grids[1];
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], y[i], zz);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            } else {
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], yy, zz);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            }
        } else {
            double *x = (double *)b->grids[0];
            double *y = NULL;
            double *z = NULL;
            double yy = 0, zz = 0;
            vtkIdType vertex;

            if (b->ndims > 2) {
                y = (double *)b->grids[1];
                z = (double *)b->grids[2];
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], y[i], z[i]);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            } else if (b->ndims > 1) {
                y = (double *)b->grids[1];
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], y[i], zz);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            } else {
                for (int i=0; i < b->nlocal; i++) {
                    vertex = i;
                    points->SetPoint(i, x[i], yy, zz);
                    ugrid->InsertNextCell(VTK_VERTEX, 1, &vertex);
                }
            }
        }

        points->Delete();

#ifdef SDF_DEBUG
        debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
        debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
        return ugrid;

    } else if (b->blocktype == SDF_BLOCKTYPE_UNSTRUCTURED_MESH) {
        if (b->populate_data) b->populate_data(h, b);

        vtkDataArray *array = vtkFloatArray::New();
        array->SetNumberOfComponents(b->ndims);
        array->SetVoidArray(b->data, b->ndims * b->npoints, 1);

        vtkPoints *points = vtkPoints::New();
        points->SetData(array);
        array->Delete();

        vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
        ugrid->SetPoints(points);
        points->Delete();

        ugrid->Allocate(b->nfaces);

        vtkIdTypeArray *nlist = vtkIdTypeArray::New();
        nlist->SetNumberOfValues(5 * b->nfaces);
        vtkIdType *nl = nlist->GetPointer(0);
        int *node = b->node_list;

        for (int i = 0; i < b->nfaces; i++) {
            *nl++ = 4;
            *nl++ = *node++;
            *nl++ = *node++;
            *nl++ = *node++;
            *nl++ = *node++;
        }

        vtkCellArray *ca = vtkCellArray::New();
        ca->SetCells(b->nfaces, nlist);
        nlist->Delete();

        ugrid->SetCells(VTK_QUAD, ca);
        ca->Delete();

        return ugrid;
    }

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
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
        debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
        debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
        return rgrid;
    }

    if (b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH) {
        vtkPoints *points  = vtkPoints::New();
        vtkStructuredGrid *sgrid = vtkStructuredGrid::New();

        points->SetNumberOfPoints(b->nlocal);
        sgrid->SetDimensions(b->local_dims);
        sgrid->SetPoints(points);

        if (b->datatype_out == SDF_DATATYPE_REAL4) {
            float *x = (float *)b->grids[0];
            float *y = NULL;
            float *z = NULL;
            float yy = 0, zz = 0;

            if (b->ndims > 2) {
                y = (float *)b->grids[1];
                z = (float *)b->grids[2];
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], y[i], z[i]);
            } else if (b->ndims > 1) {
                y = (float *)b->grids[1];
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], y[i], zz);
            } else {
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], yy, zz);
            }
        } else {
            double *x = (double *)b->grids[0];
            double *y = NULL;
            double *z = NULL;
            double yy = 0, zz = 0;

            if (b->ndims > 2) {
                y = (double *)b->grids[1];
                z = (double *)b->grids[2];
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], y[i], z[i]);
            } else if (b->ndims > 1) {
                y = (double *)b->grids[1];
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], y[i], zz);
            } else {
                for (int i=0; i < b->nlocal; i++)
                    points->SetPoint(i, x[i], yy, zz);
            }
        }

        SetUpDomainConnectivity();

#ifdef SDF_DEBUG
        debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
        debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
        return sgrid;
    }

    return NULL;
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
    debug1 << "avtSDFFileFormat::GetCurve(domain:" << domain << ", sdf_block:"
           << b << ")" << endl;

    stack_freeup_memory();

    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);

    int nlocal;

    h->current_block = mesh;
    stack_alloc(h->current_block);
    sdf_read_data(h);

    h->current_block = b;
    stack_alloc(h->current_block);
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

    sdf_block_t *b = sdf_find_block_by_name(h, varname);
    if (!b) return NULL;

    debug1 << "avtSDFFileFormat:: Found block: id:" << b->id << " for var:"
           << varname << " type " << b->blocktype << endl;

    if (b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS) {
        sdf_block_t *var;
        if (!b->data) {
            // Free any components which have already been allocated
            for (int i = 0; i < b->n_ids; i++) {
                if (b->must_read[i]) {
                    var = sdf_find_block_by_id(h, b->variable_ids[i]);
                    if (var->data) stack_free_block(var);
                }
            }
            b->done_data = 0;
            stack_alloc(b);
        }

        for (int i = 0; i < b->n_ids; i++) {
            // Fill in derived components which are not already cached
            if (b->must_read[i]) {
                var = sdf_find_block_by_id(h, b->variable_ids[i]);
                if (var) {
                    var->data = (char*)b->data + i * var->nlocal
                            * SDF_TYPE_SIZES[var->datatype_out];
                    var->dont_own_data = 1;
                    h->current_block = var;
                    stack_alloc(h->current_block);
                    sdf_read_data(h);
                }
            }
        }
#ifdef SDF_DEBUG
        debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
        debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
        return b;
    }

    if (b->data) return b;
    h->current_block = b;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE ||
            b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE) {
        stack_alloc(h->current_block);
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
                for (unsigned int i=0; i < b->ndims; i++) {
                    if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                        b->local_dims[i]--;
                    b->nlocal *= b->local_dims[i];
                }
            }

            if (!b->datatype_out)
                b->datatype_out = mesh->datatype_out;

            stack_alloc(b);
        }

        // Execute callback to fill in the derived variable
        if (b->populate_data) b->populate_data(h, b);
    } else {
        stack_alloc(h->current_block);
        sdf_read_data(h);
    }

#ifdef SDF_DEBUG
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
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

    stack_freeup_memory();

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

#ifdef SDF_DEBUG
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
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
    MPI_Comm_dup(VISIT_MPI_COMM, &comm);
    debug1 << "avtSDFFileFormat:: parallel" << endl;
#else
    debug1 << "avtSDFFileFormat:: serial" << endl;
#endif
    rank = PAR_Rank();
    ncpus = PAR_Size();
    debug1 << "avtSDFFileFormat:: rank:" << rank << ", ncpus:" << ncpus
           << ", comm:" << comm << ", filename:" << filename << endl;

    debug1 << "avtSDFFileFormat::OpenFile(0) call " << __LINE__ << endl;
    OpenFile(0);
}


void
avtSDFFileFormat::SetUpDomainConnectivity(void)
{
    debug1 << "avtSDFFileFormat::SetUpDomainConnectivity()" << endl;
    if (ncpus < 2) return;

    avtRectilinearDomainBoundaries *rdb =
        new avtRectilinearDomainBoundaries(true);

    rdb->SetNumDomains(ndomains);

    int start[3], local[3];

    for(int n=0; n < ndomains; n++) {
        int extents[6];

        sdf_get_domain_extents(h, n, start, local);

        debug1 << "avtSDFFileFormat:: Connectivity0 (" << n << ": ";
        for (int i = 0; i < 3; i++) {
            extents[2*i] = start[i];
            extents[2*i+1] = start[i] + local[i] - 1;
        debug1 << start[i] << "," << local[i] << " ";
        }
        debug1 << ")" << endl;

        debug1 << "avtSDFFileFormat:: Connectivity1 (" << n << ": ";
        for (int i = 0; i < 6; i++) debug1 << extents[i] << " ";
        debug1 << ")" << endl;

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
        stack_alloc(h->current_block);
        sdf_read_data(h);
    }
#ifdef SDF_DEBUG
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif

    sdf_block_t *v = vfm_blocks[0];
    int ndims = v->ndims;
    int nlocal = v->nlocal;

    int *material_list = new int[nlocal];
    int mixed_size = 0;
    int nmat = nm;
    char **mat_names = sblock->material_names;

    int *obdata = NULL;
    int nx, ny, ng = 0;
    if (sblock->subblock) {
        sdf_block_t *obgrp = sblock->subblock;
        sdf_block_t *ob = obgrp->subblock;
        h->current_block = ob;
        stack_alloc(h->current_block);
        sdf_read_data(h);
        obdata = (int *)ob->data;
        ng = ob->ng;
        nx = ob->local_dims[0] - 2 * ng;
        ny = ob->local_dims[1] - 2 * ng;
        nmat += obgrp->ndims;
        char **mat_namesptr, **matptr;
        mat_names = new char*[nmat];
        mat_namesptr = matptr = mat_names = sblock->material_names;
        for (int n = 0; n < nm; n++) {
            *mat_namesptr = *matptr;
            mat_namesptr++;
            matptr++;
        }
        matptr = obgrp->material_names;
        for (int n = 0; n < obgrp->ndims; n++) {
            *mat_namesptr = *matptr;
            mat_namesptr++;
            matptr++;
        }
    }

    int *mat_numbers = new int[nmat];
    for (int n = 0; n < nmat; n++) mat_numbers[n] = n + 1;

    // Fill in the pure cell array and find the size of the mixed cell arrays
    Real *vfm;
    int material_number = 0, nmats = 0;
    Real vf;
    for (int i = 0; i < nlocal; i++) {
        if (obdata) {
            if (ng) {
                // Skip ghostcells in obstacle data if they exist.
                int rem = i / nx;
                int kk = rem / ny;
                int ii = i - nx * rem;
                int jj = rem - ny * kk;
                if (ndims < 3) kk = -ng;
                int obgrp = obdata[IJK2(ii,jj,kk)];
                if (obgrp) {
                    material_list[i] = obgrp + nm;
                    continue;
                }
            } else {
                if (*obdata) {
                    material_list[i] = *obdata + nm;
                    obdata++;
                    continue;
                }
                obdata++;
            }
        }

        material_number = 0;
        nmats = 0;
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

        int idx = 0;
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

    avtMaterial *mat = new avtMaterial(nmat, mat_numbers,
            mat_names, ndims, v->local_dims, 0, material_list, mixed_size,
            mix_mat, mix_next, mix_zone, mix_vf, dom_string, 0);

    delete [] mix_vf;
    delete [] mix_next;
    delete [] mix_mat;
    delete [] mix_zone;
    delete [] mat_numbers;
    delete [] material_list;
    delete [] vfm_blocks;
    if (obdata) delete [] mat_names;

    debug1 << "avtSDFFileFormat::GetMaterial() done" << endl;

#ifdef SDF_DEBUG
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
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

    stack_freeup_memory();

    sdf_block_t *sblock = sdf_find_block_by_name(h, var);
    if (!sblock) EXCEPTION1(InvalidVariableException, var);
    h->current_block = sblock;

    debug1 << "avtSDFFileFormat:: Found block: id:" << sblock->id
           << " for material:" << var << endl;

    if (sblock->datatype_out == SDF_DATATYPE_OTHER) {
        if (sblock->variable_ids && sblock->variable_ids[0]) {
            sdf_block_t *b = sdf_find_block_by_id(h, sblock->variable_ids[0]);
            sblock->datatype_out = b->datatype_out;
        } else {
            EXCEPTION1(InvalidVariableException, sblock->id);
            return NULL;
        }
    }

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
        stack_alloc(h->current_block);
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
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
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

    stack_freeup_memory();

    sdf_block_t *sblock = sdf_find_block_by_name(h, var);
    if (!sblock) EXCEPTION1(InvalidVariableException, var);
    h->current_block = sblock;

    debug1 << "avtSDFFileFormat:: Found block: id:" << sblock->id
           << " for material:" << var << endl;

    if (sblock->datatype_out == SDF_DATATYPE_OTHER) {
        if (sblock->variable_ids && sblock->variable_ids[0]) {
            sdf_block_t *b = sdf_find_block_by_id(h, sblock->variable_ids[0]);
            sblock->datatype_out = b->datatype_out;
        } else {
            EXCEPTION1(InvalidVariableException, sblock->id);
            return NULL;
        }
    }

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
    debug1 << "avtSDFFileFormat:: SDF debug buffer: ";
    debug1 << h->dbg_buf; h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    return rv;
}
