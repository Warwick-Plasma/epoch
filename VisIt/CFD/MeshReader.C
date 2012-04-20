#include "MeshReader.h"
#include "constants.h"
#include "CommonCommands.h"
#include <vtkDoubleArray.h>
#include <DebugStream.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <avtDatabaseMetaData.h>

#define min(a, b) ((a) < (b) ? (a) : (b))


// ****************************************************************************


MeshReader::MeshReader(BlockHandler *Handler, ifstream *file, Block *Owner,
    int MaxStringLen, bool CacheOnly)
     : BlockReader(Handler, file, Owner, MaxStringLen, CacheOnly)
{
    // This is a caching reader, so ignore CacheOnly
    // Note that despite this, it only caches the metadata, not the primary data

    // Can't guarantee that we'll be in the right place, so go there
    file->seekg(this->Owner->Offset, ios::beg);
    file->read((char*)&this->MeshType, sizeof(int));
    file->read((char*)&this->Dimensions, sizeof(int));
    file->read((char*)&this->SizeOfFloat, sizeof(int));

    if (this->MeshType == MESH_CARTESIAN) {
        this->Dims = (int*)malloc(sizeof(int) * this->Dimensions);
        file->read((char*)this->Dims, sizeof(int) * this->Dimensions);

        this->Extents = malloc(this->SizeOfFloat * this->Dimensions);
        file->read((char*)this->Extents, this->SizeOfFloat * this->Dimensions);
    } else {
        file->read((char*)&this->Part_Coord_Type, sizeof(int));
        this->Dims = NULL;
        file->read((char*)&this->nPart, sizeof(long long));
        this->Extents = malloc(this->SizeOfFloat * this->Dimensions);
        file->read((char*)this->Extents, this->SizeOfFloat * this->Dimensions);
    }
}


// ****************************************************************************


MeshReader::~MeshReader()
{
    if (this->Dims) free(this->Dims);
    if (this->Extents) free(this->Extents);
}


// ****************************************************************************


void MeshReader::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    debug1 << "Populating metadata" << endl;
    char *Composite = (char*)malloc(2 * MaxStringLen + 1);
    std::string Composite_for_copy;

    memset(Composite, 0, 2 * MaxStringLen + 1);

    debug1 << "Setting composite name" << endl;
    GetCompositeName(this->Owner->Name, this->Owner->Class, Composite);

    debug1 << "Creating mesh metadata" << endl;

    if (this->MeshType == MESH_CARTESIAN && this->Dimensions == 1) {
        // 1D Cartesian meshes are meaningless (almost). They are replaced
        // by avtCurve objects in VisIt. This means that there is no mesh
        // metadata
        return;
    } else {
        int spatdim, topdim;
        avtMeshType mt;
        spatdim = this->Dimensions;

        if (this->MeshType == MESH_CARTESIAN) {
            mt = AVT_RECTILINEAR_MESH;
            topdim = spatdim;
        } else {
            mt = AVT_POINT_MESH;
            topdim = 0;
        }

        avtMeshMetaData *mmd =
            new avtMeshMetaData(Composite, 1, 0, 0, 0, spatdim, topdim, mt);
        debug1 << "Composite name is " << Composite << endl;
        debug1 << "Name length reserved" << endl;

        // std::string s = Composite;
        // mmd->name.assign(s);

        md->Add(mmd);
        debug1 << "Metadata added" << endl;

        if (this->MeshType == MESH_PARTICLE) {
            char axes_labels[6] = "X\0Y\0Z";
            char *VarComposite = (char*)malloc(3 * MaxStringLen + 1);
            avtScalarMetaData *smd = NULL;

            for (int i = 0; i < this->Dimensions; i++) {
                memset(VarComposite, 0, 3 * MaxStringLen + 1);
                GetCompositeName(axes_labels + i * 2, Composite, VarComposite);

                avtCentering cent = AVT_NODECENT;

                smd = new avtScalarMetaData(VarComposite, Composite, cent);
                md->Add(smd);
            }
            free(VarComposite);
        }
    }
    debug1 << "About to free composite name" << endl;
    free(Composite);
    debug1 << "Mesh " << Composite << " added OK" << endl;
}


// ****************************************************************************


vtkDataSet *MeshReader::GetMesh(int domain)
{
    // Spin past metadata
    file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length, ios::beg);

    debug1 << "Meta data block length for Mesh " << this->Owner->Name <<
        " is " << this->Owner->Block_MD_Length << endl;

    if (this->MeshType == MESH_CARTESIAN) {
        if (this->Dimensions > 3) return NULL;
        int Visit_Dims[3] = {1, 1, 1};
        memcpy(Visit_Dims, this->Dims, sizeof(int) * this->Dimensions);
        vtkDataArray *coords[3] = {0, 0, 0};

        for(int i = 0; i < 3; i++) {
            if (this->SizeOfFloat == 4)
                coords[i] = vtkFloatArray::New();
            else
                coords[i] = vtkDoubleArray::New();

            coords[i]->SetNumberOfTuples(Visit_Dims[i]);

            if (i < this->Dimensions) {
                void *v = coords[i]->GetVoidPointer(0);
                file->read((char*)v, this->SizeOfFloat * Visit_Dims[i]);
            } else {
                coords[i]->SetComponent(0, 0, 0.);
            }
        }

        debug1 << "Setting up grid" << endl;
        vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();

        debug1 << "Dimensions of grid are " << Visit_Dims[0] << " " <<
            Visit_Dims[1] << " " << Visit_Dims[2] << endl;

        rgrid->SetDimensions(Visit_Dims);
        rgrid->SetXCoordinates(coords[0]);
        rgrid->SetYCoordinates(coords[1]);
        rgrid->SetZCoordinates(coords[2]);
        coords[1]->Delete();
        coords[2]->Delete();
        coords[0]->Delete();

        debug1 << "About to return grid" << endl;
        return rgrid;
    } else {
        if (this->Dimensions > 3) return NULL;
        vtkPoints *points = vtkPoints::New();
        points->SetNumberOfPoints(this->nPart);
        float *pointdata = (float*)points->GetVoidPointer(0);
        float *pointnull = pointdata;

        // Visit uses a slightly mad way of representing points, so buffer
        // them in 100000 particles shouldn't tax any system too
        // much = (100000*8 ~= 800K RAM). You may want to increase this number
        // by a factor of 10-50 for parallel file systems this routine really
        // needs fixing anyway, since it probably (CHECK THIS) requires twice
        // as much RAM as is strictly should due to the presence of the points
        // and the unstructured grid objects at the same time
        long long npart_section = 100000;
        void *v = malloc(this->SizeOfFloat * npart_section);

        for (int dimswing = 0; dimswing < this->Dimensions; dimswing++) {
            long long npart_left = this->nPart;
            float temp;

            while (npart_left >0) {
                int imax = min(npart_section, npart_left);
                file->read((char*)v, this->SizeOfFloat * imax);
                for (int i = 0; i < imax; i++) {
                    if (this->SizeOfFloat == 4) {
                        // Just copy the memory over if compatible
                        temp = *((float*)v + i);
                    } else {
                        double d = *((double*)v + i);
                        temp = d;
                    }
                    // Finally copy the data into the points array
                    *pointdata = temp;
                    // Move pointer
                    pointdata += 3;
                }
                // Reduce remaining particle count
                npart_left = npart_left - npart_section;
            }
            pointdata = pointnull + dimswing + 1;
        }
        free(v);

        vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
        ugrid->SetPoints(points);
        points->Delete();
        ugrid->Allocate(this->nPart);
        vtkIdType onevertex;

        for (int i = 0; i < this->nPart; i++) {
            onevertex = i;
            ugrid->InsertNextCell(VTK_VERTEX, 1, &onevertex);
        }
        return ugrid;
    }
}


// ****************************************************************************


void *MeshReader::GetCartAxis0(int domain)
{
    file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length, ios::beg);
    void *v = malloc(this->SizeOfFloat * this->Dims[0]);
    file->read((char*)v, this->SizeOfFloat * this->Dims[0]);
    return v;
}


// ****************************************************************************


vtkDataArray *MeshReader::GetVarByName(int domain, const char *varname)
{
    debug1 << "Requesting " << varname << " by seconday interface" << endl;
    if (this->MeshType != MESH_PARTICLE) return NULL;

    char *Composite = (char*)malloc(2 * MaxStringLen + 1);
    memset(Composite, 0, 2 * MaxStringLen + 1);
    GetCompositeName(this->Owner->Name, this->Owner->Class, Composite);

    char axes_labels[6] = "X\0Y\0Z";
    char *VarComposite = (char*)malloc(3 * MaxStringLen + 1);
    vtkDataArray *Data = NULL;

    for (int i = 0; i < this->Dimensions; i++) {
        memset(VarComposite, 0, 3 * MaxStringLen + 1);
        GetCompositeName(axes_labels + i * 2, Composite, VarComposite);

        if (strcmp(varname, VarComposite) == 0) {
            file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length +
                this->nPart * this->SizeOfFloat * i, ios::beg);
            if (this->SizeOfFloat == 4) Data = vtkFloatArray::New();
            if (this->SizeOfFloat == 8) Data = vtkDoubleArray::New();
            if (!Data) return NULL;
            Data->SetNumberOfTuples(this->nPart);
            void *DataVoid = Data->GetVoidPointer(0);
            file->read((char*)DataVoid, this->nPart * this->SizeOfFloat);
            break;
        }
    }
    free(VarComposite);
    free(Composite);
    return Data;
}


// ****************************************************************************


InternalMetaData *MeshReader::GetInternalMetaData()
{
    debug1 << this->Owner->Name << " preparing metadata" << endl;
    Mesh_MetaData *MD = new Mesh_MetaData();
    MD->Dims = this->Dims;
    MD->Extents = this->Extents;

    if (this->MeshType == MESH_CARTESIAN) {
        MD->n_Elements = 1;
        for (int i = 0; i < this->Dimensions; i++) {
            MD->n_Elements *= this->Dims[i];
        }
    } else {
        MD->n_Elements = this->nPart;
    }

    MD->SizeOfFloat = this->SizeOfFloat;
    MD->Dimensions = this->Dimensions;
    MD->MeshType = this->MeshType;

    return MD;
}
