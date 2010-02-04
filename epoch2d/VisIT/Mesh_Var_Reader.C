#include "Mesh_Var_Reader.h"
#include "constants.h"
#include "CommonCommands.h"
#include <vtkDoubleArray.h>
#include <DebugStream.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>


// ****************************************************************************


Mesh_Var_Reader::Mesh_Var_Reader(BlockHandler *Handler, ifstream *file,
    Block *Owner, int MaxStringLen, bool CacheOnly)
     : BlockReader(Handler, file, Owner, MaxStringLen, CacheOnly)
{
    // This is always a caching reader, so ignore CacheOnly
    // Note that despite this, it only caches the metadata, not the primary data

    // Can't guarantee that we'll be in the right place, so go there
    file->seekg(this->Owner->Offset, ios::beg);
    file->read((char*)&this->VarType, sizeof(int));
    file->read((char*)&this->Dimensions, sizeof(int));
    file->read((char*)&this->SizeOfFloat, sizeof(int));

    this->MeshName = (char*)malloc(this->MaxStringLen + 1);
    this->MeshClass = (char*)malloc(this->MaxStringLen + 1);

    memset(this->MeshName, 0, this->MaxStringLen + 1);
    memset(this->MeshClass, 0, this->MaxStringLen + 1);

    if (this->VarType == VAR_CARTESIAN) {
        this->Dims = (int*)malloc(sizeof(int) * this->Dimensions);
        file->read((char*)this->Dims, sizeof(int) * this->Dimensions);

        // It should be possible to do something about the grid stagger using
        // the VisIt Centring commands, but I don't know how for the moment
        this->Stagger = malloc(this->SizeOfFloat * this->Dimensions);
        file->read((char*)this->Stagger, this->SizeOfFloat * this->Dimensions);

        file->read((char*)this->Extents, this->SizeOfFloat * 2);

        file->read(this->MeshName, this->MaxStringLen);
        file->read(this->MeshClass, this->MaxStringLen);
    } else if (this->VarType == VAR_PARTICLE) {
        this->Dims = NULL;
        this->Stagger = NULL;
        file->read((char*)&this->nPart, sizeof(long long));
        file->read((char*)this->Extents, this->SizeOfFloat * 2);

        file->read(this->MeshName, this->MaxStringLen);
        file->read(this->MeshClass, this->MaxStringLen);
    }

    char *MeshComposite = (char*)malloc(2 * MaxStringLen + 1);
    memset(MeshComposite, 0, 2 * MaxStringLen + 1);
    GetCompositeName(this->MeshName, this->MeshClass, MeshComposite);
    this->MeshBlock = Handler->GetBlockByComposite(MeshComposite);
    free(MeshComposite);
}


// ****************************************************************************


Mesh_Var_Reader::~Mesh_Var_Reader()
{
    free(this->Dims);
    free(this->Stagger);
    free(this->MeshName);
    free(this->MeshClass);
}


// ****************************************************************************


void Mesh_Var_Reader::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    char *Composite = (char*)malloc(2 * MaxStringLen + 1);
    memset(Composite, 0, 2 * MaxStringLen + 1);
    GetCompositeName(this->Owner->Name, this->Owner->Class, Composite);

    char *MeshComposite = (char*)malloc(2 * MaxStringLen + 1);
    memset(MeshComposite, 0, 2 * MaxStringLen + 1);
    GetCompositeName(this->MeshName, this->MeshClass, MeshComposite);

    debug1 << "Mesh associated with Block " << Owner->Name <<
        " is " << MeshName << endl;

    avtCentering cent = AVT_NODECENT;

    bool DestroyReader;
    MeshReader *R =
        (MeshReader*)this->Handler->GetReaderFromBlock(this->MeshBlock,
            DestroyReader);

    if (R) {
        Mesh_MetaData *MD = (Mesh_MetaData*)R->GetInternalMetaData();

        if (MD->MeshType == MESH_CARTESIAN && MD->Dimensions == 1) {
            avtCurveMetaData *cmd = new avtCurveMetaData(Composite);
            md->Add(cmd);
        } else {
            avtScalarMetaData *smd = new avtScalarMetaData(Composite,
                MeshComposite, cent);
            md->Add(smd);
        }
        if (DestroyReader) delete R;
    }

    free(Composite);
    free(MeshComposite);
}


// ****************************************************************************


vtkDataArray *Mesh_Var_Reader::GetVar(int domain)
{
    // Spin past metadata
    file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length, ios::beg);

    // Don't support >3D data at the moment, will have to do some slicing
    // when we do
    if (this->Dimensions > 3) return NULL;
    debug1 << "Ndims OK at " << this->Dimensions << endl;
    vtkDataArray *Data;

    if (this->SizeOfFloat == 4)
        Data = vtkFloatArray::New();
    else
        Data = vtkDoubleArray::New();

    if (this->VarType == VAR_CARTESIAN) {
        if (this->Dimensions > 1) {
            long long varlen_total = 1;
            for (int i = 0; i < this->Dimensions; i++)
                varlen_total = varlen_total * this->Dims[i];

            debug1 << "Cartesian mesh with " << varlen_total << " points found";

            Data->SetNumberOfTuples(varlen_total);
            void *DatVoid = Data->GetVoidPointer(0);
            file->read((char*)DatVoid, varlen_total * this->SizeOfFloat);
        }
        // Do nothing for 1D cartesian. They are dealt with in other ways.
    } else if (this->VarType == VAR_PARTICLE) {
        Data->SetNumberOfTuples(this->nPart);
        void *DatVoid = Data->GetVoidPointer(0);
        file->read((char*)DatVoid, this->nPart * this->SizeOfFloat);
    }

    return Data;
}


// ****************************************************************************


vtkDataSet *Mesh_Var_Reader::GetMesh(int domain)
{
    bool DestroyReader;
    MeshReader *R =
        (MeshReader*)this->Handler->GetReaderFromBlock(this->MeshBlock,
            DestroyReader);

    if (!R) return NULL;

    void *v1 = R->GetCartAxis0(domain);
    void *v2 = malloc(this->SizeOfFloat * this->Dims[0]);
    if (DestroyReader) delete R;
    file->seekg(this->Owner->Offset + this->Owner->Block_MD_Length, ios::beg);
    file->read((char*)v2, this->Dims[0] * this->SizeOfFloat);

    vtkDataArray *Data;
    if (this->SizeOfFloat == 4)
        Data = vtkFloatArray::New();
    else
        Data = vtkDoubleArray::New();

    Data->SetNumberOfComponents(3);
    Data->SetNumberOfTuples(this->Dims[0]);

    vtkPolyData *pd  = vtkPolyData::New();
    vtkPoints   *pts = vtkPoints::New();
    pd->SetPoints(pts);

    if (this->SizeOfFloat == 4) {
        for (int i = 0 ; i < this->Dims[0]; i++)
            Data->SetTuple3(i, ((float*)v1)[i], ((float*)v2)[i], 0.0);
    } else {
        for (int i = 0 ; i < this->Dims[0]; i++)
            Data->SetTuple3(i, ((double*)v1)[i], ((double*)v2)[i], 0.0);
    }

    pts->SetData(Data);
    Data->Delete();
    free(v1);
    free(v2);

    vtkCellArray *line = vtkCellArray::New();
    pd->SetLines(line);
    for (int i = 1 ; i < this->Dims[0]; i++) {
        line->InsertNextCell(2);
        line->InsertCellPoint(i - 1);
        line->InsertCellPoint(i);
    }

    pts->Delete();
    line->Delete();

    return pd;
}


// ****************************************************************************


InternalMetaData *Mesh_Var_Reader::GetInternalMetaData()
{
    debug1 << this->Owner->Name << " preparing metadata" << endl;
    Mesh_Variable_MetaData *MD = new Mesh_Variable_MetaData();
    MD->Dims = this->Dims;

    MD->Stagger = this->Stagger;
    MD->Extents = this->Extents;

    if (this->VarType == VAR_CARTESIAN) {
        MD->n_Elements = 1;
        for (int i = 0; i < this->Dimensions; i++)
            MD->n_Elements *= this->Dims[i];
    } else if (this->VarType == VAR_PARTICLE)
        MD->n_Elements = this->nPart;

    MD->MeshName = this->MeshName;
    MD->MeshClass = this->MeshClass;

    MD->SizeOfFloat = this->SizeOfFloat;
    MD->Dimensions = this->Dimensions;
    MD->VarType = this->VarType;

    return MD;
}
