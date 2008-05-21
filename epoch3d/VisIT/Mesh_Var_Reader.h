#ifndef MESH_VAR_READER_H
#define MESH_VAR_READER_H
#include "BlockReader.h"
#include "MeshReader.h"

class Mesh_Var_Reader : public BlockReader
{
 private:

    int VarType;
    int Dimensions;
    int SizeOfFloat;

    int *Dims;
    void *Stagger;
    double Extents[2];

    long long nPart;

    char *MeshName;
    char *MeshClass;
    Block *MeshBlock;
    
 public:
 virtual vtkDataArray * GetVectorVar(int domain){return NULL;}
 virtual vtkDataArray * GetVar(int domain);
 virtual vtkDataSet * GetMesh(int domain);
 virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *md);
 virtual bool Cache(){return true;}

 virtual InternalMetaData* GetInternalMetaData();
 virtual void DestroyInternalMetaData(InternalMetaData *md) {delete md;}

    Mesh_Var_Reader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen,bool CacheOnly);
    ~Mesh_Var_Reader();
};

class Mesh_Variable_MetaData  : public InternalMetaData
{
 public:

    int *Dims;
    void *Stagger;
    double *Extents;

    long long n_Elements;
    char *MeshName, *MeshClass;

    int SizeOfFloat;
    int Dimensions;
    int VarType;
};

#endif
