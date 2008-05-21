#ifndef MESHREADERH
#define MESHREADERH
#include "BlockReader.h"

class MeshReader : public BlockReader
{
 private:

    int MeshType;
    int Dimensions;
    int SizeOfFloat;

    int Part_Coord_Type;
    long long nPart;

    int *Dims;
    void *Extents;
    
 public:
 virtual vtkDataArray * GetVectorVar(int domain){return NULL;}
 virtual vtkDataArray * GetVar(int domain){return NULL;}
 virtual vtkDataSet * GetMesh(int domain);
 virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *md);
 virtual vtkDataArray * GetVarByName(int domain,const char* varname);
 virtual void* GetCartAxis0(int domain);

 virtual bool Cache(){return true;}

    MeshReader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen,bool CacheOnly);
    ~MeshReader();

    InternalMetaData* GetInternalMetaData();
    void DestroyInternalMetaData(InternalMetaData *MD) {delete MD;}
};

class Mesh_MetaData  : public InternalMetaData
{
 public:

    int *Dims;
    void *Extents;

    long long n_Elements;

    int SizeOfFloat;
    int Dimensions;
    int MeshType;
    int Part_Coord_Type;
};
#endif
