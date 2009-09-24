#ifndef STITCHED_VECTOR_READER_H
#define STITCHED_VECTOR_READER_H
#include "BlockReader.h"
#include "Mesh_Var_Reader.h"

class Stitched_Vector_Reader : public BlockReader
{
 private:

    int Dimensions;
    int SizeOfFloat;

    long long n_Elements;

    char *MeshName;
    char *MeshClass;
    Block **SubBlocks;
    
 public:
 virtual vtkDataArray * GetVectorVar(int domain);
 virtual vtkDataArray * GetVar(int domain){return NULL;}
 virtual vtkDataSet * GetMesh(int domain){return NULL;}
 virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *md);
 virtual bool Cache(){return true;}

    Stitched_Vector_Reader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen,bool CacheOnly);
    ~Stitched_Vector_Reader();
};
#endif
