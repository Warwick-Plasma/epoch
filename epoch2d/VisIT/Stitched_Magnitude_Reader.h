#ifndef STITCHED_MAGNITUDE_READER_H
#define STITCHED_MAGNITUDE_READER_H
#include "BlockReader.h"
#include "Mesh_Var_Reader.h"

class Stitched_Magnitude_Reader : public BlockReader
{
 private:

    int Dimensions;
    int SizeOfFloat;

    long long n_Elements;

    char *MeshName;
    char *MeshClass;
    Block **SubBlocks;
    
 public:
 virtual vtkDataArray * GetVar(int domain);
 virtual vtkDataArray * GetVectorVar(int domain){return NULL;}
 virtual vtkDataSet * GetMesh(int domain){return NULL;}
 virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *md);
 virtual bool Cache(){return true;}

    Stitched_Magnitude_Reader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen,bool CacheOnly);
    ~Stitched_Magnitude_Reader();
};
#endif
