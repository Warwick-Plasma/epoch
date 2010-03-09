#ifndef CFD_BLOCKREADERH
#define CFD_BLOCKREADERH

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "Block.h"
#include "BlockHandler.h"

#ifdef PARALLEL
//    #include <mpi.h>
#endif

#include <string>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <Expression.h>

#include <InvalidVariableException.h>
#include <avtSTMDFileFormat.h>

class InternalMetaData
{
public:
    int MDType;
    int CreatorType;
};

// There will be a blockhandler class so define it's existance
class BlockHandler;

class BlockReader
{
protected:
    ifstream *file;
    Block *Owner;
    bool CacheOnly;
    int MaxStringLen;
    BlockHandler *Handler;

public:
    virtual vtkDataArray *GetVectorVar(int domain) = 0;
    virtual vtkDataArray *GetVar(int domain) = 0;
    virtual vtkDataSet *GetMesh(int domain) = 0;
    virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *md) = 0;

    // These are the secondary interfaces to the reader
    // Only implement these if a single reader is supposed to register multiple
    // objects. This should only happen fairly rarely. Never implment them for
    // non-caching readers. That would be very slow

    virtual vtkDataArray *GetVectorVarByName(int domain, const char *varname)
        { return NULL; }
    virtual vtkDataArray *GetVarByName(int domain, const char *varname)
        { return NULL; }
    virtual vtkDataSet *GetMeshByName(int domain, const char *meshname)
        { return NULL; }

    virtual void FreeUpResources() {;}
    virtual bool Cache() { return false; }
    virtual void ConvertToFull() {;}

    virtual InternalMetaData *GetInternalMetaData() { return NULL; }
    virtual void DestroyInternalMetaData(InternalMetaData *MD) {;}

    virtual ~BlockReader() {;}
    BlockReader(BlockHandler *Handler, ifstream *file, Block *Owner,
        int MaxStringLength, bool CacheOnly);
};

#endif
