#ifndef CFDBLOCKHANDLERH
#define CFDBLOCKHANDLERH

#include "constants.h"
#include "CommonCommands.h"
#include "Block.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

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
#include "BlockReader.h"

class BlockHandler
{
    friend class BlockReader;

protected:

    Block *HeadBlock, *TailBlock;
    ifstream serial_in;
    int Header_Offset;
    int Block_Header_Size;
    int MaxStringLen;
    int nBlocks;
    int File_Version, File_Revision;
    int Endianness;
    int JobID1;
    int JobID2;

    // These are strictly data from blocks, but put them here since they're
    // Used a lot and can only meaningfully be present once

    int Cycle;
    double Time;

    void ClearBlockChain();
    Block *AddBlockToChain(long long offset);

public:

    bool Open(const char *filename);
    void *GetAuxiliaryData(const char *var, int domain, const char *type,
        void *, DestructorFunction &df) { df = NULL; return NULL; }
    int GetCycle(void) const { return this->Cycle; }
    double GetTime(void) { return this->Time; }
    bool GetCycleAndTime(void);

    vtkDataArray *GetVectorVar(int domain, const char *varname);
    vtkDataArray *GetVar(int domain, const char *varname);
    vtkDataSet *GetMesh(int domain, const char *meshname);
    void PopulateDatabaseMetaData(avtDatabaseMetaData *md);
    void FreeUpResources() {;}
    bool BuildBlockChains(void);

    BlockReader *GetBlockReader(Block *B, bool CacheTest);

    long long GetOffsetByNameandClass(const char *Name, const char *Class)
        { return -1; }
    long long GetOffsetByComposite(const char *Name) { return -1; }
    Block *GetBlockByComposite(const char *Name);
    BlockReader *GetReaderFromBlock(Block *Block, bool &DestroyReader);

    BlockHandler() { this->HeadBlock=NULL; }
    ~BlockHandler() { this->ClearBlockChain(); }
};

#endif
