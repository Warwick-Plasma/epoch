#include "BlockHandler.h"
#include "DebugStream.h"
#include "MeshReader.h"
#include "Mesh_Var_Reader.h"
#include "Stitched_Vector_Reader.h"
#include "Stitched_Magnitude_Reader.h"
#include <InvalidFilesException.h>


// ****************************************************************************


bool BlockHandler::Open(const char *filename)
{
    if (serial_in.is_open()) serial_in.close();
    this->ClearBlockChain();

    serial_in.open(filename, ios::in &ios::binary);
    if (!serial_in.good()) return false;

    char CFD[4];
    memset(CFD, 0, 4);
    serial_in.read(CFD, 3);
    debug1 << "Reading CFD Header" << endl;

    if (strcmp(CFD, "CFD") != 0) {
        debug1 << "Invalid file (No CFD header marker)" << endl;
        serial_in.close();
        return false;
    }

    serial_in.read((char*)&this->Header_Offset, sizeof(int));
    serial_in.read((char*)&this->Block_Header_Size, sizeof(int));

    serial_in.read((char*)&this->File_Version, sizeof(int));
    serial_in.read((char*)&this->File_Revision, sizeof(int));

    debug1 << this->Header_Offset << endl;

    debug1 << "Version " << File_Version << " Revision " <<
        File_Revision << endl;

    // Can't open the file with higher version number than reader
    if (File_Version < CFD_VERSION) {
        debug1 << "Invalid file (File version too high)        " << endl;
        serial_in.close();
        return false;
    }

    // However, don't care about revision number because this is a reader
    // It will include code for earlier version and just ignore additional data

    serial_in.read((char*)&this->MaxStringLen, sizeof(int));
    serial_in.read((char*)&this->nBlocks, sizeof(int));

    serial_in.read((char*)&this->Endianness, sizeof(int));
    if (Endianness == 0x0f0e0201 || Endianness == 0x02010f0e) {
        debug1 << "Wrong endianness in CFD file \"" << filename << "\"." <<endl;
        debug1 << "Please convert" << endl;
        cerr << "Wrong endianness in CFD file \"" << filename << "\"." << endl;
        cerr << "Please convert" << endl;
        EXCEPTION1(InvalidFilesException,
                   "Wrong endianness in CFD file. Please convert");
        serial_in.close();
        return false;
    }

    if (File_Revision > 0 || File_Version > 1) {
        serial_in.read((char*)&this->JobID1, sizeof(int));
        serial_in.read((char*)&this->JobID2, sizeof(int));
        serial_in.read((char*)&this->Cycle, sizeof(int));
        serial_in.read((char*)&this->Time, sizeof(double));
    } else {
        GetCycleAndTime();
    }

    if (this->nBlocks <= 0) {
        debug1 << "Invalid file (File is empty)" << endl;
    }

    debug1 << "Block contains " << this->nBlocks << " blocks " << endl;
    return true;
}

bool BlockHandler::GetCycleAndTime(void)
{
    long long offset = this->Header_Offset;
    long long type_offset = offset;
    long long block_length;
    int type;

    for (int i = 0; i < this->nBlocks; i++) {
        type_offset = offset + 2*this->MaxStringLen;
        serial_in.seekg(type_offset, ios::beg);
        serial_in.read((char*)&type, sizeof(int));
        serial_in.read((char*)&block_length, sizeof(long long));
        serial_in.read((char*)&block_length, sizeof(long long));

        if (type == TYPE_SNAPSHOT) {
            serial_in.read((char*)&this->Cycle, sizeof(int));
            serial_in.read((char*)&this->Time, sizeof(double));
            return true;
        }
        offset += this->Block_Header_Size + block_length;
    }
    return false;
}

bool BlockHandler::BuildBlockChains(void)
{
    debug1 << "Now building block chains" << endl;

    ClearBlockChain();

    long long offset = this->Header_Offset;
    for (int i = 0; i < this->nBlocks; i++) {
        Block *Current = AddBlockToChain(offset);
        if (!Current) {
            ClearBlockChain();
            serial_in.close();
            return false;
        }
        offset = offset + this->Block_Header_Size + Current->Block_Length;
        debug1 << "Offset " << MaxStringLen << endl;
        debug1 << "Read block " << Current->Name << endl;
    }

    MaxStringLen = this->MaxStringLen;
    debug1 << "Handler OK" << endl;
    return true;
}


// ****************************************************************************


void BlockHandler::ClearBlockChain()
{
    Block *Head = this->HeadBlock, *Next;
    while (Head) {
        if (Head->Name) free(Head->Name);
        if (Head->Class) free(Head->Class);
        if (Head->Reader) delete Head->Reader;
        Next = (Block*)Head->Next;
        delete Head;
        Head = Next;
    }
    this->HeadBlock = NULL;
    this->TailBlock = NULL;
}


// ****************************************************************************


Block *BlockHandler::AddBlockToChain(long long offset)
{
    Block *B = new Block();
    // Already dealing with the block header here, so no need to keep
    // How to get to it
    B->Offset = offset + this->Block_Header_Size;
    B->Name = (char*)malloc(sizeof(char) * this->MaxStringLen + 1);
    B->Class = (char*)malloc(sizeof(char) * this->MaxStringLen + 1);
    memset(B->Name, 0, this->MaxStringLen + 1);
    memset(B->Class, 0, this->MaxStringLen + 1);
    B->Reader = NULL;
    B->Prev = NULL;
    B->Next = NULL;

    serial_in.seekg(offset, ios::beg);
    serial_in.read(B->Name, this->MaxStringLen);
    serial_in.read(B->Class, this->MaxStringLen);

    serial_in.read((char*)&B->Type, sizeof(int));
    serial_in.read((char*)&B->Block_MD_Length, sizeof(long long));
    serial_in.read((char*)&B->Block_Length, sizeof(long long));

    if (B->Type == TYPE_SNAPSHOT) {
        // Now know that we're just reading in the snapshot data, so just
        // read it
        serial_in.read((char*)&this->Cycle, sizeof(int));
        serial_in.read((char*)&this->Time, sizeof(double));
       // Keeping the block is a real waste of space, but for now, what the hell
    }

    BlockReader *R = GetBlockReader(B, true);

    if (R && R->Cache()) {
        B->Reader = R;
        R->ConvertToFull();
        debug1 << "Block " << B->Name << " is caching reader" << endl;
    } else {
        delete R;
    }

    if (!this->HeadBlock) {
        // First block in list so just add it directly
        this->HeadBlock = B;
        this->TailBlock = B;
        return B;
    }

    this->TailBlock->Next = B;
    B->Prev = this->TailBlock;
    this->TailBlock = B;

    return B;
}


// ****************************************************************************


BlockReader *BlockHandler::GetBlockReader(Block *Block, bool CacheTest)
{
    debug1 << "Dealing with block " << Block->Name << " of type " <<
        Block->Type << endl;

    if (Block->Type == TYPE_MESH) {
        debug1 << "Mesh found" << endl;
        return new MeshReader(this, &serial_in, Block, MaxStringLen, CacheTest);
    } else if (Block->Type == TYPE_MESH_VARIABLE) {
        debug1 << "Mesh_Variable found" << endl;
        return new Mesh_Var_Reader(this, &serial_in, Block, MaxStringLen,
            CacheTest);
    } else if (Block->Type == TYPE_STITCHED_VECTOR) {
        debug1 << "Stitched Vector Variable found" << endl;
        return new Stitched_Vector_Reader(this, &serial_in, Block,
            MaxStringLen, CacheTest);
    } else if (Block->Type == TYPE_STITCHED_MAGNITUDE) {
        debug1 << "Stitched Magnitude Variable found" << endl;
        return new Stitched_Magnitude_Reader(this, &serial_in, Block,
            MaxStringLen, CacheTest);
    }

    return NULL;
}


// ****************************************************************************


Block *BlockHandler::GetBlockByComposite(const char *CompositeName)
{
    Block *B = this->HeadBlock;
    char *CompositeBlock = (char*)malloc(this->MaxStringLen + 1);

    while(B) {
        memset(CompositeBlock, 0, this->MaxStringLen + 1);
        GetCompositeName(B->Name, B->Class, CompositeBlock);
        if (strcmp(CompositeBlock, CompositeName) == 0) return B;
        B = (Block*)B->Next;
    }

    return NULL;
}


// ****************************************************************************


BlockReader *BlockHandler::GetReaderFromBlock(Block *Block, bool &DestroyReader)
{
    DestroyReader = false;
    if (Block->Reader) return Block->Reader;

    DestroyReader = true;
    BlockReader *R = GetBlockReader(Block, false);

    return R;
}


// ****************************************************************************


void BlockHandler::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    Block *B = this->HeadBlock;
    BlockReader *R = NULL;
    bool DestroyBlock;
    char buf[32];

    snprintf(buf, 32, "JobId: %d.%d", this->JobID1, this->JobID2);
    md->SetDatabaseComment(buf);

    debug1 << "Populating Metadata" << endl;

    while (B) {
        R = GetReaderFromBlock(B, DestroyBlock);
        debug1 << "Block Read OK " << B->Type << " " << B->Name << endl;
        if (R) {
            debug1 << "Dealing with Block " << B->Name << endl;
            // debug1 << R->Cache() << endl;
            R->PopulateDatabaseMetaData(md);
            if (DestroyBlock) delete R;
        } else {
            debug1 << "There was no Block reader for block " << B->Name <<
                " of type " << B->Type <<
                " while attempting to populate metadata" << endl;
        }

        B = (Block*)B->Next;
    }
}


// ****************************************************************************


vtkDataSet *BlockHandler::GetMesh(int domain, const char *meshname)
{
    debug1 << "Requesting mesh " << meshname << endl;
    Block *B = GetBlockByComposite(meshname);
    BlockReader *R;
    bool DestroyReader;
    vtkDataSet *Data = NULL;

    if (B) {
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetMesh(domain);
            if (DestroyReader) delete R;
        }
        debug1 << "Got mesh data from block " << B->Name << endl;
    }

    // If we're here then none of the blocks responded on their primary
    // interface. Try them on secondary interfaces

    B = this->HeadBlock;
    while (B && !Data) {
        debug1 << "Getting data from secondary interface" << endl;
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetMeshByName(domain, meshname);
            if (DestroyReader) delete R;
        }
        // If any block has responded by filling the "Data" object then
        // return it
        B = (Block*)B->Next;
    }

    return Data;
}


// ****************************************************************************


vtkDataArray *BlockHandler::GetVar(int domain, const char *varname)
{
    debug1 << "Requesting variable " << varname << endl;
    Block *B = GetBlockByComposite(varname);
    BlockReader *R;
    bool DestroyReader;
    vtkDataArray *Data = NULL;

    if (B) {
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetVar(domain);
            if (DestroyReader) delete R;
        }
    }

    // If we're here then none of the blocks responded on their primary
    // interface. Try them on secondary interfaces

    B = this->HeadBlock;
    while (B && !Data) {
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetVarByName(domain, varname);
            if (DestroyReader) delete R;
        }
        // If any block has responded by filling the "Data" object then
        // return it
        B = (Block*)B->Next;
    }

    return Data;
}


// ****************************************************************************


vtkDataArray *BlockHandler::GetVectorVar(int domain, const char *varname)
{
    Block *B = GetBlockByComposite(varname);
    BlockReader *R;
    bool DestroyReader;
    vtkDataArray *Data = NULL;

    if (B) {
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetVectorVar(domain);
            if (DestroyReader) delete R;
        }
    }

    // If we're here then none of the blocks responded on their primary
    // interface. Try them on secondary interfaces

    B = this->HeadBlock;
    while (B && !Data) {
        R = GetReaderFromBlock(B, DestroyReader);
        if (R) {
            Data = R->GetVectorVarByName(domain, varname);
            if (DestroyReader) delete R;
        }
        // If any block has responded by filling the "Data" object then
        // return it
        B = (Block*)B->Next;
    }

    return Data;
}
