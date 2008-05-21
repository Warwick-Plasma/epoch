#include "BlockReader.h"
BlockReader::BlockReader(BlockHandler *Handler,ifstream *file,Block* Owner,int MaxStringLen, bool CacheOnly)
{
    this->Handler=Handler;
    this->file=file;
    this->Owner=Owner;
    this->CacheOnly=CacheOnly;
    this->MaxStringLen=MaxStringLen;
}
