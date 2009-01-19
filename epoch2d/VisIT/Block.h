#ifndef BLOCKH
#define BLOCKH

class BlockReader;

typedef struct Block_struct
{
//The name and class of the block
char *Name,*Class;
//The type of the block
int Type;
//The offset in the file of the block
long long Offset, Block_Length,Block_MD_Length;
//The next and previous entries in the linked list
void *Next, *Prev;

BlockReader *Reader;
} Block;
#endif
