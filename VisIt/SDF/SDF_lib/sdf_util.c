#include <string.h>
#include <sdf.h>
#include "sdf_util.h"


sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id)
{
    sdf_block_t *current, *b;
    size_t len;
    int i;

    if (!h || !h->blocklist || !id)
        return NULL;

    current = h->blocklist;
    len = strlen(id) + 1;
    for (i=0; i < h->nblocks; i++) {
        b = current;
        if (memcmp(id, b->id, len) == 0) return b;
        current = b->next;
    }

    return NULL;
}


sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name)
{
    sdf_block_t *current, *b;
    size_t len;
    int i;

    if (!h || !h->blocklist || !name)
        return NULL;

    current = h->blocklist;
    len = strlen(name) + 1;
    for (i=0; i < h->nblocks; i++) {
        b = current;
        if (memcmp(name, b->name, len) == 0) return b;
        current = b->next;
    }

    return NULL;
}
