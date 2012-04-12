#include <stdlib.h>
#include "sdf.h"


sdf_block_t *sdf_callback_grid_component(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
    b->data = mesh->grids[b->nm];
    return b;
}


int sdf_add_derived_blocks(sdf_file_t *h)
{
    sdf_block_t *b, *cur, *last_block, *append, *append_head;
    int i, len1, len2, nappend = 0;
    char *str;

    cur = h->current_block;
    append = append_head = calloc(1, sizeof(sdf_block_t));

    last_block = b = h->blocklist;
    while (b) {
        if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
            for (i = 0; i < b->ndims; i++) {
                append->next = calloc(1, sizeof(sdf_block_t));

                nappend++;
                append = append->next;
                append->next = NULL;

                len1 = strlen(b->id);
                len2 = strlen(b->dim_labels[i]);
                str = (char*)malloc(len1 + len2 + 2);
                memcpy(str, b->id, len1);
                str[len1] = '/';
                memcpy(str+len1+1, b->dim_labels[i], len2);
                str[len1+len2+1] = '\0';
                append->id = str;

                len1 = strlen(b->name);
                str = (char*)malloc(len1 + len2 + 2);
                memcpy(str, b->name, len1);
                str[len1] = '/';
                memcpy(str+len1+1, b->dim_labels[i], len2);
                str[len1+len2+1] = '\0';
                append->name = str;

                SDF_SET_ENTRY_ID(append->units, b->dim_units[i]);
                SDF_SET_ENTRY_ID(append->mesh_id, b->id);
                append->nm = i;
                append->ndims = 1;
                append->variable_ids = calloc(append->ndims, sizeof(char*));
                SDF_SET_ENTRY_ID(append->variable_ids[0], b->id);
                append->populate_data = sdf_callback_grid_component;
                append->done_header = 1;
                append->blocktype = SDF_BLOCKTYPE_POINT_DERIVED;
            }
        }
        last_block = b;
        b = b->next;
    }

    if (nappend) {
        last_block->next = append_head->next;
        h->nblocks += nappend;
    }

    h->current_block = cur;

    return 0;
}
