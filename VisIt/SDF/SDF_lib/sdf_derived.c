#include <stdlib.h>
#include "sdf.h"


sdf_block_t *sdf_callback_grid_component(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);

    if (!b->grids) b->grids = calloc(1, sizeof(float*));
    b->data = b->grids[0] = mesh->grids[b->nm];
    b->datatype_out = mesh->datatype_out;
    if (b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED)
        b->local_dims[0] = b->nlocal = mesh->nlocal;
    else
        b->local_dims[0] = mesh->local_dims[b->nm];
    return b;
}


int sdf_add_derived_blocks(sdf_file_t *h)
{
    sdf_block_t *b, *cur, *append, *append_head, *append_tail;
    int i, len1, len2, nappend = 0;
    char *str, *name1, *name2;
    char *grid_ids[] = { "x", "y", "z" };

    cur = h->current_block;
    append = append_head = calloc(1, sizeof(sdf_block_t));

    b = h->blocklist;
    while (b) {
        if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH
            || b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
            for (i = 0; i < b->ndims; i++) {
                // Add 1d arrays for each coordinate dimension of the
                // particles. (ie. all the x's, all the y's, all the z's).
                // These are required in order to perform scatter plots.
                append->next = calloc(1, sizeof(sdf_block_t));

                nappend++;
                append = append_tail = append->next;
                append->next = NULL;

                if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
                    name1 = b->id;
                    name2 = grid_ids[i];
                    len1 = strlen(name1);
                    len2 = strlen(name2);
                    str = (char*)malloc(len1 + len2 + 2);
                    memcpy(str, name1, len1);
                    str[len1] = '/';
                    memcpy(str+len1+1, name2, len2);
                    str[len1+len2+1] = '\0';
                    append->id = str;

                    append->blocktype = SDF_BLOCKTYPE_POINT_DERIVED;
                    name2 = b->dim_labels[i];
                } else {
                    append->id = NULL;
                    SDF_SET_ENTRY_ID(append->id, grid_ids[i]);

                    append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
                    name2 = grid_ids[i];
                    append->dont_display = 1;
                }

                name1 = b->name;
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 2);
                memcpy(str, name1, len1);
                str[len1] = '/';
                memcpy(str+len1+1, name2, len2);
                str[len1+len2+1] = '\0';
                append->name = str;

                SDF_SET_ENTRY_ID(append->units, b->dim_units[i]);
                SDF_SET_ENTRY_ID(append->mesh_id, b->id);
                append->nm = i;
                append->ndims = 1;
                append->n_ids = 1;
                append->variable_ids = calloc(append->n_ids, sizeof(char*));
                SDF_SET_ENTRY_ID(append->variable_ids[0], b->id);
                append->must_read = calloc(append->n_ids, sizeof(char*));
                append->must_read[0] = 1;
                append->populate_data = sdf_callback_grid_component;
                append->done_header = 1;
                // Hack to prevent storage being allocated for this variable.
                append->dont_allocate = 1;
            }
        }
        b = b->next;
    }

    if (nappend) {
        h->tail->next = append_head->next;
        h->tail = append_tail;
        h->nblocks += nappend;
    }

    h->current_block = cur;

    return 0;
}
