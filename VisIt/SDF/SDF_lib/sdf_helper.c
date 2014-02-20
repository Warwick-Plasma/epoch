#include <string.h>
#include <sdf.h>
#include "sdf_helper.h"
#include "stack_allocator.h"


int sdf_helper_read_data(sdf_file_t *h, sdf_block_t *b)
{
    int i;
    sdf_block_t *block;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED ||
            b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {
        for (i = 0; i < b->n_ids; i++) {
            // Fill in derived components which are not already cached
            if (b->must_read[i]) {
                block = sdf_find_block_by_id(h, b->variable_ids[i]);
                if (block && !block->data) sdf_helper_read_data(h, block);
            }
        }

        // Allocate derived variable data if required
        if (!b->data && !b->dont_allocate) {
            block = sdf_find_block_by_id(h, b->mesh_id);
            b->ndims = block->ndims;
            memcpy(b->local_dims, block->local_dims, b->ndims*sizeof(int));

            if (b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED) {
                b->nelements_local = block->dims[0];
            } else {
                b->nelements_local = 1;
                for (i = 0; i < b->ndims; i++) {
                    if (b->stagger == SDF_STAGGER_CELL_CENTRE)
                        b->local_dims[i]--;
                    b->nelements_local *= b->local_dims[i];
                }
            }

            if (!b->datatype_out)
                b->datatype_out = block->datatype_out;

            stack_alloc(b);
        }

        // Execute callback to fill in the derived variable
        if (b->populate_data) b->populate_data(h, b);
    }

    stack_alloc(b);

    h->current_block = b;
    return sdf_read_data(h);
}
