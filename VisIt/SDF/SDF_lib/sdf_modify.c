#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sdf.h>
#ifdef PARALLEL
#include <mpi.h>
#endif



int sdf_modify_array(sdf_file_t *h, sdf_block_t *b, void *data)
{
    if (!b) return 1;

    return sdf_modify_array_section(h, b, data, NULL, NULL);
}



int sdf_modify_array_section(sdf_file_t *h, sdf_block_t *b, void *data,
        uint64_t *starts, uint64_t *ends)
{
    uint64_t *data_starts, *data_ends;
    int i, j, k, alloc_starts = 0, alloc_ends = 0, combined_writes = 0;
    int loop_counts[3] = {1, 1, 1};
    uint64_t offset_starts[3] = {0}, offset_ends[3] = {0};
    uint64_t length, count, offset, data_offset;
    int errcode = 1, sz;

    // Sanity checks
    if (!b || b->data_location < h->first_block_location
            || b->data_length <= 0) return errcode;

    assert(sizeof(*starts) == sizeof(*b->dims));
    assert(b->ndims < 5);

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_PLAIN_VARIABLE:
    case SDF_BLOCKTYPE_POINT_VARIABLE:
    //case SDF_BLOCKTYPE_ARRAY:
    //case SDF_BLOCKTYPE_PLAIN_DERIVED:
    //case SDF_BLOCKTYPE_POINT_DERIVED:
        break;
    default:
        return errcode;
    }
    if (!b->dims) return errcode;

    if (starts) {
        data_starts = starts;
    } else {
        data_starts = calloc(b->ndims, sizeof(*data_starts));
        alloc_starts = 1;
    }

    if (ends) {
        data_ends = ends;
    } else {
        data_ends = malloc(b->ndims * sizeof(*data_ends));
        alloc_ends = 1;
        memcpy(data_ends, b->dims, b->ndims * sizeof(*data_ends));
    }

    // Sanity check
    // data_starts[i] < 0 not required since it is unsigned
    for (i = 0; i < b->ndims; i++)
        if (data_ends[i] > b->dims[i] ||
                data_starts[i] > data_ends[i]) goto cleanup;

    // First check for any writes which can be combined
    sz = length = SDF_TYPE_SIZES[b->datatype];
    for (i = 0; i < b->ndims; i++) {
        count = data_ends[i] - data_starts[i];
        length *= count;
        if (count == b->dims[i]) {
            combined_writes++;
            sz *= b->dims[i];
        } else
            break;
    }

    for (i = combined_writes, j = 0; i < b->ndims; i++, j++) {
        offset_starts[j] = data_starts[i] * sz;
        offset_ends[j] = (b->dims[i] - data_ends[i]) * sz;
        sz *= b->dims[i];
    }

    for (i = combined_writes+1, j = 0; i < b->ndims; i++, j++)
        loop_counts[j] = data_ends[i] - data_starts[i];

    // Dumb implementation, but I can't think of anything clever right now
    // Loop over a maximum of 3-dimensions writing array chunks
    // (allows a 4D array to be written)

    offset = data_offset = 0;
    for (k = 0; k < loop_counts[2]; k++) {
        offset += offset_starts[2];
        for (j = 0; j < loop_counts[1]; j++) {
            offset += offset_starts[1];
            for (i = 0; i < loop_counts[0]; i++) {
                offset += offset_starts[0];
                memcpy((char*)b->data + offset,
                       (char*)data + data_offset, length);
                // Only need to update file if this block is present
                if (b->block_start > 0)
                    sdf_write_at(h, b->data_location+offset,
                                 (char*)data+data_offset, length);
                data_offset += length;
                offset += length + offset_ends[0];
            }
            offset += offset_ends[1];
        }
        offset += offset_ends[2];
    }

    errcode = 0;
cleanup:

    if (alloc_starts) free(data_starts);
    if (alloc_ends) free(data_ends);

    return errcode;
}



int sdf_modify_array_element(sdf_file_t *h, sdf_block_t *b, void *data,
        uint64_t *index)
{
    uint64_t index1[4] = {0};
    int i;

    if (!b) return 1;

    for (i = 0; i < b->ndims; i++) index1[i] = index[i] + 1;

    return sdf_modify_array_section(h, b, data, index, index1);
}
