#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <sdf.h>
#ifdef PARALLEL
#include <mpi.h>
#endif



static int sdf_read_inline_block_locations(sdf_file_t *h)
{
    sdf_block_t *b;
    uint64_t next_location;

    h->use_summary = 0;
    h->tmp_flag = 1; // set tmp_flag to force reread
    sdf_read_summary(h);
    h->tmp_flag = 0;

    b = h->blocklist;
    h->start_location = h->current_location = 0;
    next_location = h->first_block_location;
    while (b) {
        b->inline_block_start = next_location;
        SDF_READ_ENTRY_INT8(next_location);
        b->inline_next_block_location = next_location;
        h->current_location += h->block_header_length + b->info_length - 8;
        b = b->next;
    }
    free(h->buffer);
    h->buffer = NULL;

    h->start_location = h->summary_location;
    h->inline_metadata_read = 1;

    return 0;
}



static int sdf_read_summary_block_locations(sdf_file_t *h)
{
    sdf_block_t *b;
    uint64_t next_location;

    h->use_summary = 1;
    h->tmp_flag = 1; // set tmp_flag to force reread
    sdf_read_summary(h);
    h->tmp_flag = 0;

    b = h->blocklist;
    h->start_location = h->current_location = h->summary_location;
    next_location = h->summary_location;
    while (b) {
        b->inline_block_start = b->block_start;
        b->inline_next_block_location = b->next_block_location;
        b->block_start = next_location;
        SDF_READ_ENTRY_INT8(next_location);
        h->current_location += h->block_header_length + b->info_length - 8;
        b = b->next;
    }
    free(h->buffer);
    h->buffer = NULL;

    h->start_location = 0;
    h->summary_metadata_read = 1;

    return 0;
}



static void sdf_modify_rewrite_header(sdf_file_t *h)
{
    off_t offset, zero = 0;

    if (h->rank != h->rank_master) return;

    offset = 4 * SOI4 + h->id_length;

    if (h->blocklist)
        sdf_write_at(h, offset, &h->blocklist->inline_block_start, SOI8);
    else
        sdf_write_at(h, offset, &zero, SOI8);
    offset += SOI8;

    sdf_write_at(h, offset, &h->summary_location, SOI8); offset += SOI8;

    sdf_write_at(h, offset, &h->summary_size, SOI4); offset += SOI4;

    sdf_write_at(h, offset, &h->nblocks_file, SOI4); offset += SOI4;
}



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
                if (b->in_file)
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



static int modify_remove_block(sdf_file_t *h, sdf_block_t *block, int freeit)
{
    if (!block) return 0;

    if (block->in_file) {
        if (!h->inline_metadata_read && !h->inline_metadata_invalid)
            sdf_read_inline_block_locations(h);

        if (!h->summary_metadata_read && !h->summary_metadata_invalid)
            sdf_read_summary_block_locations(h);

        if (block->prev) {
            block->prev->next_block_location = block->next_block_location;
            block->prev->next_block_modified = 1;
        } else {
            h->first_block_modified = 1;
        }
        h->nblocks_file--;
    }

    if (block->next) {
        if (block->prev) {
            block->prev->next = block->next;
            block->next->prev = block->prev;
        } else {
            block->next->prev = NULL;
            h->blocklist = block->next;
        }
    } else {
        if (block->prev) {
            block->prev->next = NULL;
            h->tail = block->prev;
        } else
            h->blocklist = NULL;
    }

    if (block == h->last_block_in_file) {
        h->can_truncate = 1;
        if (block->prev) {
            h->last_block_in_file = block->prev;
            block->prev->next_block_location = 0;
        }
    }

    if (freeit) sdf_free_block(h, block);
    h->metadata_modified = 1;
    h->nblocks--;

    return 0;
}



int sdf_modify_remove_block(sdf_file_t *h, sdf_block_t *block)
{
    return modify_remove_block(h, block, 1);
}



int sdf_modify_remove_block_id(sdf_file_t *h, const char *id)
{
    sdf_block_t *b = sdf_find_block_by_id(h, id);

    if (b)
        return sdf_modify_remove_block(h, b);
    else
        return -1;
}



int sdf_modify_remove_block_name(sdf_file_t *h, const char *name)
{
    sdf_block_t *b = sdf_find_block_by_name(h, name);

    if (b)
        return sdf_modify_remove_block(h, b);
    else
        return -1;
}



int sdf_modify_rewrite_metadata(sdf_file_t *h)
{
    sdf_block_t *b, *next;
    int use_summary;

    if (!h->metadata_modified) return 0;
    if (h->rank != h->rank_master) return 0;

    next = b = h->blocklist;

    use_summary = h->use_summary;
    h->use_summary = 0;
    if (!h->first_block_modified) {
        // Find first modified block
        while (next) {
            b = next;
            next = b->next;
            if (!b->in_file) continue;
            if (b->next_block_modified) break;
            if (b->rewrite_metadata) {
                h->current_block = b;
                b->done_header = b->done_info = 0;
                sdf_write_meta(h);
                b->rewrite_metadata = 0;
            }
        }
    }
    next = b;

    if (!h->blocklist) {
        h->summary_size = 0;
        h->summary_location = h->first_block_location;
    }

    // Now rewrite remaining inline block locations
    while (next) {
        b = next;
        next = b->next;
        if (!b->in_file) continue;
        if (b->rewrite_metadata) {
            h->current_block = b;
            b->done_header = b->done_info = 0;
            sdf_write_meta(h);
            b->rewrite_metadata = 0;
        } else {
            if (b->next) {
                sdf_write_at(h, b->inline_block_start,
                        &b->next->inline_block_start, SOI8);
            } else {
                sdf_write_at(h, b->inline_block_start,
                        &h->summary_location, SOI8);
            }
        }
    }
    h->use_summary = use_summary;

    sdf_write_new_summary(h);

    sdf_modify_rewrite_header(h);

    h->metadata_modified = h->first_block_modified = 0;

    sdf_flush(h);

    return 0;
}
