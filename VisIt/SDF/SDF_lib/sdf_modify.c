#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <sdf.h>
#include "sdf_output.h"
#include "sdf_input.h"
#include "sdf_control.h"
#ifdef PARALLEL
#include <mpi.h>
#endif

/**
 * @defgroup modify
 * @brief Routines for modifying an SDF file in-place
 */

#define COPY_ENTRY(copy, original, count) do { \
        size_t _len; \
        if ((copy)) { \
            _len = (count) * sizeof(*(copy)); \
            (copy) = malloc(_len); \
            memcpy((copy), (original), _len); \
        } \
    } while(0)

#define COPY_ENTRY_NAME(copyname, count) \
   COPY_ENTRY((copy->copyname), (original->copyname), (count))

#define COPY_ENTRY_STR(copyob, originalob) do { \
        size_t _len; \
        if ((copyob)) { \
            _len = strlen((copyob)) * sizeof(char) + 1; \
            (copyob) = malloc(_len); \
            memcpy((copyob), (originalob), _len); \
        } \
    } while(0)

#define COPY_ENTRY_STRING(copyname) \
    COPY_ENTRY_STR(copy->copyname, original->copyname)

#define COPY_ENTRY_STRING_ARRAY(copyname, count) do { \
        size_t _len, _i; \
        char **_ptr; \
        if (copy->copyname) { \
            _len = (count) * sizeof(char*); \
            _ptr = copy->copyname = malloc(_len); \
            for (_i = 0; _i < (count); _i++) {\
                COPY_ENTRY_STR(*_ptr, original->copyname); \
                _ptr++; \
            } \
        } \
    } while(0)



static int sdf_read_inline_block_locations(sdf_file_t *h)
{
    sdf_block_t *b;
    int64_t next_location;

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
    int64_t next_location;

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



static int copy_block(sdf_block_t *copy, const sdf_block_t *original)
{
    const sdf_block_t *b = original;

    memcpy(copy, original, sizeof(*original));

    COPY_ENTRY_NAME(extents, 2*b->ndims);
    if (b->blocktype == SDF_BLOCKTYPE_STATION) {
        COPY_ENTRY_NAME(dim_mults, b->nvariables);
        COPY_ENTRY_STRING_ARRAY(dim_units, b->nvariables);
        COPY_ENTRY_STRING_ARRAY(variable_ids, b->nvariables);
        COPY_ENTRY_STRING_ARRAY(material_names, b->nvariables);
    } else {
        COPY_ENTRY_NAME(dim_mults, b->ndims);
        COPY_ENTRY_STRING_ARRAY(dim_units, b->ndims);
        COPY_ENTRY_STRING_ARRAY(variable_ids, b->ndims);
        COPY_ENTRY_STRING_ARRAY(material_names, b->ndims);
    }
    COPY_ENTRY_NAME(station_x, b->nstations);
    COPY_ENTRY_NAME(station_y, b->nstations);
    COPY_ENTRY_NAME(station_z, b->nstations);
    COPY_ENTRY_NAME(station_nvars, b->nstations);
    COPY_ENTRY_NAME(variable_types, b->nvariables);
    COPY_ENTRY_NAME(station_move, b->nstations);
    COPY_ENTRY_STRING(id);
    COPY_ENTRY_STRING(units);
    COPY_ENTRY_STRING(mesh_id);
    COPY_ENTRY_STRING(material_id);
    COPY_ENTRY_STRING(vfm_id);
    COPY_ENTRY_STRING(obstacle_id);
    COPY_ENTRY_STRING(name);
    COPY_ENTRY_STRING(material_name);
    COPY_ENTRY_STRING_ARRAY(dim_labels, b->ndims);
    COPY_ENTRY_STRING_ARRAY(station_ids, b->nstations);
    COPY_ENTRY_STRING_ARRAY(station_names, b->nstations);
    copy->nelements_blocks = NULL;
    copy->data_length_blocks = NULL;
    copy->dims_in = NULL;
    copy->station_index = NULL;
    copy->station_id = NULL;
    copy->must_read = NULL;
    copy->node_list = NULL;
    copy->boundary_cells = NULL;
    copy->grids = NULL;
    copy->data = NULL;
    copy->next = NULL;
    copy->prev = NULL;
    copy->subblock = NULL;
    copy->subblock2 = NULL;
    copy->block_start = copy->next_block_location = copy->data_location
            = copy->inline_block_start = copy->inline_next_block_location = 0;
    copy->in_file = 0;

    return 0;
}



static sdf_block_t *append_block_to_blocklist(sdf_file_t *h, sdf_block_t *b)
{
    h->tail->next = b;
    b->prev = h->tail;
    b->next = NULL;
    h->tail = b;
    h->nblocks++;

    return b;
}



/**
 * @ingroup modify
 * @{
 */
int sdf_modify_array(sdf_file_t *h, sdf_block_t *b, void *data)
{
    if (!b) return 1;

    return sdf_modify_array_section(h, b, data, NULL, NULL);
}



// Should be updated to mimic Python's array slice notation
int sdf_modify_array_section(sdf_file_t *h, sdf_block_t *b, void *data,
        int64_t *starts, int64_t *ends)
{
    int64_t *data_starts, *data_ends;
    int i, n, alloc_starts = 0, alloc_ends = 0, combined_writes = 0;
    int *loop_counts = NULL, *idx = NULL;
    int64_t *offset_starts = NULL, *offset_ends = NULL;
    int64_t length, count, offset, nwrites, data_location;
    int errcode = 1, sz, ndims;
    char *copy_data_ptr, *block_data_ptr;

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

    nwrites = 1;
    ndims = b->ndims - combined_writes;
    if (ndims > 0) {
        offset_starts = malloc(ndims * sizeof(*offset_starts));
        offset_ends   = malloc(ndims * sizeof(*offset_ends));
        loop_counts   = malloc(ndims * sizeof(*loop_counts));
        idx           = malloc(ndims * sizeof(*idx));

        n = combined_writes;
        for (i = 0; i < ndims; i++) {
            offset_starts[i] = data_starts[n] * sz;
            offset_ends[i] = (b->dims[n] - data_ends[n]) * sz;
            sz *= b->dims[n];
            n++;
            if (n < ndims)
                loop_counts[i] = data_ends[n] - data_starts[n];
            else
                loop_counts[i] = 1;
            nwrites *= loop_counts[i];
            idx[i] = 0;
        }
    }

    // Loop over an arbitrary number of dimensions writing array chunks

    offset = 0;
    copy_data_ptr = data;
    block_data_ptr = b->data;
    data_location = b->data_location;

    for (n = 0; n < nwrites; n++) {
        for (i = 0; i < ndims; i++) {
            offset += offset_starts[i];
            if (idx[i] != 0)
                break;
        }

        block_data_ptr += offset;
        memcpy(block_data_ptr, copy_data_ptr, length);

        // Only need to update file if this block is present
        if (b->in_file) {
            data_location += offset;
            sdf_write_at(h, data_location, copy_data_ptr, length);
        }

        copy_data_ptr += length;
        offset = length;

        for (i = 0; i < ndims; i++) {
            offset += offset_ends[i];
            idx[i]++;
            if (idx[i] != loop_counts[i]) break;
            idx[i] = 0;
        }
    }

    errcode = 0;
cleanup:

    if (alloc_starts)  free(data_starts);
    if (alloc_ends)    free(data_ends);

    if (offset_starts) free(offset_starts);
    if (offset_ends)   free(offset_ends);
    if (loop_counts)   free(loop_counts);
    if (idx)           free(idx);

    return errcode;
}



int sdf_modify_array_element(sdf_file_t *h, sdf_block_t *b, void *data,
        int64_t *index)
{
    int64_t index1[4] = {0};
    int i;

    if (!b) return 1;

    for (i = 0; i < b->ndims; i++) index1[i] = index[i] + 1;

    return sdf_modify_array_section(h, b, data, index, index1);
}



int sdf_modify_add_block(sdf_file_t *h, sdf_block_t *block)
{
    sdf_block_t *b, *next;
    int64_t extent = h->first_block_location;
    int64_t sz;

    append_block_to_blocklist(h, block);

    block->in_file = 0;
    block->rewrite_metadata = h->metadata_modified = 1;

    if (!h->inline_metadata_read && !h->inline_metadata_invalid)
        sdf_read_inline_block_locations(h);

    if (!h->summary_metadata_read && !h->summary_metadata_invalid)
        sdf_read_summary_block_locations(h);

    // First find the furthest extent of the inline metadata and/or data
    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;
        if (!b->in_file) continue;

        sz = b->inline_block_start + h->block_header_length + b->info_length;
        if (sz > extent) extent = sz;

        sz = b->data_location + b->data_length;
        if (sz > extent) extent = sz;
    }

    b = h->last_block_in_file;
    block->block_start = b->next_block_location =
            b->block_start + h->block_header_length + b->info_length;
    block->inline_block_start = b->inline_next_block_location = extent;
    block->data_location = extent + h->block_header_length +
            block->info_length;
    h->summary_location = block->inline_next_block_location =
            block->data_location + block->data_length;
    block->rewrite_metadata = 1;
    block->in_file = 1;
    h->last_block_in_file = block;

    h->metadata_modified = 1;
    h->nblocks_file++;

    return 0;
}



int sdf_modify_add_block_copy(sdf_file_t *h, sdf_block_t *copy)
{
    sdf_block_t *new = malloc(sizeof(sdf_block_t));

    copy_block(new, copy);

    return sdf_modify_add_block(h, new);
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



int sdf_modify_add_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material)
{
    sdf_block_t *b = stitched;
    char **new_material_names, **new_variable_ids;
    int64_t info_length, ndims;
    int i, errcode = 1;

    if (!stitched || !material) return errcode;

    // Check that the blocks are valid for this operation
    switch (stitched->blocktype) {
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL:
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATVAR:
        break;
    default:
        return errcode;
    }

    if (material->blocktype != SDF_BLOCKTYPE_PLAIN_VARIABLE) return errcode;

    sdf_modify_add_block_copy(h, material);

    ndims = b->ndims;
    b->ndims++;

    new_variable_ids   = malloc(b->ndims * sizeof(char*));

    for (i=0; i < ndims; i++) {
        new_variable_ids[i] = strdup(b->variable_ids[i]);
        free(b->variable_ids[i]);
    }

    new_variable_ids[i] = strdup(material->id);
    free(b->variable_ids);
    b->variable_ids = new_variable_ids;

    info_length = SOI4 + (b->ndims + 2) * h->id_length;

    if (stitched->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL
            || stitched->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL) {
        new_material_names = malloc((b->ndims+1) * sizeof(char*));

        for (i=0; i < ndims; i++) {
            new_material_names[i] = strdup(b->material_names[i]);
            free(b->material_names[i]);
        }

        new_material_names[i] = strdup(material->name);
        free(b->material_names);
        b->material_names = new_material_names;
        b->blocktype = SDF_BLOCKTYPE_STITCHED_MATERIAL;

        info_length = SOI4 + (b->ndims + 1) * h->id_length
            + b->ndims * h->string_length;
    } else
        b->blocktype = SDF_BLOCKTYPE_STITCHED_MATVAR;

    // Allow for the possibility that info_length was padded
    if (info_length > b->info_length) {
        b->info_length = info_length;
        modify_remove_block(h, b, 0);
        sdf_modify_add_block(h, b);
    }

    return 0;
}



int sdf_modify_remove_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material)
{
    int i, n, len1, len2, matnum = -1;
    char **names, **ids;
    int errcode = 1;

    if (!stitched || !material) return errcode;

    // Check that the blocks are valid for this operation
    switch (stitched->blocktype) {
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL:
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATVAR:
        break;
    default:
        return errcode;
    }

    if (material->blocktype != SDF_BLOCKTYPE_PLAIN_VARIABLE) return errcode;

    // Find the material to be removed
    len1 = strlen(material->id) + 1;
    for (i = 0; i < stitched->ndims; i++) {
        len2 = strlen(stitched->variable_ids[i]) + 1;
        if (len1 == len2) {
            if (!memcmp(material->id, stitched->variable_ids[i], len1)) {
                matnum = i;
                break;
            }
        }
    }

    if (matnum < 0) return errcode;

    ids = malloc((stitched->ndims-1) * sizeof(char*));

    for (i = 0, n = 0; i < stitched->ndims; i++) {
        if (i == matnum)
            free(stitched->variable_ids[i]);
        else
            ids[n++] = stitched->variable_ids[i];
    }

    free(stitched->variable_ids);
    stitched->variable_ids = ids;

    stitched->blocktype = SDF_BLOCKTYPE_STITCHED_MATVAR;

    if (stitched->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL
            || stitched->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL) {
        names = malloc((stitched->ndims-1) * sizeof(char*));

        for (i = 0, n = 0; i < stitched->ndims; i++) {
            if (i == matnum)
                free(stitched->material_names[i]);
            else
                ids[n++] = stitched->material_names[i];
        }

        free(stitched->material_names);
        stitched->material_names = names;

        stitched->blocktype = SDF_BLOCKTYPE_STITCHED_MATERIAL;
    }

    stitched->ndims--;
    stitched->rewrite_metadata = 1;

    sdf_modify_remove_block(h, material);

    return 0;
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
/** @} */
