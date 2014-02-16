#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sdf.h>
#ifdef PARALLEL
#include <mpi.h>
#endif


#define SDF_COMMON_INFO() do { \
    if (!h->current_block || !h->current_block->done_header) { \
      if (h->rank == h->rank_master) { \
        fprintf(stderr, "*** ERROR ***\n"); \
        fprintf(stderr, "SDF block header has not been read." \
                " Ignoring call.\n"); \
      } \
      return 1; \
    } \
    b = h->current_block; \
    if (b->done_info) return 0; \
    h->current_location = b->block_start + h->block_header_length; \
    b->done_info = 1; } while(0)


static int sdf_read_constant(sdf_file_t *h);
static int sdf_read_stitched(sdf_file_t *h);
static int sdf_read_stitched_material(sdf_file_t *h);
static int sdf_read_stitched_matvar(sdf_file_t *h);
static int sdf_read_stitched_species(sdf_file_t *h);
static int sdf_read_stitched_obstacle_group(sdf_file_t *h);
static int sdf_read_array(sdf_file_t *h);
static int sdf_read_array_info(sdf_file_t *h);
static int sdf_read_cpu_split_info(sdf_file_t *h);
static int sdf_read_run_info(sdf_file_t *h);
static void build_summary_buffer(sdf_file_t *h);
static int sdf_read_next_block_header(sdf_file_t *h);



void sdf_trim(char *str)
{
    int i, len = strlen(str);
    char *ptr = str + len - 1;

    for (i=0, ptr=str+len-1; i<len && *ptr==' '; i++, ptr--)
        *ptr = '\0';
}



static inline int sdf_get_next_block(sdf_file_t *h)
{
    if (h->blocklist) {
        if (!h->current_block)
            h->current_block = h->blocklist;
        else if (h->current_block->next)
            h->current_block = h->current_block->next;
        else {
            sdf_block_t *block = malloc(sizeof(sdf_block_t));
            memset(block, 0, sizeof(sdf_block_t));
            if (h->use_summary)
                block->block_start = h->tail->next_block_location;
            else
                block->block_start = h->current_location;
            h->tail->next = block;
            h->tail->next->prev = h->tail;
            h->current_block = h->tail = block;
        }
    } else {
        sdf_block_t *block = malloc(sizeof(sdf_block_t));
        memset(block, 0, sizeof(sdf_block_t));
        if (h->use_summary)
            block->block_start = h->summary_location;
        else
            block->block_start = 0;
        h->blocklist = h->tail = h->current_block = block;
    }

    return 0;
}



int sdf_read_bytes(sdf_file_t *h, char *buf, int buflen)
{
#ifdef PARALLEL
    return MPI_File_read(h->filehandle, buf, buflen, MPI_BYTE,
            MPI_STATUS_IGNORE);
#else
    return (1 != fread(buf, buflen, 1, h->filehandle));
#endif
}



int sdf_read_header(sdf_file_t *h)
{
    int buflen;

    h->indent = 0;

    if (h->done_header) return 1;

    buflen = h->first_block_location;
    h->buffer = malloc(buflen);

    h->current_location = h->start_location = 0;

    if (h->rank == h->rank_master) {
        sdf_seek(h);
        sdf_read_bytes(h, h->buffer, buflen);
    }
    sdf_broadcast(h, h->buffer, buflen);

    // If this isn't SDF_MAGIC then this isn't a SDF file;
    if (memcmp(h->buffer, SDF_MAGIC, 4) != 0) {
        sdf_close(h);
        return -1;
    }

    h->current_location += 4;

    SDF_READ_ENTRY_INT4(h->endianness);
    if (h->endianness == 0x0f0e0201) h->swap = 1;

    SDF_READ_ENTRY_INT4(h->file_version);
    if (h->file_version > SDF_VERSION) {
        sdf_close(h);
        return -1;
    }

    SDF_READ_ENTRY_INT4(h->file_revision);

    SDF_READ_ENTRY_ID(h->code_name);

    SDF_READ_ENTRY_INT8(h->first_block_location);

    SDF_READ_ENTRY_INT8(h->summary_location);

    SDF_READ_ENTRY_INT4(h->summary_size);

    SDF_READ_ENTRY_INT4(h->nblocks);

    SDF_READ_ENTRY_INT4(h->block_header_length);

    SDF_READ_ENTRY_INT4(h->step);

    SDF_READ_ENTRY_REAL8(h->time);

    SDF_READ_ENTRY_INT4(h->jobid1);

    SDF_READ_ENTRY_INT4(h->jobid2);

    SDF_READ_ENTRY_INT4(h->string_length);

    SDF_READ_ENTRY_INT4(h->code_io_version);

    SDF_READ_ENTRY_LOGICAL(h->restart_flag);

    SDF_READ_ENTRY_LOGICAL(h->other_domains);

    free(h->buffer);
    h->buffer = NULL;

    h->current_location = h->first_block_location;
    h->done_header = 1;

    if (h->summary_size == 0) h->use_summary = 0;

    return 0;
}



int sdf_read_summary(sdf_file_t *h)
{
    if (h->blocklist) {
        h->current_block = h->blocklist;
        return 0;
    }

    if (!h->done_header) sdf_read_header(h);
    h->current_block = NULL;

    // Read the whole summary block into a temporary buffer on rank 0
    if (h->use_summary > 0) {
        h->current_location = h->start_location = h->summary_location;
        h->buffer = malloc(h->summary_size);

        if (h->rank == h->rank_master) {
            sdf_seek(h);
            sdf_read_bytes(h, h->buffer, h->summary_size);
        }
    } else {
        build_summary_buffer(h);
    }

    // Send the temporary buffer to all processors
    sdf_broadcast(h, h->buffer, h->summary_size);

    return 0;
}



int sdf_read_blocklist(sdf_file_t *h)
{
    int i;
#ifdef PARALLEL
    int fix;
    sdf_block_t *b, *next, *mesh;
#endif

    sdf_read_summary(h);

    // Construct the metadata blocklist using the contents of the buffer
    for (i = 0; i < h->nblocks; i++) {
        SDF_DPRNT("\n");
        sdf_read_block_info(h);
    }

    free(h->buffer);

    h->buffer = NULL;
    h->current_block = h->blocklist;

#ifdef PARALLEL
    // Hack to fix cartesian blocks whose mesh sizes don't match the stagger
    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;
        if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE && b->stagger) {
            fix = 0;
            mesh = sdf_find_block_by_id(h, b->mesh_id);
            for (i = 0; i < b->ndims; i++) {
                if (b->const_value[i] && b->dims[i] != mesh->dims[i]) {
                    fix = 1;
                    b->const_value[i] = 0;
                }
            }
            if (fix) {
                // Re-calculate per block parallel factorisation
                sdf_factor(h);
            }
        }
    }
#endif
    return 0;
}



int sdf_read_block_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int ret = 0;

    sdf_read_next_block_header(h);
    b = h->current_block;
    if (b->done_info) return 0;

    h->indent += 2;
    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH
            || b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH)
        ret = sdf_read_plain_mesh_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        ret = sdf_read_point_mesh_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE)
        ret = sdf_read_plain_variable_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        ret = sdf_read_point_variable_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_CONSTANT)
        ret = sdf_read_constant(h);
    else if (b->blocktype == SDF_BLOCKTYPE_ARRAY)
        ret = sdf_read_array_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT)
        ret = sdf_read_cpu_split_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_RUN_INFO)
        ret = sdf_read_run_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS
            || b->blocktype == SDF_BLOCKTYPE_STITCHED_TENSOR
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_TENSOR)
        ret = sdf_read_stitched(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL)
        ret = sdf_read_stitched_material(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATVAR
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATVAR)
        ret = sdf_read_stitched_matvar(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_SPECIES
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_SPECIES)
        ret = sdf_read_stitched_species(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP)
        ret = sdf_read_stitched_obstacle_group(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STATION)
        ret = sdf_read_station_info(h);

    return ret;
}



int sdf_read_data(sdf_file_t *h)
{
    sdf_block_t *b;

    b = h->current_block;

    if (b->populate_data) {
        b->populate_data(h, b);
        return 0;
    } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH)
        return sdf_read_plain_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH)
        return sdf_read_lagran_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        return sdf_read_point_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE)
        return sdf_read_plain_variable(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        return sdf_read_point_variable(h);
    else if (b->blocktype == SDF_BLOCKTYPE_ARRAY
            || b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT)
        return sdf_read_array(h);

    return 1;
}



// Read the block header into the current block
static int sdf_read_next_block_header(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;

    if (!h->done_header) {
        if (h->rank == h->rank_master) {
            fprintf(stderr, "*** ERROR ***\n");
            fprintf(stderr, "SDF header has not been read. Ignoring call.\n");
        }
        return 1;
    }

    sdf_get_next_block(h);
    b = h->current_block;

    if (!b) {
        if (h->rank == h->rank_master) {
            fprintf(stderr, "*** ERROR ***\n");
            fprintf(stderr, "SDF block not initialised. Ignoring call.\n");
        }
        return 1;
    }

    if (b->done_header) {
        h->current_location = b->block_start + h->block_header_length;
        return 0;
    }

    h->indent = 2;

    if (h->use_summary)
        h->current_location = b->block_start;

    SDF_READ_ENTRY_INT8(b->next_block_location);

    SDF_READ_ENTRY_INT8(b->data_location);

    SDF_READ_ENTRY_ID(b->id);

    SDF_READ_ENTRY_INT8(b->data_length);

    SDF_READ_ENTRY_TYPE(blocktype);

    SDF_READ_ENTRY_TYPE(datatype);

    SDF_READ_ENTRY_INT4(b->ndims);

    SDF_READ_ENTRY_STRING(b->name);

    // Older versions of the file did not contain the block
    // info length in the header.
    if (h->file_version + h->file_revision > 1)
        SDF_READ_ENTRY_INT4(b->info_length);

    if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE
            || b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        b->stagger = SDF_STAGGER_VERTEX;
    else
        b->stagger = SDF_STAGGER_CELL_CENTRE;
    for (i = 0; i < 3; i++) b->dims[i] = 1;

    if (b->blocktype == SDF_BLOCKTYPE_STATION)
        h->station_file = 1;

    b->done_header = 1;
    h->current_location = b->block_start + h->block_header_length;

    b->datatype_out = b->datatype;
    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8)
        b->datatype_out = SDF_DATATYPE_REAL4;
#ifdef PARALLEL
    switch (b->datatype) {
    case(SDF_DATATYPE_REAL4):
        b->mpitype = MPI_FLOAT;
        break;
    case(SDF_DATATYPE_REAL8):
        b->mpitype = MPI_DOUBLE;
        break;
    case(SDF_DATATYPE_INTEGER4):
        b->mpitype = MPI_INT;
        break;
    case(SDF_DATATYPE_INTEGER8):
        b->mpitype = MPI_LONG_LONG;
        break;
    case(SDF_DATATYPE_CHARACTER):
        b->mpitype = MPI_CHAR;
        break;
    case(SDF_DATATYPE_LOGICAL):
        b->mpitype = MPI_CHAR;
        break;
    }
    b->mpitype_out = b->mpitype;
#endif

    return 0;
}



// Read all the block metadata sections and copy them into
// one contiguous buffer (h->buffer).
static void build_summary_buffer(sdf_file_t *h)
{
    int i, buflen;
    uint64_t data_location, block_location, next_block_location;
    uint32_t info_length;
    char *bufptr;
    char skip_summary;

    struct list_entry {
        void *buffer;
        int len;
        struct list_entry *next;
    } *blockbuf_head, *blockbuf;

    h->current_location = h->first_block_location;

    // Read the file and build the buffer on rank zero.
    if (h->rank == h->rank_master) {
        next_block_location = block_location = h->current_location;

        blockbuf_head = blockbuf = calloc(1,sizeof(*blockbuf));
        buflen = 0;
        h->nblocks = 0;
        skip_summary = (h->summary_location &&
                h->summary_location != h->first_block_location);
        // Read the block metadata into a temporary linked list structure
        while (1) {
            if (skip_summary &&
                h->current_location >= h->summary_location) break;

            sdf_seek(h);

            // Read the fixed length block header
            blockbuf->len = h->block_header_length;
            blockbuf->buffer = malloc(blockbuf->len);
            i = sdf_read_bytes(h, blockbuf->buffer, blockbuf->len);
            if (i != 0) break;

            memcpy(&next_block_location, blockbuf->buffer, sizeof(uint64_t));
            memcpy(&data_location, (char*)blockbuf->buffer+8, sizeof(uint64_t));

            if (h->swap) {
                _SDF_BYTE_SWAP64(next_block_location);
                _SDF_BYTE_SWAP64(data_location);
            }

            // Older versions of the file did not contain the block
            // info length in the header.
            if (h->file_version + h->file_revision > 1) {
                memcpy(&info_length, (char*)blockbuf->buffer+132,
                       sizeof(uint32_t));
                if (h->swap) _SDF_BYTE_SWAP32(info_length);
            } else {
                if (data_location > block_location)
                    info_length = (uint32_t)(data_location
                        - block_location) - h->block_header_length;
                else
                    info_length = (uint32_t)(next_block_location
                        - block_location) - h->block_header_length;
            }

            // Read the block specific metadata if it exists
            if (info_length > 0) {
                blockbuf->next = calloc(1,sizeof(*blockbuf));
                blockbuf = blockbuf->next;

                blockbuf->len = info_length;
                blockbuf->buffer = malloc(blockbuf->len);
                i = sdf_read_bytes(h, blockbuf->buffer, blockbuf->len);
                if (i != 0) break;
            }

            buflen += h->block_header_length + info_length;

            h->nblocks++;

            blockbuf->next = calloc(1,sizeof(*blockbuf));
            blockbuf = blockbuf->next;

            if (h->current_location > next_block_location) break;

            h->current_location = block_location = next_block_location;
        }

        if (blockbuf->buffer) {
            free(blockbuf->buffer);
            blockbuf->buffer = NULL;
        }

        // Copy the contents of the linked list into a single contiguous buffer.
        bufptr = h->buffer = malloc(buflen);
        while (blockbuf_head) {
            blockbuf = blockbuf_head;

            if (blockbuf->buffer == NULL) {
                free(blockbuf);
                break;
            }

            memcpy(bufptr, blockbuf->buffer, blockbuf->len);
            bufptr += blockbuf->len;

            blockbuf_head = blockbuf->next;

            free(blockbuf->buffer);
            free(blockbuf);
        }
    }

    // Send the temporary buffer length to all processors
    sdf_broadcast(h, &buflen, sizeof(buflen));

    // Allocate the buffer on all non rank zero processors
    if (h->rank != h->rank_master)
        h->buffer = malloc(buflen);

    h->summary_size = buflen;
    h->current_location = 0;
}



static int sdf_read_stitched(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_material(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matnames  ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_matvar(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ID(b->material_id);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_species(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - matname   CHARACTER(string_length)
    // - specnames ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ID(b->material_id);

    SDF_READ_ENTRY_STRING(b->material_name);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_obstacle_group(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger         INTEGER(i4)
    // - obstacle_id     CHARACTER(id_length)
    // - vfm_id          CHARACTER(id_length)
    // - obstacle_names  ndims*CHARACTER(string_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->obstacle_id);

    SDF_READ_ENTRY_ID(b->vfm_id);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);

    b->done_data = 1;

    return 0;
}



static int sdf_read_constant(sdf_file_t *h)
{
    sdf_block_t *b;

    // Metadata is
    // - value     TYPE_SIZE

    b = h->current_block;
    SDF_READ_ENTRY_CONST(b->const_value);

    b->stagger = SDF_STAGGER_VERTEX;

    return 0;
}



static int sdf_array_datatype(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

#ifdef PARALLEL
    int sizes[SDF_MAXDIMS];
    for (n=0; n < b->ndims; n++) sizes[n] = (int)b->dims[n];

    sdf_factor(h);

    MPI_Type_create_subarray(b->ndims, sizes, b->local_dims, b->starts,
        MPI_ORDER_FORTRAN, b->mpitype, &b->distribution);
    MPI_Type_commit(&b->distribution);
#else
    for (n=0; n < b->ndims; n++) b->local_dims[n] = (int)b->dims[n];
#endif
    for (n=b->ndims; n < 3; n++) b->local_dims[n] = 1;

    b->nelements_local = 1;
    for (n=0; n < b->ndims; n++) b->nelements_local *= b->local_dims[n];

    return 0;
}



static int sdf_read_run_info(sdf_file_t *h)
{
    sdf_block_t *b;
    struct run_info *run;

    run = calloc(1, sizeof(*run));

    // Metadata is
    // - version   INTEGER(i4)
    // - revision  INTEGER(i4)
    // - commit_id CHARACTER(string_length)
    // - sha1sum   CHARACTER(string_length)
    // - compmac   CHARACTER(string_length)
    // - compflag  CHARACTER(string_length)
    // - defines   INTEGER(i8)
    // - compdate  INTEGER(i4)
    // - rundate   INTEGER(i4)
    // - iodate    INTEGER(i4)
    // - minor_rev INTEGER(i4)

    SDF_COMMON_INFO();
    SDF_READ_ENTRY_INT4(run->version);
    SDF_READ_ENTRY_INT4(run->revision);
    SDF_READ_ENTRY_STRING(run->commit_id);
    SDF_READ_ENTRY_STRING(run->sha1sum);
    SDF_READ_ENTRY_STRING(run->compile_machine);
    SDF_READ_ENTRY_STRING(run->compile_flags);
    SDF_READ_ENTRY_INT8(run->defines);
    SDF_READ_ENTRY_INT4(run->compile_date);
    SDF_READ_ENTRY_INT4(run->run_date);
    SDF_READ_ENTRY_INT4(run->io_date);
    if (h->file_version == 1 && h->file_revision < 2)
        run->minor_rev = 0;
    else
        SDF_READ_ENTRY_INT4(run->minor_rev);

    h->current_block->data = run;
    b->done_data = 1;

    return 0;
}



static int sdf_read_array_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    uint32_t dims_in[SDF_MAXDIMS];
    uint32_t *dims_ptr = dims_in;

    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    b->nelements_local = 1;
    for (i = 0; i < b->ndims; i++) {
        b->local_dims[i] = b->dims[i] = dims_in[i];
        b->nelements_local *= b->dims[i];
    }

    return 0;
}



static int sdf_read_cpu_split_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    uint32_t dims_in[SDF_MAXDIMS];
    uint32_t *dims_ptr = dims_in;

    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_INT4(b->geometry);
    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    for (i = 0; i < b->ndims; i++) {
        b->local_dims[i] = b->dims[i] = dims_in[i];
    }
    if (b->geometry == 1 || b->geometry == 4) {
        b->nelements_local = 0;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local += b->dims[i];
    } else if (b->geometry == 2) {
        b->nelements_local = (int)(b->dims[0] * (b->dims[1] + 1));
        if (b->ndims > 2)
            b->nelements_local += b->dims[0] * b->dims[1] * b->dims[2];
    } else if (b->geometry == 3) {
        b->nelements_local = 1;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local *= b->dims[i];
    }

    return 0;
}



static int sdf_read_array(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    char *p;
    int n, i, count;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_array_info(h);

    h->current_location = b->data_location;

    n = SDF_TYPE_SIZES[b->datatype] * b->nelements_local;
    if (h->mmap) {
        b->data = h->mmap + h->current_location;
    } else {
        if (b->data) free(b->data);
        b->data = malloc(n);
        sdf_seek(h);
        sdf_read_bytes(h, b->data, n);
        if (h->swap) sdf_convert_array_to_float(h, &b->data, n);
    }

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n<b->ndims; n++) SDF_DPRNT("%i ",b->local_dims[n]);
        SDF_DPRNT("\n  ");
        if (b->datatype_out == SDF_DATATYPE_CHARACTER) {
            p = b->data;
            for (n=0; n < b->local_dims[1]; n++) {
                p = (char*)b->data + n * b->local_dims[0];
                count = 0;
                for (i=0; i < b->local_dims[0]; i++) {
                    if (*p == '\0') break;
                    if (*p != ' ') count++;
                    p++;
                }
                SDF_DPRNT("c*%i[%i] ", b->local_dims[0], n);
                p = (char*)b->data + n * b->local_dims[0];
                for (i=0; i < count; i++) {
                    SDF_DPRNT("%c", *p);
                    p++;
                }
                SDF_DPRNT("\n  ");
            }
        } else {
            SDF_DPRNTar(b->data, b->nelements_local);
        }
    }

    b->done_data = 1;

    return 0;
}
