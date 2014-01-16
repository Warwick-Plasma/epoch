#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sdf.h>
#ifdef PARALLEL
#include <mpi.h>
#endif

static int write_header(sdf_file_t *h);
static int write_constant(sdf_file_t *h);
static int write_stitched(sdf_file_t *h);
static int write_stitched_material(sdf_file_t *h);
static int write_stitched_matvar(sdf_file_t *h);
static int write_run_info_meta(sdf_file_t *h);
static int write_data(sdf_file_t *h);
int sdf_safe_write_id(sdf_file_t *h, char *string);
int sdf_safe_write_string(sdf_file_t *h, char *string);



static int sdf_write_bytes(sdf_file_t *h, void *buf, int buflen)
{
#ifdef PARALLEL
    return MPI_File_write(h->filehandle, buf, buflen, MPI_BYTE,
            MPI_STATUS_IGNORE);
#else
    return (1 != fwrite(buf, buflen, 1, h->filehandle));
#endif
}



static int sdf_flush(sdf_file_t *h)
{
    //errcode += sdf_update(h);
#ifdef PARALLEL
    return MPI_File_sync(h->filehandle);
#else
    return fflush(h->filehandle);
#endif
}



#define MIN(a,b) (((a) < (b)) ? (a) : (b))

static int safe_copy_string(char *s1, char *s2)
{
    int len1, len2;

    len1 = strlen(s1);
    len2 = strlen(s2);

    memset(s2, 0, len2);
    memcpy(s2, s1, len1);

    return 0;
}



static int write_block(sdf_file_t *h)
{
    sdf_block_t *b;

    b = h->current_block;
    b->block_start = h->current_location;

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_CONSTANT:
        write_constant(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
        write_stitched_material(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
        write_stitched_matvar(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_TENSOR:
        write_stitched(h);
        break;
    case SDF_BLOCKTYPE_RUN_INFO:
        write_run_info_meta(h);
        break;
    case SDF_BLOCKTYPE_ARRAY:
    case SDF_BLOCKTYPE_PLAIN_MESH:
        write_data(h);
        break;
    default:
        printf("ignored id: %s\n", b->id);
    }
#if 0
    SDF_BLOCKTYPE_PLAIN_MESH:
    SDF_BLOCKTYPE_POINT_MESH,
    SDF_BLOCKTYPE_PLAIN_VARIABLE,
    SDF_BLOCKTYPE_POINT_VARIABLE,
    //SDF_BLOCKTYPE_CONSTANT,
    //SDF_BLOCKTYPE_ARRAY,
    SDF_BLOCKTYPE_RUN_INFO,
    SDF_BLOCKTYPE_SOURCE,
    SDF_BLOCKTYPE_STITCHED_TENSOR,
    //SDF_BLOCKTYPE_STITCHED_MATERIAL,
    SDF_BLOCKTYPE_STITCHED_MATVAR,
    SDF_BLOCKTYPE_STITCHED_SPECIES,
    SDF_BLOCKTYPE_SPECIES,
    SDF_BLOCKTYPE_PLAIN_DERIVED,
    SDF_BLOCKTYPE_POINT_DERIVED,
    SDF_BLOCKTYPE_CONTIGUOUS_TENSOR,
    SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL,
    SDF_BLOCKTYPE_CONTIGUOUS_MATVAR,
    SDF_BLOCKTYPE_CONTIGUOUS_SPECIES,
    SDF_BLOCKTYPE_CPU_SPLIT,
    SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP,
    SDF_BLOCKTYPE_UNSTRUCTURED_MESH,
    SDF_BLOCKTYPE_STITCHED,
    SDF_BLOCKTYPE_CONTIGUOUS,
    SDF_BLOCKTYPE_LAGRANGIAN_MESH,
    SDF_BLOCKTYPE_STATION,
    SDF_BLOCKTYPE_STATION_DERIVED,
#endif
    return 0;
}



int sdf_write(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b;

    errcode = write_header(h);
    sdf_flush(h);

    b = h->blocklist;
    while (b) {
        h->current_block = b;
        b->done_header = b->done_info = b->done_data = 0;
        errcode += write_block(h);
 sdf_flush(h);
        b = b->next;
    }

    return errcode;
}



static int write_header(sdf_file_t *h)
{
    uint32_t int4;
    int errcode;
    char padding[6];

    errcode = 0;

    if (h->done_header) {
        if (h->rank == h->rank_master) {
            printf("*** WARNING ***\n");
            printf("SDF header already written. Ignoring extra call.\n");
        }
        return errcode;
    }

    // header length - must be updated if sdf_write_header changes
    h->first_block_location = SDF_HEADER_LENGTH;
    // block header length - must be updated if sdf_write_block_header changes
    h->block_header_length = SDF_BLOCK_HEADER_LENGTH;

    // Currently no blocks written
    h->nblocks = 0;
    h->summary_location = h->first_block_location;
    h->summary_size = 0;
    h->current_location = 0;
    //h->data_location = 0;

    if (h->rank == h->rank_master) {
        errcode += sdf_seek(h);

        // Write the header
        errcode += sdf_write_bytes(h, SDF_MAGIC, 4);

        int4 = SDF_ENDIANNESS;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        int4 = SDF_VERSION;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        int4 = SDF_REVISION;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        errcode += sdf_safe_write_id(h, h->code_name);

        errcode += sdf_write_bytes(h, &h->first_block_location, SOI8);

        // Must be consistent with the c_summary_offset in sdf_common
        errcode += sdf_write_bytes(h, &h->summary_location, SOI8);

        errcode += sdf_write_bytes(h, &h->summary_size, SOI4);

        errcode += sdf_write_bytes(h, &h->nblocks, SOI4);

        errcode += sdf_write_bytes(h, &h->block_header_length, SOI4);

        errcode += sdf_write_bytes(h, &h->step, SOI4);

        errcode += sdf_write_bytes(h, &h->time, SOF8);

        errcode += sdf_write_bytes(h, &h->jobid1, SOI4);

        errcode += sdf_write_bytes(h, &h->jobid2, SOI4);

        errcode += sdf_write_bytes(h, &h->string_length, SOI4);

        errcode += sdf_write_bytes(h, &h->code_io_version, SOI4);

        errcode += sdf_write_bytes(h, &h->restart_flag, 1);

        errcode += sdf_write_bytes(h, &h->other_domains, 1);

        memset(padding, 0, 6);
        errcode += sdf_write_bytes(h, padding, 6);
    }

    h->current_location = h->first_block_location;
    h->done_header = 1;

    return 0;
}



int sdf_write_header(sdf_file_t *h, char *code_name, int code_io_version,
        int step, double time, char restart, int jobid1, int jobid2)
{
    int errcode;

    errcode = 0;

    // header length - must be updated if sdf_write_header changes
    h->first_block_location = SDF_HEADER_LENGTH;
    // block header length - must be updated if sdf_write_block_header changes
    h->block_header_length = SDF_BLOCK_HEADER_LENGTH;

    if (h->code_name) free(h->code_name);
    errcode += safe_copy_string(code_name, h->code_name);

    h->step = step;
    h->time = time;
    h->restart_flag = restart;
    h->jobid1 = jobid1;
    h->jobid2 = jobid2;

    return write_header(h);
}



static int write_block_header(sdf_file_t *h)
{
    sdf_block_t *b;
    int errcode, block_info_length;

    errcode = 0;

    b = h->current_block;
    if (b->done_header) return errcode;

    b->data_location = b->block_start + b->info_length;
    b->next_block_location = b->data_location + b->data_length;

    // If this routine is changed then the value of h->block_header_length
    // must be changed accordingly in sdf_write_header

    if (h->rank == h->rank_master) {
        h->current_location = b->block_start;

        errcode += sdf_seek(h);

        // Write the block header
        errcode += sdf_write_bytes(h, &b->next_block_location, SOI8);

        errcode += sdf_write_bytes(h, &b->data_location, SOI8);

        errcode += sdf_safe_write_id(h, b->id);

        errcode += sdf_write_bytes(h, &b->data_length, SOI8);

        errcode += sdf_write_bytes(h, &b->blocktype, SOI4);

        errcode += sdf_write_bytes(h, &b->datatype, SOI4);

        errcode += sdf_write_bytes(h, &b->ndims, SOI4);

        errcode += sdf_safe_write_string(h, b->name);

        block_info_length = b->info_length - h->block_header_length;
        errcode += sdf_write_bytes(h, &block_info_length, SOI4);
    }

    h->current_location = b->block_start + h->block_header_length;
    b->done_header = 1;

    return errcode;
}



size_t trimwhitespace(const char *str_in, char *str_out, size_t len)
{
    size_t out_size;
    const char *end;

    if (len == 0) return 0;

    // Trim leading space
    while (isspace(*str_in)) str_in++;

    if (*str_in == 0) { // All spaces?
        *str_out = 0;
        return 1;
    }

    // Trim trailing space
    end = str_in + strlen(str_in) - 1;
    while (end > str_in && isspace(*end)) end--;
    end++;

    // Set output size to minimum of trimmed string length and
    // buffer size minus 1
    out_size = (end - str_in) < len-1 ? (end - str_in) : len-1;

    // Copy trimmed string and add null terminator
    memcpy(str_out, str_in, out_size);
    str_out[out_size] = 0;

    return out_size;
}



int sdf_safe_write_string_len(sdf_file_t *h, char *string, int length)
{
    char output[length];
    int len_s;

    len_s = trimwhitespace(string, output, length);

    if (len_s > length && h->rank == h->rank_master) {
        printf("*** WARNING ***\n");
        printf("Output string \"%s\" has been truncated.", output);
    }

    if (len_s+1 < length) output[len_s+1] = 0;

    return sdf_write_bytes(h, output, length);
}



int sdf_safe_write_string(sdf_file_t *h, char *string)
{
    return sdf_safe_write_string_len(h, string, h->string_length);
}



int sdf_safe_write_id(sdf_file_t *h, char *string)
{
    int length = SDF_ID_LENGTH;
    return sdf_safe_write_string_len(h, string, length);
}



static int write_constant(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b = h->current_block;

    b->nelements = b->ndims = 1;

    b->info_length = h->block_header_length + SDF_TYPE_SIZES[b->datatype];
    b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write data (in metadata section)
    if (h->rank == h->rank_master)
        errcode += sdf_write_bytes(h, &b->const_value,
                                   SDF_TYPE_SIZES[b->datatype]);

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_array_meta(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    b->nelements = 1;
    for (i=0; i < b->ndims; i++)
        b->nelements *= b->dims[i];

    // Metadata is
    // - dims      ndims*INTEGER(i4)

    b->info_length = h->block_header_length + b->ndims * SOI4;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    if (h->rank == h->rank_master)
        errcode += sdf_write_bytes(h, b->dims, b->ndims * SOI4);

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



int write_mesh_meta(sdf_file_t *h)
{
    int errcode, i;
    uint32_t uint4;
    sdf_block_t *b = h->current_block;

    b->nelements = 0;
    for (i=0; i < b->ndims; i++)
        b->nelements += b->dims[i];

    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - dims        INTEGER(i4), DIMENSION(ndims)

    b->info_length = h->block_header_length + (b->ndims + 1) * SOI4 +
        (3 * b->ndims) * SOF8 + 2 * b->ndims * SDF_ID_LENGTH;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, b->dim_mults, b->ndims * SOF8);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_labels[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_units[i]);

        errcode += sdf_write_bytes(h, &b->geometry, SOI4);

        errcode += sdf_write_bytes(h, b->extents, 2 * b->ndims * SOF8);

        for (i=0; i < b->ndims; i++) {
            uint4 = b->dims[i];
            errcode += sdf_write_bytes(h, &uint4, SOI4);
        }
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_stitched(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 +
            (b->ndims + 1) * SDF_ID_LENGTH;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_material(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - material_names ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 +
            (b->ndims + 1) * SDF_ID_LENGTH + b->ndims * h->string_length;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_string(h, b->material_names[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_matvar(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 +
            (b->ndims + 2) * SDF_ID_LENGTH;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_safe_write_id(h, b->material_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_run_info_meta(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b = h->current_block;
    struct run_info *run = b->data;

    b->datatype = SDF_DATATYPE_OTHER;
    b->blocktype = SDF_BLOCKTYPE_RUN_INFO;
    b->ndims = 1;

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

    b->info_length = h->block_header_length + 6 * SOI4 + SOI8 +
            4 * h->string_length;
    b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &run->version, SOI4);

        errcode += sdf_write_bytes(h, &run->revision, SOI4);

        errcode += sdf_safe_write_string(h, run->commit_id);

        errcode += sdf_safe_write_string(h, run->sha1sum);

        errcode += sdf_safe_write_string(h, run->compile_machine);

        errcode += sdf_safe_write_string(h, run->compile_flags);

        errcode += sdf_write_bytes(h, &run->defines, SOI8);

        errcode += sdf_write_bytes(h, &run->compile_date, SOI4);

        errcode += sdf_write_bytes(h, &run->run_date, SOI4);

        errcode += sdf_write_bytes(h, &run->io_date, SOI4);

        errcode += sdf_write_bytes(h, &run->minor_rev, SOI4);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_meta(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_ARRAY:
        return write_array_meta(h);
        break;
    case SDF_BLOCKTYPE_PLAIN_MESH:
        return write_mesh_meta(h);
        break;
    default:
        printf("ignored id: %s\n", b->id);
    }

    return 1;
}



static int write_data(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b;

    errcode = 0;
    b = h->current_block;

    // Write header
    errcode += write_meta(h);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += sdf_seek(h);

        // Actual array
        if (b->data) {
            errcode += sdf_write_bytes(h, b->data, b->data_length);
        } else if (b->grids) {
            for (i=0; i < b->ndims; i++)
                errcode += sdf_write_bytes(h, b->grids[i],
                        b->dims[i] * SDF_TYPE_SIZES[b->datatype]);
        }
    }

    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

    return errcode;
}



#if 0
int sdf_update(sdf_file_t *h)
{
    int offset;
    int errcode;
    int int4;

    // No open file or not writing
    IF (h->filehandle == -1 || !h->writing) return 0;

    // Update summary and nblocks info
    if (h->rank == h->rank_master) {
        if (h->error_code != 0) {
            h->nblocks = -h->error_code;
            h->summary_location = 0;
            h->summary_size = 0;
        }
        if (h->summary_location != h->summary_location_wrote) {
            offset = SDF_SUMMARY_OFFSET;
            errcode += sdf_write_at(h, offset, &h->summary_location, 8);
            h->summary_location_wrote = h->summary_location;
        }
        if (h->summary_size != h->summary_size_wrote) {
            offset = SDF_SUMMARY_OFFSET + 8;
            errcode += sdf_write_at(h, offset, &h->summary_size, 4);
            h->summary_size_wrote = h->summary_size;
        }
        if (h->nblocks != h->nblocks_wrote) {
            offset = SDF_SUMMARY_OFFSET + 12;
            errcode += sdf_write_at(h, offset, &h->nblocks, 4);
            h->nblocks_wrote = h->nblocks;
        }
        if (h->step != h->step_wrote) {
            offset = SDF_SUMMARY_OFFSET + 20;
            int4 = INT(h->step,i4);
            errcode += sdf_write_at(h, offset, &int4, 4);
            h->step_wrote = h->step;
        }
        if (h->time != h->time_wrote) {
            offset = SDF_SUMMARY_OFFSET + 24;
            errcode += sdf_write_at(h, offset, &h->time, 8);
            h->time_wrote = h->time;
        }
    }

}


int sdf_write_block_header(h, id, name);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    TYPE(sdf_block_type), POINTER :: b;

    b => h->current_block;

    // If this routine is changed then the value of h->block_header_length
    // must be changed accordingly in sdf_open_clobber

    if (! h->done_header) {
        if (h->rank == h->rank_master) {
              PRINT*,'*** WARNING ***';
              PRINT*,'SDF header not yet written. Ignoring write call.';
        }
        return 0;
    }

    if (b->done_header) {
        if (h->rank == h->rank_master) {
              PRINT*,'*** WARNING ***';
              PRINT*,'SDF block header already written. Ignoring extra call.';
        }
        return 0;
    }

    errcode += safe_copy_unique_id(h, b, id);
    errcode += safe_copy_string(name, b->name);

    errcode += write_block_header(h);

    h->nblocks = h->nblocks + 1_4;

}



int sdf_safe_write_string_len(h, string, length);

    TYPE(sdf_file_handle) :: h;
    char *string;
    int length;
    char *output;
    int len_s, errcode;

    len_s = LEN_TRIM(string);

    if (len_s .GT. length .AND. h->rank == h->rank_master) {
        PRINT*, '*** WARNING ***';
        PRINT*, 'Output string "' // TRIM(string) // '" has been truncated';
    }

    // Thisint expects that the record marker is in place and that
    // the view is set correctly. errcode += it only on the node which is doing the
    // writing. You still have to advance the file pointer yourself on all nodes

    output = ' ';
    output(1:MIN(length, len_s)) = string(1:MIN(length, len_s));

    // If this isn't the full string length then tag in a ACHAR(0) to help
    // With C++ string handling
    IF (len_s + 1 .LT. length) output(len_s+1:length) = ACHAR(0);

    errcode += MPI_FILE_WRITE(h->filehandle, output, length,
              MPI_CHARACTER, MPI_STATUS_IGNORE, errcode);

}



int sdf_safe_write_string(h, string, length_in);

    TYPE(sdf_file_handle) :: h;
    char *string;
    int length_in;
    int length;

    if (PRESENT(length_in)) {
        length = length_in;
    } else {
        length = h->string_length;
    }

    errcode += sdf_safe_write_string_len(h, string, length);

}



int sdf_safe_write_id(h, string);

    TYPE(sdf_file_handle) :: h;
    char *string;

    errcode += sdf_safe_write_string_len(h, string, INT(c_id_length));

}



  FUNCTION sdf_string_lowercase(string_in) RESULT(string_out);

    char *lwr = 'abcdefghijklmnopqrstuvwxyz';
    char *upr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    char *string_in;
    char *string_out;
    int i, idx;

    string_out = string_in;

    for (i = 1, LEN(string_out)) {
        idx = INDEX(upr, string_out(i:i));
        IF (idx != 0) string_out(i:i) = lwr(idx:idx);
    }

  END FUNCTION sdf_string_lowercase;



int sdf_safe_string_composite(h, string1, string2, output_string);

    TYPE(sdf_file_handle) :: h;
    char *string1, string2;
    char *output_string;
    int len1, len2;

    len1 = LEN_TRIM(string1);
    len2 = LEN_TRIM(string2);
    olen = LEN(output_string);

    output_string = '';

    if (olen < len1 + 1) {
        output_string(1:olen) = string1(1:olen);
    ELSE if (olen < len1 + len2 + 1) {
        output_string(1:len1) = string1(1:len1);
        output_string(len1+1:len1+1) = '/';
        output_string(len1+2:olen) = string2(1:olen-len1-1);
    } else {
        output_string(1:len1) = string1(1:len1);
        output_string(len1+1:len1+1) = '/';
        output_string(len1+2:len1+len2+1) = string2(1:len2);
    }

}



int write_run_info_minor(h, version, revision, minor_rev, commit_id,
        sha1sum, compile_machine, compile_flags, defines, compile_date,
        run_date, rank_write);

    TYPE(sdf_file_handle) :: h;
    int version, revision, minor_rev;
    char *commit_id, sha1sum;
    char *compile_machine, compile_flags;
    int defines;
    int compile_date, run_date;
    int rank_write;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    IF (! ASSOCIATED(b->run)) ALLOCATE(b->run);
    b->run->version = version;
    b->run->revision = revision;
    b->run->minor_rev = minor_rev;
    errcode += safe_copy_string(commit_id, b->run->commit_id);
    errcode += safe_copy_string(sha1sum, b->run->sha1sum);
    errcode += safe_copy_string(compile_machine, b->run->compile_machine);
    errcode += safe_copy_string(compile_flags, b->run->compile_flags);
    b->run->defines = defines;
    b->run->compile_date = compile_date;
    b->run->run_date = run_date;
    b->run->io_date = get_unix_time();

    errcode += write_run_info_meta(h, 'run_info', 'Run_info');

    h->rank_master = h->default_rank;

}



int write_run_info_old(h, version, revision, commit_id,
        sha1sum, compile_machine, compile_flags, defines, compile_date,
        run_date, rank_write);

    TYPE(sdf_file_handle) :: h;
    int version, revision;
    char *commit_id, sha1sum;
    char *compile_machine, compile_flags;
    int defines;
    int compile_date, run_date;
    int rank_write;

    errcode += write_run_info_minor(h, version, revision, 0, commit_id,
        sha1sum, compile_machine, compile_flags, defines, compile_date,
        run_date, rank_write);

}



int write_1d_array_integer_spec(h, id, name, n1, array, rank_write);

    int ndims = 1;
    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int n1;
    int array;
    int rank_write;
    int errcode;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;
    b->ndims = ndims;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = n1;

    // Write header

    errcode += write_array_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        errcode += MPI_FILE_WRITE(h->filehandle, array, n1, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

}



int write_1d_array_integer(h, id, name, array, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int array;
    int rank_write;
    int n1;

    n1 = SIZE(array,1);
    errcode += write_1d_array_integer_spec(h, id, name, n1, array, rank_write);

}



int sdf_write_stitched(h, id, name, mesh_id, stagger,
        variable_ids, ndims, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *variable_ids(:);
    int ndims;
    int data_length;
    TYPE(sdf_block_type), POINTER :: b;

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(data_length)) {
              b->data_length = data_length;
              b->blocktype = SDF_BLOCKTYPE_CONTIGUOUS;
        } else {
              b->data_length = 0;
              b->blocktype = SDF_BLOCKTYPE_STITCHED;
        }
    }

    errcode += write_stitched(h, id, name, mesh_id, stagger, variable_ids, ndims,
              data_length);

}



int sdf_write_stitched_tensor(h, id, name, mesh_id, stagger,
        variable_ids, ndims, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *variable_ids(:);
    int ndims;
    int data_length;
    TYPE(sdf_block_type), POINTER :: b;

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(data_length)) {
              b->data_length = data_length;
              b->blocktype = SDF_BLOCKTYPE_CONTIGUOUS_TENSOR;
        } else {
              b->data_length = 0;
              b->blocktype = SDF_BLOCKTYPE_STITCHED_TENSOR;
        }
    }

    errcode += write_stitched(h, id, name, mesh_id, stagger, variable_ids, ndims,
              data_length);

}



int sdf_write_stitched_tensor_mat(h, id, name, mesh_id, stagger,
        variable_ids, material_names, ndims_in, nmat_in, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *variable_ids(:);
    char *material_names(:);
    int ndims_in, nmat_in;
    int data_length;
    int i, j, ndims, nmat;
    int maxstring = 512;
    char *temp_name;
    char *ids;
    char *new_variable_ids;

    if (PRESENT(ndims_in)) {
        ndims = ndims_in;
    } else {
        ndims = INT(SIZE(variable_ids),i4);
    }

    if (PRESENT(nmat_in)) {
        nmat = nmat_in;
    } else {
        nmat = INT(SIZE(material_names),i4);
    }

    ALLOCATE(ids(nmat));

    for (i = 1,nmat) {
        if (LEN_TRIM(material_names(i)) == 0) {
              ids(i) = '';
        } else {
              errcode += sdf_safe_string_composite(h, id,
                        sdf_string_lowercase(material_names(i)), ids(i));
        }
    }

    ALLOCATE(new_variable_ids(ndims));
    for (i = 1,nmat) {
        IF (LEN_TRIM(material_names(i)) == 0) CYCLE;
        for (j = 1,ndims) {
              errcode += sdf_safe_string_composite(h, variable_ids(j),
                        sdf_string_lowercase(material_names(i)), new_variable_ids(j));
        }
        errcode += sdf_safe_string_composite(h, name, material_names(i), temp_name);
        errcode += sdf_write_stitched_tensor(h, ids(i), temp_name, mesh_id,
                      stagger, new_variable_ids, ndims, data_length);
    }

    DEALLOCATE(ids);
    DEALLOCATE(new_variable_ids);

}



int sdf_write_stitched_material(h, id, name, mesh_id, stagger,
        material_names, variable_ids, ndims, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *material_names(:), variable_ids(:);
    int ndims;
    int data_length;
    int i, errcode;
    TYPE(sdf_block_type), POINTER :: b;

    if (h->blocktype == SDF_BLOCKTYPE_CONTIGUOUS) {
        errcode += sdf_write_stitched(h, id, name, mesh_id, stagger,
                      variable_ids, ndims, data_length);
        return 0;
    }

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(ndims)) {
              b->ndims = ndims;
        } else {
              b->ndims = INT(SIZE(variable_ids),i4);
        }
    }

    b => h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - material_names ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 + (b->ndims + 1) * SDF_ID_LENGTH
              + b->ndims * h->string_length;

    // Write header
    if (PRESENT(id)) {
        b->stagger = stagger;
        if (PRESENT(data_length)) {
              b->data_length = data_length;
              b->blocktype = SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL;
        } else {
              b->data_length = 0;
              b->blocktype = SDF_BLOCKTYPE_STITCHED_MATERIAL;
        }
        errcode += safe_copy_id(h, mesh_id, b->mesh_id);
        errcode += sdf_write_block_header(h, id, name);
        ALLOCATE(b->material_names(b->ndims));
        ALLOCATE(b->variable_ids(b->ndims));
        for (i = 1, b->ndims) {
              errcode += safe_copy_string(material_names(i), b->material_names(i));
              errcode += safe_copy_id(h, variable_ids(i), b->variable_ids(i));
        }
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        // Write metadata
        errcode += MPI_FILE_WRITE(h->filehandle, b->stagger, 1, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_string(h, b->material_names(i));
        }

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_id(h, b->variable_ids(i));
        }
    }

    h->rank_master = h->default_rank;
    if (b->data_length .GT. 0) {
        h->current_location = b->next_block_location;
    } else {
        h->current_location = b->block_start + b->info_length;
    }
    b->done_info = 1;
    b->done_data = 1;

}



int sdf_write_stitched_matvar(h, id, name, mesh_id, stagger,
        material_id, variable_ids, ndims, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *material_id, variable_ids(:);
    int ndims;
    int data_length;
    int i, errcode;
    TYPE(sdf_block_type), POINTER :: b;

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(ndims)) {
              b->ndims = ndims;
        } else {
              b->ndims = INT(SIZE(variable_ids),i4);
        }
    }

    b => h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 + (b->ndims + 2) * SDF_ID_LENGTH;

    // Write header
    if (PRESENT(id)) {
        b->stagger = stagger;
        if (PRESENT(data_length)) {
              b->data_length = data_length;
              b->blocktype = SDF_BLOCKTYPE_CONTIGUOUS_MATVAR;
        } else {
              b->data_length = 0;
              b->blocktype = SDF_BLOCKTYPE_STITCHED_MATVAR;
        }
        errcode += safe_copy_id(h, mesh_id, b->mesh_id);
        errcode += safe_copy_id(h, material_id, b->material_id);
        errcode += sdf_write_block_header(h, id, name);
        ALLOCATE(b->variable_ids(b->ndims));
        for (i = 1, b->ndims) {
              errcode += safe_copy_id(h, variable_ids(i), b->variable_ids(i));
        }
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        // Write metadata
        errcode += MPI_FILE_WRITE(h->filehandle, b->stagger, 1, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_safe_write_id(h, b->material_id);

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_id(h, b->variable_ids(i));
        }
    }

    h->rank_master = h->default_rank;
    if (b->data_length .GT. 0) {
        h->current_location = b->next_block_location;
    } else {
        h->current_location = b->block_start + b->info_length;
    }
    b->done_info = 1;
    b->done_data = 1;

}



int sdf_write_stitched_species(h, id, name, mesh_id, stagger,
        material_id, material_name, specnames, variable_ids, ndims, data_length);

    TYPE(sdf_file_handle) :: h;
    char *id, name, mesh_id;
    int stagger;
    char *material_id, material_name;
    char *specnames(:), variable_ids(:);
    int ndims;
    int data_length;
    int i, errcode;
    TYPE(sdf_block_type), POINTER :: b;

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(ndims)) {
              b->ndims = ndims;
        } else {
              b->ndims = INT(SIZE(variable_ids),i4);
        }
    }

    b => h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - matname   CHARACTER(string_length)
    // - specnames ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 + (b->ndims + 2) * SDF_ID_LENGTH
              + (b->ndims + 1) * h->string_length;

    // Write header
    if (PRESENT(id)) {
        b->stagger = stagger;
        if (PRESENT(data_length)) {
              b->data_length = data_length;
              b->blocktype = SDF_BLOCKTYPE_CONTIGUOUS_SPECIES;
        } else {
              b->data_length = 0;
              b->blocktype = SDF_BLOCKTYPE_STITCHED_SPECIES;
        }
        errcode += safe_copy_id(h, mesh_id, b->mesh_id);
        errcode += safe_copy_id(h, material_id, b->material_id);
        errcode += safe_copy_string(material_name, b->material_name);
        ALLOCATE(b->material_names(b->ndims));
        ALLOCATE(b->variable_ids(b->ndims));
        for (i = 1, b->ndims) {
              errcode += safe_copy_string(specnames(i), b->material_names(i));
              errcode += safe_copy_id(h, variable_ids(i), b->variable_ids(i));
        }
        errcode += sdf_write_block_header(h, id, name);
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        // Write metadata
        errcode += MPI_FILE_WRITE(h->filehandle, b->stagger, 1, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_safe_write_id(h, b->material_id);

        errcode += sdf_safe_write_string(h, b->material_name);

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_string(h, b->material_names(i));
        }

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_id(h, b->variable_ids(i));
        }
    }

    h->rank_master = h->default_rank;
    if (b->data_length .GT. 0) {
        h->current_location = b->next_block_location;
    } else {
        h->current_location = b->block_start + b->info_length;
    }
    b->done_info = 1;
    b->done_data = 1;

}



int sdf_write_stitched_obstacle_group(h, id, name, obstacle_id,
        vfm_id, stagger, obstacle_names, ndims, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name, obstacle_id, vfm_id;
    int stagger;
    char *obstacle_names(:);
    int ndims;
    int rank_write;
    int i, errcode;
    TYPE(sdf_block_type), POINTER :: b;

    if (PRESENT(id)) {
        errcode += sdf_get_next_block(h);
        b => h->current_block;
        if (PRESENT(ndims)) {
              b->ndims = ndims;
        } else {
              b->ndims = INT(SIZE(obstacle_names),i4);
        }
    }

    b => h->current_block;

    b->datatype = SDF_DATATYPE_OTHER;
    b->blocktype = SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    // Metadata is
    // - stagger              INTEGER(i4)
    // - obstacle_id    CHARACTER(id_length)
    // - vfm_id               CHARACTER(id_length)
    // - obstacle_names ndims*CHARACTER(string_length)

    b->info_length = h->block_header_length + SOI4 + 2 * SDF_ID_LENGTH
              + b->ndims * h->string_length;
    b->data_length = 0;

    // Write header
    if (PRESENT(id)) {
        b->stagger = stagger;
        errcode += safe_copy_id(h, obstacle_id, b->obstacle_id);
        errcode += safe_copy_id(h, vfm_id, b->vfm_id);
        errcode += sdf_write_block_header(h, id, name);
        ALLOCATE(b->material_names(b->ndims));
        for (i = 1, b->ndims) {
              errcode += safe_copy_string(obstacle_names(i), b->material_names(i));
        }
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        // Write metadata
        errcode += MPI_FILE_WRITE(h->filehandle, b->stagger, 1, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);

        errcode += sdf_safe_write_id(h, b->obstacle_id);
        errcode += sdf_safe_write_id(h, b->vfm_id);

        for (i = 1, b->ndims) {
              errcode += sdf_safe_write_string(h, b->material_names(i));
        }
    }

    h->rank_master = h->default_rank;
    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_constant_meta(h, id, name);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int errcode, var_len;
    TYPE(sdf_block_type), POINTER :: b;

    b => h->current_block;

    b->blocktype = SDF_BLOCKTYPE_CONSTANT;
    b->ndims = 1;

    var_len = 1;
    b->nelements = var_len;
    b->info_length = h->block_header_length + SDF_TYPE_SIZES[b->datatype];
    b->data_length = 0;

    // Write header
    if (PRESENT(id)) {
        errcode += sdf_write_block_header(h, id, name);
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        // Write data (in metadata section)
        errcode += MPI_FILE_WRITE(h->filehandle, b->const_value, var_len, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_constant_integer(h, id, name, value, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int value;
    int rank_write;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->const_value(1:h->SOI) = TRANSFER(value, b->const_value(1:h->SOI));

    errcode += write_constant_meta(h, id, name);

    h->rank_master = h->default_rank;

}



int write_constant_logical(h, id, name, value, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    LOGICAL, INTENT(IN) :: value;
    int rank_write;
    char *cvalue;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = 1;
    b->datatype = SDF_DATATYPE_LOGICAL;
    b->mpitype = MPI_CHARACTER;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    if (value) {
        cvalue = ACHAR(1);
    } else {
        cvalue = ACHAR(0);
    }

    b->const_value(1:1) = TRANSFER(cvalue, b->const_value(1:1));

    errcode += write_constant_meta(h, id, name);

    h->rank_master = h->default_rank;

}



int sdf_write_source_code(h, id, name, array, last, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    char *array;
    char *last;
    int rank_write;
    int i, sz;
    int errcode, len1, len2;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = 1;
    b->datatype = SDF_DATATYPE_CHARACTER;
    b->mpitype = MPI_CHARACTER;
    b->blocktype = SDF_BLOCKTYPE_SOURCE;
    b->ndims = 0;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    if (h->rank == h->rank_master) {
        sz   = SIZE(array);
        len1 = LEN(array);
        len2 = LEN(last);
        b->info_length = h->block_header_length;
        b->data_length = sz*len1 + len2;
    }

    errcode += MPI_BCAST(b->data_length, 1, MPI_INTEGER8, 0, h->comm, errcode);

    // Write header
    errcode += sdf_write_block_header(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Write data
        for (i = 1, sz) {
              errcode += MPI_FILE_WRITE(h->filehandle, array(i), len1,
                        b->mpitype, MPI_STATUS_IGNORE, errcode);
        }
        errcode += MPI_FILE_WRITE(h->filehandle, last, len2,
                      b->mpitype, MPI_STATUS_IGNORE, errcode);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_cpu_split_meta(h, id, name);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int errcode, ndims;
    TYPE(sdf_block_type), POINTER :: b;

    b => h->current_block;

    b->blocktype = SDF_BLOCKTYPE_CPU_SPLIT;
    ndims = b->ndims;

    // Metadata is
    // - type        INTEGER(i4)
    // - dims        ndims*INTEGER(i4)

    b->info_length = h->block_header_length + (b->ndims + 1) * SOI4;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header

    if (PRESENT(id)) {
        errcode += sdf_write_block_header(h, id, name);
    } else {
        errcode += write_block_header(h);
    }

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_WRITE(h->filehandle, b->geometry, 1, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);
        errcode += MPI_FILE_WRITE(h->filehandle, b->dims, ndims, MPI_INTEGER4,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

}



int write_cpu_split_part(h, id, name, npart, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int npart(:);
    int rank_write;
    int errcode;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = SOI8;
    b->datatype = SDF_DATATYPE_INTEGER8;
    b->mpitype = MPI_INTEGER8;
    b->geometry = 4;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = SIZE(npart);
    b->ndims = 1;
    b->nelements = b->dims(1);

    // Write header

    errcode += write_cpu_split_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        errcode += MPI_FILE_WRITE(h->filehandle, npart, b->dims(1), b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_cpu_split_1d_spec(h, id, name, ndim1, nmax1, ndim2, nmax2,
                      ndim3, nmax3, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int ndim1, nmax1(:);
    int ndim2, nmax2(:), ndim3, nmax3(:);
    int rank_write;
    int errcode, i;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;
    b->geometry = 1;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = ndim1 - 1;
    b->ndims = 1;
    if (PRESENT(nmax2)) {
        b->dims(2) = ndim2 - 1;
        b->ndims = 2;
    }
    if (PRESENT(nmax3)) {
        b->dims(3) = ndim3 - 1;
        b->ndims = 3;
    }

    b->nelements = 0;
    for (i = 1,b->ndims) {
        b->nelements = b->nelements + b->dims(i);
    }

    // Write header

    errcode += write_cpu_split_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        errcode += MPI_FILE_WRITE(h->filehandle, nmax1, b->dims(1), b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
        if (b->ndims .GT. 1) {
              errcode += MPI_FILE_WRITE(h->filehandle, nmax2, b->dims(2), b->mpitype,
                        MPI_STATUS_IGNORE, errcode);
        }
        if (b->ndims .GT. 2) {
              errcode += MPI_FILE_WRITE(h->filehandle, nmax3, b->dims(3), b->mpitype,
                        MPI_STATUS_IGNORE, errcode);
        }
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_cpu_split_1d(h, id, name, nmax1, nmax2, nmax3, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int nmax1(:);
    int nmax2(:), nmax3(:);
    int rank_write;
    int ndim1, ndim2, ndim3;

    ndim1 = SIZE(nmax1);
    if (PRESENT(nmax3)) {
        ndim2 = SIZE(nmax2);
        ndim3 = SIZE(nmax3);
        errcode += write_cpu_split_1d_spec(h, id, name, ndim1, nmax1, ndim2, nmax2,
                      ndim3, nmax3, rank_write);
    ELSE if (PRESENT(nmax2)) {
        ndim2 = SIZE(nmax2);
        errcode += write_cpu_split_1d_spec(h, id, name, ndim1, nmax1, ndim2, nmax2);
    } else {
        errcode += write_cpu_split_1d_spec(h, id, name, ndim1, nmax1);
    }

}



int write_cpu_split_mix(h, id, name, nmax1, nmax2, nmax3, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int nmax1(:), nmax2(:,:);
    int nmax3(:,:,:);
    int rank_write;
    int errcode, n1, n2, n3, npt;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;
    b->geometry = 2;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    n1 = SIZE(nmax1);
    n2 = SIZE(nmax2,1);
    b->dims(1) = n1;
    b->dims(2) = n2;

    b->ndims = 2;
    b->nelements = n1 * (1 + n2);
    if (PRESENT(nmax3)) {
        n3 = SIZE(nmax3,1);
        b->dims(3) = n3;
        b->ndims = 3;
        b->nelements = b->nelements + n1 * n2 * n3;
    }

    // Write header

    errcode += write_cpu_split_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        npt = n1;
        errcode += MPI_FILE_WRITE(h->filehandle, nmax1, npt, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
        npt = npt * n2;
        errcode += MPI_FILE_WRITE(h->filehandle, nmax2, npt, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
        if (b->ndims .GT. 2) {
              npt = npt * n3;
              errcode += MPI_FILE_WRITE(h->filehandle, nmax3, npt, b->mpitype,
                        MPI_STATUS_IGNORE, errcode);
        }
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_cpu_split_3d(h, id, name, cpu_splits, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int cpu_splits(:,:);
    int rank_write;
    int errcode, npt;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;
    b->geometry = 3;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = SIZE(cpu_splits,1);
    b->dims(2) = SIZE(cpu_splits,2);
    b->ndims = 2;
    npt = b->dims(1) * b->dims(2);

    b->nelements = npt;

    // Write header

    errcode += write_cpu_split_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        errcode += MPI_FILE_WRITE(h->filehandle, cpu_splits, npt, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_info = 1;
    b->done_data = 1;

}



int write_2d_array_integer_spec(h, id, name, n1, n2, array, rank_write);

    int ndims = 2;
    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int n1, n2;
    int array;
    int rank_write;
    int errcode, var_len, i;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = INT(h->SOI,i4);
    b->datatype = h->datatype_integer;
    b->mpitype = h->mpitype_integer;
    b->ndims = ndims;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = n1;
    b->dims(2) = n2;

    // Write header

    errcode += write_array_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        if (n1 == SIZE(array,1)) {
              var_len = INT(b->nelements);
              errcode += MPI_FILE_WRITE(h->filehandle, array, var_len, b->mpitype,
                        MPI_STATUS_IGNORE, errcode);
        } else {
              for (i = 1,n2) {
                      errcode += MPI_FILE_WRITE(h->filehandle, array(1,i), n1, b->mpitype,
                          MPI_STATUS_IGNORE, errcode);
              }
        }
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

}



int write_2d_array_integer(h, id, name, array, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int array;
    int rank_write;
    int n1, n2;

    n1 = SIZE(array,1);
    n2 = SIZE(array,2);
    errcode += write_2d_array_integer_spec(h, id, name, n1, n2, array, rank_write);

}



int write_1d_array_logical_spec(h, id, name, n1, array, rank_write);

    int ndims = 1;
    TYPE(sdf_file_handle) :: h;
    char *id, name;
    int n1;
    LOGICAL, DIMENSION(:), INTENT(IN) :: array;
    int rank_write;
    char *carray;
    int errcode, i;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = 1;
    b->datatype = SDF_DATATYPE_LOGICAL;
    b->mpitype = MPI_CHARACTER;
    b->ndims = ndims;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    b->dims(1) = n1;

    // Write header

    errcode += write_array_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        ALLOCATE(carray(n1));
        for (i = 1,n1) {
              if (array(i)) {
                      carray(i) = ACHAR(1);
              } else {
                      carray(i) = ACHAR(0);
              }
        }

        // Actual array
        errcode += MPI_FILE_WRITE(h->filehandle, carray, n1, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);

        DEALLOCATE(carray);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

}



int write_1d_array_logical(h, id, name, array, rank_write);

    TYPE(sdf_file_handle) :: h;
    char *id, name;
    LOGICAL, DIMENSION(:), INTENT(IN) :: array;
    int rank_write;
    int n1;

    n1 = SIZE(array,1);
    errcode += write_1d_array_logical_spec(h, id, name, n1, array, rank_write);

}



int write_2d_array_character(h, id, name, array, rank_write);

    int ndims = 2;
    TYPE(sdf_file_handle) :: h;
    char *id, name;
    char *array;
    int rank_write;
    int errcode, var_len, n1, n2;
    TYPE(sdf_block_type), POINTER :: b;

    errcode += sdf_get_next_block(h);
    b => h->current_block;

    b->type_size = 1;
    b->datatype = SDF_DATATYPE_CHARACTER;
    b->mpitype = MPI_CHARACTER;
    b->ndims = ndims;

    IF (PRESENT(rank_write)) h->rank_master = rank_write;

    n1 = LEN(array);
    n2 = SIZE(array);
    b->dims(1) = n1;
    b->dims(2) = n2;

    // Write header

    errcode += write_array_meta(h, id, name);

    h->current_location = b->data_location;

    if (h->rank == h->rank_master) {
        errcode += MPI_FILE_SEEK(h->filehandle, h->current_location, MPI_SEEK_SET,
                      errcode);

        // Actual array
        var_len = INT(b->nelements);
        errcode += MPI_FILE_WRITE(h->filehandle, array, var_len, b->mpitype,
                      MPI_STATUS_IGNORE, errcode);
    }

    h->rank_master = h->default_rank;
    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

}

#endif
