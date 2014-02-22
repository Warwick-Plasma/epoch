#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include "sdf.h"
#include "sdf_list_type.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

int metadata, contents, debug, single, use_mmap, ignore_summary, ascii_header;
int exclude_variables, derived, index_offset;
int64_t array_ndims, *array_starts, *array_ends, *array_strides;
int slice_direction, slice_dim[3];
char *output_file;
struct id_list {
  char *id;
  struct id_list *next;
} *variable_ids, *variable_last_id;

int nrange;

struct range_type {
  unsigned int start, end;
} *range_list;

struct slice_list {
    char *data;
    int datatype, nelements, sz;
    struct slice_list *next;
} *slice_head, *slice_tail;


int close_files(sdf_file_t *h);


void usage(int err)
{
    fprintf(stderr, "usage: sdffilter [options] <sdf_filename>\n");
    fprintf(stderr, "\noptions:\n\
  -h --help            Show this usage message\n\
  -n --no-metadata     Don't show metadata blocks (shown by default)\n\
  -c --contents        Show block's data content\n\
  -s --single          Convert block data to single precision\n\
  -v --variable=id     Find the block with id matching 'id'\n\
  -x --exclude=id      Exclude the block with id matching 'id'\n\
  -m --mmap            Use mmap'ed file I/O\n\
  -i --no-summary      Ignore the metadata summary\n\
  -a --array-section=s Read in the specified array section. The array section\n\
                       's' mimics Python's slicing notation.\n\
  -d --derived         Add derived blocks\n\
  -I --c-indexing      Array indexing starts from 1 by default. If this flag\n\
                       is used then the indexing starts from 0.\n\
  -1 --1dslice=arg     Output 1D slice as a multi-column gnuplot file.\n\
                       The argument is 1,2 or 3 integers separated by commas.\n\
  -H --no-ascii-header When writing multi-column ascii data, a header is\n\
                       included for use by gnuplot or other plotting\n\
                       utilities. This flag disables the header.\n\
");
/*
  -o --output          Output filename\n\
  -D --debug           Show the contents of the debug buffer\n\
*/

    exit(err);
}


int range_sort(const void *v1, const void *v2)
{
    struct range_type *a = (struct range_type *)v1;
    struct range_type *b = (struct range_type *)v2;

    return (a->start - b->start);
}


void parse_1d_slice(char *slice)
{
    int i, len = strlen(slice);
    int done_direction, done_dim1, dim;
    char *old, *ptr;

    done_direction = done_dim1 = 0;
    slice_direction = 0;

    for (i = 0, old = ptr = slice; i < len+1; i++, ptr++) {
        if (*ptr == ',' || *ptr == '\0') {
            if (done_dim1) {
                dim = strtol(old, NULL, 10);
                if (dim > 0) dim -= index_offset;
                if (slice_direction == 2)
                    slice_dim[1] = dim;
                else
                    slice_dim[2] = dim;
            } else if (done_direction) {
                dim = strtol(old, NULL, 10);
                if (dim > 0) dim -= index_offset;
                if (slice_direction == 0)
                    slice_dim[1] = dim;
                else
                    slice_dim[0] = dim;
                done_dim1 = 1;
            } else {
                slice_direction = strtol(old, NULL, 10) - index_offset;
                if (slice_direction < 0 || slice_direction > 2) {
                    fprintf(stderr, "ERROR: invalid slice direction.\n");
                    exit(1);
                }
                done_direction = 1;
            }
            old = ptr + 1;
        }
    }

    if (array_starts) free(array_starts);
    if (array_ends) free(array_ends);
    if (array_strides) free(array_strides);

    array_ndims = 3;

    array_starts  = calloc(array_ndims, sizeof(*array_starts));
    array_ends    = malloc(array_ndims * sizeof(*array_ends));
    array_strides = malloc(array_ndims * sizeof(*array_strides));

    for (i = 0; i < array_ndims; i++) {
        array_strides[i] = 1;
        if (i == slice_direction) {
            array_starts[i] = 0;
            array_ends[i] = INT64_MAX;
        } else {
            array_starts[i] = slice_dim[i];
            array_ends[i] = array_starts[i] + 1;
        }
    }
}


void parse_array_section(char *array_section)
{
    int ndim, i, len = strlen(array_section), done_start, done_end;
    char *ptr, *old;

    if (array_starts) free(array_starts);
    if (array_ends) free(array_ends);
    if (array_strides) free(array_strides);

    array_ndims = 1;
    for (i = 0, ptr = array_section; i < len; i++, ptr++)
        if (*ptr == ',') array_ndims++;

    array_starts  = calloc(array_ndims, sizeof(*array_starts));
    array_ends    = malloc(array_ndims * sizeof(*array_ends));
    array_strides = malloc(array_ndims * sizeof(*array_strides));

    for (i = 0; i < array_ndims; i++)
        array_strides[i] = 1;

    done_start = done_end = ndim = 0;
    for (i = 0, old = ptr = array_section; i < len+1; i++, ptr++) {
        if (*ptr == ':' || *ptr == ',' || *ptr == '\0') {
            if (done_end) {
                array_strides[ndim] = strtol(old, NULL, 10);
                if (array_strides[ndim] == 0) array_strides[ndim] = 1;
                if (array_strides[ndim] < 0) {
                    fprintf(stderr, "ERROR: negative stride values not"
                                    " supported.\n");
                    exit(1);
                }
            } else if (done_start) {
                if (ptr - old > 0) {
                    array_ends[ndim] = strtol(old, NULL, 10);
                    if (array_ends[ndim] > 0) array_ends[ndim] -= index_offset;
                } else
                    array_ends[ndim] = INT64_MAX;
                done_end = 1;
            } else {
                array_starts[ndim] = strtol(old, NULL, 10);
                if (array_starts[ndim] > 0) array_starts[ndim] -= index_offset;
                array_ends[ndim] = array_starts[ndim] + 1;
                done_start = 1;
            }
            old = ptr + 1;
            if (*ptr == ',') {
                done_start = done_end = 0;
                ndim++;
            }
        }
    }
}


char *parse_args(int *argc, char ***argv)
{
    char *ptr, *file = NULL;
    int c, i, err, range, sz, nrange_max, got_include, got_exclude;
    struct range_type *range_tmp;
    struct stat statbuf;
    static struct option longopts[] = {
        { "1dslice",       required_argument, NULL, '1' },
        { "array-section", required_argument, NULL, 'a' },
        { "contents",      no_argument,       NULL, 'c' },
        { "derived",       no_argument,       NULL, 'd' },
        { "help",          no_argument,       NULL, 'h' },
        { "no-ascii-header",no_argument,      NULL, 'H' },
        { "no-summary",    no_argument,       NULL, 'i' },
        { "c-indexing",    no_argument,       NULL, 'I' },
        { "mmap",          no_argument,       NULL, 'm' },
        { "no-metadata",   no_argument,       NULL, 'n' },
        { "single",        no_argument,       NULL, 's' },
        { "variable",      required_argument, NULL, 'v' },
        { "exclude",       required_argument, NULL, 'x' },
        { NULL,            0,                 NULL,  0  }
        //{ "debug",         no_argument,       NULL, 'D' },
        //{ "output",        required_argument, NULL, 'o' },
    };

    metadata = debug = index_offset = 1;
    ascii_header = 1;
    contents = single = use_mmap = ignore_summary = exclude_variables = 0;
    derived = 0;
    slice_direction = -1;
    variable_ids = NULL;
    variable_last_id = NULL;
    output_file = NULL;
    array_starts = array_ends = array_strides = NULL;
    array_ndims = nrange_max = nrange = 0;
    sz = sizeof(struct range_type);

    got_include = got_exclude = 0;

    while ((c = getopt_long(*argc, *argv,
            "1:a:cdhHiImnsv:x:", longopts, NULL)) != -1) {
        switch (c) {
        case '1':
            contents = 1;
            parse_1d_slice(optarg);
            break;
        case 'a':
            contents = 1;
            parse_array_section(optarg);
            break;
        case 'c':
            contents = 1;
            break;
        case 'd':
            derived = 1;
            break;
        case 'h':
            usage(0);
            break;
        case 'H':
            ascii_header = 0;
            break;
        case 'i':
            ignore_summary = 1;
            break;
        case 'I':
            index_offset = 0;
            break;
        case 'm':
            use_mmap = 1;
            break;
        case 'n':
            metadata = 0;
            break;
        case 'o':
            if (output_file) free(output_file);
            output_file = malloc(strlen(optarg)+1);
            memcpy(output_file, optarg, strlen(optarg)+1);
            break;
        case 's':
            single = 1;
            break;
        case 'v':
        case 'x':
            err = 0;
            if (c == 'v') {
                if (got_exclude) err = 1;
                got_include = 1;
            } else {
                if (got_include) err = 1;
                got_exclude = 1;
                exclude_variables = 1;
            }
            if (err) {
                fprintf(stderr, "ERROR: cannot both include and "
                        "exclude variables.\n");
                exit(1);
            }
            if (*optarg >= '0' && *optarg <= '9') {
                ptr = optarg;
                range = 0;
                while (ptr < optarg + strlen(optarg) + 1) {
                    if (range) {
                        i = (int)strtol(ptr, &ptr, 10);
                        if (i == 0)
                            range_list[nrange-1].end = -1;
                        else if (i < range_list[nrange-1].start)
                            nrange--;
                        else
                            range_list[nrange-1].end = i;
                        range = 0;
                    } else {
                        nrange++;
                        // Grow array if necessary
                        if (nrange > nrange_max) {
                            if (nrange_max == 0) {
                                nrange_max = 128;
                                range_list = malloc(nrange_max * sz);
                            } else {
                                i = 2 * nrange_max;

                                range_tmp = malloc(i * sz);
                                memcpy(range_tmp, range_list, nrange_max * sz);
                                free(range_list);
                                range_list = range_tmp;

                                nrange_max = i;
                            }
                        }

                        i = (int)strtol(ptr, &ptr, 10);
                        range_list[nrange-1].start = i;
                        range_list[nrange-1].end = i;
                        if (*ptr == '-') range = 1;
                    }

                    ptr++;
                }
            } else {
                if (!variable_ids) {
                    variable_last_id =
                            variable_ids = malloc(sizeof(*variable_ids));
                } else {
                    variable_last_id->next = malloc(sizeof(*variable_ids));
                    variable_last_id = variable_last_id->next;
                }
                variable_last_id->next = NULL;
                variable_last_id->id = malloc(strlen(optarg)+1);
                memcpy(variable_last_id->id, optarg, strlen(optarg)+1);
            }
            break;
        default:
            usage(1);
        }
    }

    if ((optind+1) == *argc) {
        file = (*argv)[optind];
        err = lstat(file, &statbuf);
        if (err) {
            fprintf(stderr, "Error opening file %s\n", file);
            exit(1);
        }
    } else {
        fprintf(stderr, "No file specified\n");
        usage(1);
    }

    if (nrange > 0) {
        // Sanitize range list
        qsort(range_list, nrange, sz, &range_sort);
        for (i=1; i < nrange; ) {
            if (range_list[i].start <= range_list[i-1].end+1) {
                if (range_list[i].end > range_list[i-1].end)
                    range_list[i-1].end = range_list[i].end;
                memcpy(range_list+i, range_list+i+1, (nrange-i) * sz);
                nrange--;
            } else
                i++;
        }

        // Shrink array
        range_tmp = malloc(nrange * sz);
        memcpy(range_tmp, range_list, nrange * sz);
        free(range_list);
        range_list = range_tmp;
    }

    return file;
}



static int pretty_print_slice(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *cur;
    static sdf_block_t *mesh = NULL;
    int *idx, *fac, dim[3];
    int i, n, rem, sz, print, errcode = 0;
    char *ptr, *dptr;
    float r4;
    double r8;

    if (b->blocktype != SDF_BLOCKTYPE_PLAIN_VARIABLE &&
            b->blocktype != SDF_BLOCKTYPE_PLAIN_DERIVED) return errcode;

    if (slice_direction >= b->ndims) return 0;

    if (!mesh) {
        mesh = sdf_find_block_by_id(h, b->mesh_id);
        if (!mesh) {
            fprintf(stderr, "ERROR: unable to find mesh.\n");
            exit(1);
        }

        cur = h->current_block;
        h->current_block = mesh;

        sdf_read_data(h);

        h->current_block = cur;
        if (!mesh->grids) {
            fprintf(stderr, "ERROR: unable to read mesh.\n");
            exit(1);
        }

        sz = SDF_TYPE_SIZES[mesh->datatype_out];

        slice_head = slice_tail = calloc(1, sizeof(*slice_head));
        slice_tail->datatype = mesh->datatype_out;
        slice_tail->sz = SDF_TYPE_SIZES[slice_tail->datatype];
        slice_tail->nelements = b->local_dims[slice_direction];
        dptr = slice_tail->data = calloc(slice_tail->nelements, sz);

        ptr = mesh->grids[slice_direction];

        // Center grid if necessary
        if (slice_tail->nelements == (mesh->local_dims[slice_direction]-1)) {
            if (slice_tail->datatype == SDF_DATATYPE_REAL4) {
                for (n = 0; n < slice_tail->nelements; n++) {
                    r4 = 0.5 * (*((float*)ptr) + *((float*)ptr+1));
                    memcpy(dptr, &r4, sz);
                    dptr += sz;
                    ptr += sz;
                }
            } else if (slice_tail->datatype == SDF_DATATYPE_REAL8) {
                for (n = 0; n < slice_tail->nelements; n++) {
                    r8 = 0.5 * (*((double*)ptr) + *((double*)ptr+1));
                    memcpy(dptr, &r8, sz);
                    dptr += sz;
                    ptr += sz;
                }
            }
        } else {
            for (n = 0; n < slice_tail->nelements; n++) {
                memcpy(dptr, ptr, sz);
                dptr += sz;
                ptr += sz;
            }
        }

        if (ascii_header) {
            printf("# 1D array slice through (");
            for (n = 0; n < 3; n++) {
                if (n != 0) printf(",");
                if (n == slice_direction)
                    printf(":");
                else
                    printf("%i", slice_dim[n]+index_offset);
            }
            printf(")\n#\n# %s\t%s\t(%s)\n", mesh->id,
                   mesh->dim_labels[slice_direction],
                   mesh->dim_units[slice_direction]);
        }
    }

    idx = malloc(b->ndims * sizeof(*idx));
    fac = malloc(b->ndims * sizeof(*fac));

    rem = 1;
    for (i = 0; i < b->ndims; i++) {
        fac[i] = rem;
        rem *= b->local_dims[i];
    }

    for (i = 0; i < b->ndims; i++) {
        if (slice_dim[i] >= b->local_dims[i]) {
            errcode = 1;
            fprintf(stderr, "ERROR: slice dimension lies outside array.\n");
            goto cleanup;
        } 
        dim[i] = slice_dim[i];
        if (dim[i] < 0) dim[i] += b->local_dims[i];
    }

    sz = SDF_TYPE_SIZES[b->datatype_out];

    slice_tail = slice_tail->next = calloc(1, sizeof(*slice_tail));
    slice_tail->datatype = b->datatype_out;
    slice_tail->sz = SDF_TYPE_SIZES[slice_tail->datatype];
    slice_tail->nelements = b->local_dims[slice_direction];
    dptr = slice_tail->data = calloc(slice_tail->nelements, sz);

    ptr = b->data;
    for (n = 0; n < b->nelements_local; n++) {
        rem = n;
        for (i = b->ndims-1; i >= 0; i--) {
            idx[i] = rem / fac[i];
            rem -= idx[i] * fac[i];
        }

        print = 1;
        for (i = 0; i < b->ndims; i++) {
            if (i == slice_direction) continue;
            if (idx[i] == dim[i]) continue;
            print = 0;
        }

        if (print) {
            memcpy(dptr, ptr, sz);
            dptr += sz;
        }
        ptr += sz;
    }

    if (ascii_header)
        printf("# %s\t%s\t(%s)\n", b->id, b->name, b->units);

cleanup:
    if (idx) free(idx);
    if (fac) free(fac);

    return errcode;
}


static void pretty_print_slice_finish(void)
{
    int i;
    struct slice_list *sl;
    char *ptr;

    if (!slice_head) return;

    if (ascii_header) printf("#\n");

    for (i = 0; i < slice_head->nelements; i++) {
        sl = slice_head;
        while (sl) {
            if (sl != slice_head) printf("    ");
            ptr = sl->data + i * sl->sz;
            switch (sl->datatype) {
            case SDF_DATATYPE_INTEGER4:
                printf("%i", *((uint32_t*)ptr));
                break;
            case SDF_DATATYPE_INTEGER8:
                printf("%llu", *((uint64_t*)ptr));
                break;
            case SDF_DATATYPE_REAL4:
                printf("%14.6E", *((float*)ptr));
                break;
            case SDF_DATATYPE_REAL8:
                printf("%14.6E", *((double*)ptr));
                break;
            }
            sl = sl->next;
        }
        printf("\n");
    }
}


static void pretty_print(sdf_file_t *h, sdf_block_t *b, int idnum)
{
    int *idx, *fac, *printed, *starts = NULL, *ends = NULL;
    int i, n, rem, sz, left, digit, idx0, min_ndims, print;
    char *ptr;
    static const int fmtlen = 32;
    char **fmt;

    if (slice_direction != -1) {
        pretty_print_slice(h, b);
        return;
    }

    idx = malloc(b->ndims * sizeof(*idx));
    fac = malloc(b->ndims * sizeof(*fac));
    fmt = malloc(b->ndims * sizeof(*fmt));
    printed = malloc(b->ndims * sizeof(*printed));

    min_ndims = MIN(array_ndims, b->ndims);

    if (min_ndims > 0) {
        starts = malloc(array_ndims * sizeof(*starts));
        ends   = malloc(array_ndims * sizeof(*ends));
    }

    rem = 1;
    for (i = 0; i < b->ndims; i++) {
        left = b->local_dims[i];
        fac[i] = rem;
        rem *= left;
        digit = 0;
        while (left) {
            left /= 10;
            digit++;
        }
        if (!digit) digit = 1;
        fmt[i] = malloc(fmtlen * sizeof(**fmt));
        if (i == 0)
            snprintf(fmt[i], fmtlen, "%i %%%i.%ii", idnum, digit, digit);
        else
            snprintf(fmt[i], fmtlen, ",%%%i.%ii", digit, digit);
    }

    for (i = 0; i < min_ndims; i++) {
        starts[i] = array_starts[i];
        ends[i] = array_ends[i];
        if (starts[i] < 0) {
            starts[i] += b->local_dims[i];
            if (ends[i] == 0) ends[i] = b->local_dims[i];
        }
        if (ends[i] < 0) ends[i] += b->local_dims[i];
        printed[i] = -array_strides[i];
    }

    sz = SDF_TYPE_SIZES[b->datatype_out];

    ptr = b->data;
    for (n = 0; n < b->nelements_local; n++) {
        rem = n;
        for (i = b->ndims-1; i >= 0; i--) {
            idx0 = idx[i] = rem / fac[i];
            if (i < array_ndims) idx[i] += starts[i];
            rem -= idx0 * fac[i];
        }

        print = 1;
        for (i = 0; i < min_ndims; i++) {
            if (idx[i] < starts[i]) {
                print = 0;
                break;
            }
            if (idx[i] >= ends[i]) {
                print = 0;
                break;
            }
            if ((idx[i] - printed[i]) > 0 &&
                (idx[i] - printed[i]) < array_strides[i]) {
                print = 0;
                break;
            }
        }
        if (print) {
            for (i = 0; i < b->ndims; i++) {
                printf(fmt[i], idx[i]+index_offset);
                printed[i] = idx[i];
            }
            switch (b->datatype_out) {
            case SDF_DATATYPE_INTEGER4:
                printf(":  %i\n", *((uint32_t*)ptr));
                break;
            case SDF_DATATYPE_INTEGER8:
                printf(":  %llu\n", *((uint64_t*)ptr));
                break;
            case SDF_DATATYPE_REAL4:
                printf(":  %12.6E\n", *((float*)ptr));
                break;
            case SDF_DATATYPE_REAL8:
                printf(":  %12.6E\n", *((double*)ptr));
                break;
            }
        }
        ptr += sz;
    }

    free(idx);
    free(fac);
    free(printed);
    if (starts) free(starts);
    if (ends) free(ends);

    for (i = 0; i < b->ndims; i++) free(fmt[i]);
    free(fmt);
}


int main(int argc, char **argv)
{
    char *file = NULL;
    int i, n, block, err, found, idx, len, range_start;
    int nelements_max;
    sdf_file_t *h, *oh;
    sdf_block_t *b, *next, *mesh, *mesh0;
    list_t *station_blocks;
    comm_t comm;

    file = parse_args(&argc, &argv);

#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
#else
    comm = 0;
#endif

    h = sdf_open(file, comm, SDF_READ, use_mmap);
    if (!h) {
        fprintf(stderr, "Error opening file %s\n", file);
        return 1;
    }
    h->use_float = single;
    h->print = debug;
    if (ignore_summary) h->use_summary = 0;

    sdf_read_header(h);
    h->current_block = NULL;

    // If nblocks is negative then the file is corrupt
    if (h->nblocks < 0) {
        block = (-h->nblocks) / 64;
        err = -h->nblocks - 64 * block;
        fprintf(stderr, "Error code %s found at block %i\n",
                sdf_error_codes_c[err], block);
        if (output_file) free(output_file);
        output_file = NULL;
        //return 1;
    }

    if (output_file) {
        oh = sdf_open(output_file, comm, SDF_WRITE, 0);
        sdf_close(oh);
    }

    if (!metadata && !contents) return close_files(h);

    if (derived)
        sdf_read_blocklist_all(h);
    else
        sdf_read_blocklist(h);

    list_init(&station_blocks);

    nelements_max = 0;
    range_start = 0;
    mesh0 = NULL;
    found = 1;
    next = h->blocklist;
    for (i = 0, idx = 1; next; i++, idx++) {
        h->current_block = b = next;
        next = b->next;

        if (nrange > 0 || variable_ids) found = 0;

        for (n = range_start; n < nrange; n++) {
            if (idx < range_list[n].start)
                break;
            if (idx <= range_list[n].end) {
                found = 1;
                break;
            }
            range_start++;
        }

        if (found == 0 && variable_ids) {
            variable_last_id = variable_ids;
            while (variable_last_id) {
                if (!memcmp(b->id, variable_last_id->id,
                        strlen(variable_last_id->id)+1)) {
                    found = 1;
                    break;
                }
                variable_last_id = variable_last_id->next;
            }
        }

        if (exclude_variables) {
            if (found) continue;
        } else {
            if (!found) continue;
        }

        if (metadata && slice_direction == -1) {
            printf("%4i id: %s\n", i+1, b->id);
        }

        if (!contents) continue;

        switch (b->blocktype) {
        case SDF_BLOCKTYPE_PLAIN_DERIVED:
            if (b->station_id) {
                sdf_read_data(h);
                mesh = sdf_find_block_by_id(h, b->mesh_id);
                if (!mesh) continue;
                if (mesh->nelements > nelements_max) {
                    nelements_max = mesh->nelements;
                    mesh0 = mesh;
                }

                if (!mesh->done_data) {
                    h->current_block = mesh;
                    sdf_read_data(h);
                }

                list_append(station_blocks, b);
            }
        case SDF_BLOCKTYPE_PLAIN_VARIABLE:
            sdf_read_data(h);
            pretty_print(h, b, idx);
            break;
        }
    }

    pretty_print_slice_finish();

    if (mesh0 && (variable_ids || nrange > 0)) {
        if (ascii_header) {
            printf("# Stations Time History File\n#\n");
            // This gives garbage output
            //printf("# %s\t%s\t(%s)\n", mesh0->id, mesh0->name, mesh0->units);
            printf("# time\tTime\t(%s)\n", mesh0->units);
        }

        nelements_max = 0;
        b = list_start(station_blocks);
        for (i = 0; i < station_blocks->count; i++) {
            len = strlen(b->station_id);
            if (ascii_header)
                printf("# %s\t%s\t(%s)\n", &b->id[len+1],
                       &b->name[len+1], b->units);
            idx = b->opt + b->nelements - mesh0->opt;
            if (idx > nelements_max)
                nelements_max = idx;
            b = list_next(station_blocks);
        }

        if (ascii_header) printf("#\n");

        for (n = 0; n < nelements_max; n++) {
            printf("%12.6g", ((float*)mesh0->data)[n]);

            b = list_start(station_blocks);
            for (i = 0; i < station_blocks->count; i++) {
                idx = n + mesh0->opt - b->opt;
                if (idx >= 0 && idx < b->nelements)
                    printf("    %12.6g", ((float*)b->data)[idx]);
                else
                    printf("    %12.6g", 0.0);
                b = list_next(station_blocks);
            }

            printf("\n");
        }
    }

    list_destroy(&station_blocks);
    if (range_list) free(range_list);

    return close_files(h);
}


int close_files(sdf_file_t *h)
{
    sdf_close(h);
#ifdef PARALLEL
    MPI_Finalize();
#endif

    return 0;
}
