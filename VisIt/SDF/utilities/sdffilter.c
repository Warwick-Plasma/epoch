#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include "sdf.h"
#include "sdf_list_type.h"
#include "sdf_helper.h"
#include "stack_allocator.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

int metadata, contents, debug, single, use_mmap, ignore_summary, ascii_header;
int exclude_variables, derived, index_offset, element_count;
int just_id, verbose_metadata;
int format_rowindex, format_index, format_number;
int64_t array_ndims, *array_starts, *array_ends, *array_strides;
int slice_direction, slice_dim[3];
char *output_file;
char *format_float, *format_int, *format_space;
//static char *default_float = "%3.2fE%d";
static char *default_float = "%14.6E";
static char *default_int   = "%" PRIi64;
static char *default_space = "    ";
static char *default_indent = "  ";
static char indent[64];

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

static char width_fmt[16];
#define SET_WIDTH(string) do { \
        int _l = strlen((string)); \
        snprintf(width_fmt, 16, "%%-%is", _l); \
    } while(0)

#define PRINT(name,variable,fmt) do { \
        if (!(variable)) break; \
        printf(indent, 1); \
        printf(width_fmt, (name)); \
        printf(" " fmt, (variable)); \
        printf("\n"); \
    } while(0)

#define PRINTAR(name,array,fmt,len) do { \
        int _i; \
        if (!(array)) break; \
        printf(indent, 1); \
        printf(width_fmt, (name)); \
        printf(" (" fmt, (array)[0]); \
        for (_i = 1; _i < (len); _i++) printf("," fmt, (array)[_i]); \
        printf(")\n"); \
    } while(0)

#define PRINTDAR(name,array,fmt,len) do { \
        int _i; \
        if (!(array)) break; \
        printf(indent, 1); \
        printf(width_fmt, (name)); \
        printf(" (" fmt, sdf_datatype_c[(array)[0]]); \
        for (_i = 1; _i < (len); _i++) \
            printf("," fmt, sdf_datatype_c[(array)[_i]]); \
        printf(")\n"); \
    } while(0)


int close_files(sdf_file_t *h);


void usage(int err)
{
    fprintf(stderr, "usage: sdffilter [options] <sdf_filename>\n");
    fprintf(stderr, "\noptions:\n\
  -h --help            Show this usage message\n\
  -n --no-metadata     Don't show metadata blocks (shown by default)\n\
  -j --just-id         Only show ID and number for metadata blocks\n\
  -l --less-verbose    Print metadata less verbosely\n\
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
  -C --count=n         When pretty-printing array contents, write 'n'\n\
                       elements per line.\n\
  -F --format-float=f  Use specified format for printing floating-point array\n\
                       contents.\n\
  -N --format-int=f    Use specified format for printing integer array\n\
                       contents.\n\
  -S --format-space=f  Use specified spacing between array elements.\n\
  -K --format-number   Show block number before each row of array elements.\n\
  -R --format-rowindex Show array indices before each row of array elements.\n\
  -J --format-index    Show array indices before each array element.\n\
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
        { "count",         required_argument, NULL, 'C' },
        { "derived",       no_argument,       NULL, 'd' },
        { "help",          no_argument,       NULL, 'h' },
        { "no-ascii-header",no_argument,      NULL, 'H' },
        { "format-float",  required_argument, NULL, 'F' },
        { "no-summary",    no_argument,       NULL, 'i' },
        { "c-indexing",    no_argument,       NULL, 'I' },
        { "just-id",       no_argument,       NULL, 'j' },
        { "format-index",  no_argument,       NULL, 'J' },
        { "format-number", no_argument,       NULL, 'K' },
        { "less-verbose",  no_argument,       NULL, 'l' },
        { "mmap",          no_argument,       NULL, 'm' },
        { "no-metadata",   no_argument,       NULL, 'n' },
        { "format-int",    required_argument, NULL, 'N' },
        { "format-rowindex",no_argument,      NULL, 'R' },
        { "single",        no_argument,       NULL, 's' },
        { "format-space",  required_argument, NULL, 'S' },
        { "variable",      required_argument, NULL, 'v' },
        { "exclude",       required_argument, NULL, 'x' },
        { NULL,            0,                 NULL,  0  }
        //{ "debug",         no_argument,       NULL, 'D' },
        //{ "output",        required_argument, NULL, 'o' },
    };

    metadata = debug = index_offset = element_count = verbose_metadata = 1;
    ascii_header = 1;
    contents = single = use_mmap = ignore_summary = exclude_variables = 0;
    derived = format_rowindex = format_index = format_number = just_id = 0;
    slice_direction = -1;
    variable_ids = NULL;
    variable_last_id = NULL;
    output_file = NULL;
    array_starts = array_ends = array_strides = NULL;
    array_ndims = nrange_max = nrange = 0;
    sz = sizeof(struct range_type);

    format_int = malloc(strlen(default_int)+1);
    memcpy(format_int, default_int, strlen(default_int)+1);

    format_float = malloc(strlen(default_float)+1);
    memcpy(format_float, default_float, strlen(default_float)+1);

    format_space = malloc(strlen(default_space)+1);
    memcpy(format_space, default_space, strlen(default_space)+1);

    got_include = got_exclude = 0;

    while ((c = getopt_long(*argc, *argv,
            "1:a:cC:dF:hHiIjJKlmnN:RsS:v:x:", longopts, NULL)) != -1) {
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
        case 'C':
            element_count = strtol(optarg, NULL, 10);
            if (element_count < 1) element_count = 1;
            break;
        case 'd':
            derived = 1;
            break;
        case 'F':
            free(format_float);
            format_float = malloc(strlen(optarg)+1);
            memcpy(format_float, optarg, strlen(optarg)+1);
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
        case 'j':
            just_id = 1;
            break;
        case 'J':
            format_index = 1;
            break;
        case 'K':
            format_number = 1;
            break;
        case 'l':
            verbose_metadata = 0;
            break;
        case 'm':
            use_mmap = 1;
            break;
        case 'n':
            metadata = 0;
            break;
        case 'N':
            free(format_int);
            format_int = malloc(strlen(optarg)+1);
            memcpy(format_int, optarg, strlen(optarg)+1);
            break;
        case 'o':
            if (output_file) free(output_file);
            output_file = malloc(strlen(optarg)+1);
            memcpy(output_file, optarg, strlen(optarg)+1);
            break;
        case 'R':
            format_rowindex = 1;
            break;
        case 's':
            single = 1;
            break;
        case 'S':
            free(format_space);
            format_space = malloc(strlen(optarg)+1);
            memcpy(format_space, optarg, strlen(optarg)+1);
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

    special_format = 0;
    if (format_float[strlen(format_float)-1] == 'd')
        special_format = 1;

    return file;
}


static int set_array_section(sdf_block_t *b)
{
    return sdf_block_set_array_section(b, array_ndims, array_starts,
                                       array_ends, array_strides);
}


static int pretty_print_slice(sdf_file_t *h, sdf_block_t *b)
{
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

        set_array_section(mesh);
        sdf_helper_read_data(h, mesh);

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
    int exponent, special_format = 0;
    struct slice_list *sl;
    char *ptr;
    double r8;

    if (!slice_head) return;

    if (ascii_header) printf("#\n");

    special_format = 0;
    if (format_float[strlen(format_float)-1] == 'd')
        special_format = 1;

    for (i = 0; i < slice_head->nelements; i++) {
        sl = slice_head;
        while (sl) {
            if (sl != slice_head) printf(format_space,1);
            ptr = sl->data + i * sl->sz;
            switch (sl->datatype) {
            case SDF_DATATYPE_INTEGER4:
                printf(format_int, *((uint32_t*)ptr));
                break;
            case SDF_DATATYPE_INTEGER8:
                printf(format_int, *((uint64_t*)ptr));
                break;
            case SDF_DATATYPE_REAL4:
                if (special_format) {
                    r8 = *((float*)ptr);
                    if (r8 == 0)
                        exponent = 0;
                    else {
                        exponent = (int)floor(log10(fabs(r8)));
                        r8 *= pow(10, -1.0 * exponent);
                    }
                    printf(format_float, r8, exponent);
                } else
                    printf(format_float, *((float*)ptr));
                break;
            case SDF_DATATYPE_REAL8:
                if (special_format) {
                    r8 = *((double*)ptr);
                    if (r8 == 0)
                        exponent = 0;
                    else {
                        exponent = (int)floor(log10(fabs(r8)));
                        r8 *= pow(10, -1.0 * exponent);
                    }
                    printf(format_float, r8, exponent);
                } else
                    printf(format_float, *((double*)ptr));
                break;
            }
            sl = sl->next;
        }
        printf("\n");
    }
}


static void pretty_print_mesh(sdf_file_t *h, sdf_block_t *b, int idnum)
{
    int *idx, *fac;
    int i, n, rem, sz, left, digit, ncount, dim;
    int exponent, special_format = 0;
    char *ptr;
    static const int fmtlen = 32;
    char **fmt;
    double r8;

    if (slice_direction != -1) {
        pretty_print_slice(h, b);
        return;
    }

    idx = malloc(b->ndims * sizeof(*idx));
    fac = malloc(b->ndims * sizeof(*fac));
    fmt = malloc(b->ndims * sizeof(*fmt));

    rem = 1;
    for (i = 0; i < b->ndims; i++) {
        if (b->array_starts)
            left = b->array_ends[i] - b->array_starts[i];
        else
            left = b->local_dims[i];
        fac[i] = rem;
        rem *= left;
        digit = 0;
        if (b->array_ends)
            left = b->array_ends[i] + index_offset - 1;
        while (left) {
            left /= 10;
            digit++;
        }
        if (!digit) digit = 1;
        if (format_rowindex || format_index) {
            ptr = fmt[i] = malloc(fmtlen * sizeof(**fmt));
            if (i != 0) *ptr++ = ',';
            sz = snprintf(ptr, fmtlen-2, "%%%i.%ii", digit, digit);
            if (i == b->ndims-1) {
                ptr += sz;
                *ptr++ = ')';
                *ptr++ = '\0';
            }
        } else
            fmt[i] = calloc(1, sizeof(**fmt));
    }

    special_format = 0;
    if (format_float[strlen(format_float)-1] == 'd')
        special_format = 1;

    sz = SDF_TYPE_SIZES[b->datatype_out];

    ncount = 0;
    dim = 0;
    if (b->array_starts) {
        for (i = 0; i < b->ndims; i++)
            idx[i] = b->array_starts[i];
        for (i = 0; i < b->ndims; i++) {
            if (b->array_ends[i] > b->array_starts[i])
                break;
            dim++;
        }
    }
    ptr = b->grids[dim];

    for (n = 0; n < b->nelements_local; n++) {
        ncount++;
        if (ncount == 1) {
            if (format_number) printf("%i ", idnum);
        } else {
            if (format_index) printf(" ");
        }

        if ((ncount ==1 && format_rowindex) || format_index) {
            for (i = 0; i < b->ndims; i++) {
                if (i == dim) {
                    printf(fmt[i], idx[i]+index_offset);
                } else {
                    if (i != 0) printf(",");
                    printf("0");
                    if (i == b->ndims-1) printf(")");
                }
            }
        }

        if (ncount != 1 && format_index) printf(format_space,1);

        switch (b->datatype_out) {
        case SDF_DATATYPE_INTEGER4:
            printf(format_int, *((uint32_t*)ptr));
            break;
        case SDF_DATATYPE_INTEGER8:
            printf(format_int, *((uint64_t*)ptr));
            break;
        case SDF_DATATYPE_REAL4:
            if (special_format) {
                r8 = *((float*)ptr);
                if (r8 == 0)
                    exponent = 0;
                else {
                    exponent = (int)floor(log10(fabs(r8)));
                    r8 *= pow(10, -1.0 * exponent);
                }
                printf(format_float, r8, exponent);
            } else
                printf(format_float, *((float*)ptr));
            break;
        case SDF_DATATYPE_REAL8:
            if (special_format) {
                r8 = *((double*)ptr);
                if (r8 == 0)
                    exponent = 0;
                else {
                    exponent = (int)floor(log10(fabs(r8)));
                    r8 *= pow(10, -1.0 * exponent);
                }
                printf(format_float, r8, exponent);
            } else
                printf(format_float, *((double*)ptr));
            break;
        }
        idx[dim]++;
        ptr += sz;
        if (b->array_ends && idx[dim] >= b->array_ends[dim]) {
            idx[dim] = 0;
            dim++;
            ptr = b->grids[dim];
            ncount = element_count;
        } else if (idx[dim] >= b->local_dims[dim]) {
            idx[dim] = 0;
            dim++;
            ptr = b->grids[dim];
            ncount = element_count;
        }
        if (ncount == element_count) {
            printf("\n");
            ncount = 0;
        }
    }
    if (ncount) printf("\n");

    free(idx);
    free(fac);

    for (i = 0; i < b->ndims; i++) free(fmt[i]);
    free(fmt);
}


static void pretty_print(sdf_file_t *h, sdf_block_t *b, int idnum)
{
    int *idx, *fac;
    int i, n, rem, sz, left, digit, ncount, idx0;
    int exponent, special_format = 0;
    char *ptr;
    static const int fmtlen = 32;
    char **fmt;
    double r8;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH ||
            b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
        pretty_print_mesh(h, b, idnum);
        return;
    }

    if (slice_direction != -1) {
        pretty_print_slice(h, b);
        return;
    }

    idx = malloc(b->ndims * sizeof(*idx));
    fac = malloc(b->ndims * sizeof(*fac));
    fmt = malloc(b->ndims * sizeof(*fmt));

    rem = 1;
    for (i = 0; i < b->ndims; i++) {
        if (b->array_starts)
            left = b->array_ends[i] - b->array_starts[i];
        else
            left = b->local_dims[i];
        fac[i] = rem;
        rem *= left;
        digit = 0;
        if (b->array_ends)
            left = b->array_ends[i] + index_offset - 1;
        while (left) {
            left /= 10;
            digit++;
        }
        if (!digit) digit = 1;
        if (format_rowindex || format_index) {
            fmt[i] = malloc(fmtlen * sizeof(**fmt));
            if (i == 0)
                snprintf(fmt[i], fmtlen, "%%%i.%ii", digit, digit);
            else if (i == b->ndims-1)
                snprintf(fmt[i], fmtlen, ",%%%i.%ii)", digit, digit);
            else
                snprintf(fmt[i], fmtlen, ",%%%i.%ii", digit, digit);
        } else
            fmt[i] = calloc(1, sizeof(**fmt));
    }

    special_format = 0;
    if (format_float[strlen(format_float)-1] == 'd')
        special_format = 1;

    sz = SDF_TYPE_SIZES[b->datatype_out];

    ptr = b->data;
    ncount = 0;
    for (n = 0; n < b->nelements_local; n++) {
        rem = n;
        for (i = b->ndims-1; i >= 0; i--) {
            idx0 = idx[i] = rem / fac[i];
            if (b->array_starts) idx[i] += b->array_starts[i];
            rem -= idx0 * fac[i];
        }

        ncount++;
        if (ncount == 1) {
            if (format_number)
                printf("%i ", idnum);
            for (i = 0; i < b->ndims; i++)
                printf(fmt[i], idx[i]+index_offset);
        } else {
            if (format_index) {
                printf(" ");
                for (i = 0; i < b->ndims; i++)
                    printf(fmt[i], idx[i]+index_offset);
            }
            printf(format_space,1);
        }

        switch (b->datatype_out) {
        case SDF_DATATYPE_INTEGER4:
            printf(format_int, *((uint32_t*)ptr));
            break;
        case SDF_DATATYPE_INTEGER8:
            printf(format_int, *((uint64_t*)ptr));
            break;
        case SDF_DATATYPE_REAL4:
            if (special_format) {
                r8 = *((float*)ptr);
                if (r8 == 0)
                    exponent = 0;
                else {
                    exponent = (int)floor(log10(fabs(r8)));
                    r8 *= pow(10, -1.0 * exponent);
                }
                printf(format_float, r8, exponent);
            } else
                printf(format_float, *((float*)ptr));
            break;
        case SDF_DATATYPE_REAL8:
            if (special_format) {
                r8 = *((double*)ptr);
                if (r8 == 0)
                    exponent = 0;
                else {
                    exponent = (int)floor(log10(fabs(r8)));
                    r8 *= pow(10, -1.0 * exponent);
                }
                printf(format_float, r8, exponent);
            } else
                printf(format_float, *((double*)ptr));
            break;
        }
        if (ncount == element_count) {
            printf("\n");
            ncount = 0;
        }
        ptr += sz;
    }
    if (ncount) printf("\n");

    free(idx);
    free(fac);

    for (i = 0; i < b->ndims; i++) free(fmt[i]);
    free(fmt);
}


static void print_metadata_plain_mesh(sdf_block_t *b)
{
    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SET_WIDTH("dim_labels:");
    PRINTAR("dim_mults:", b->dim_mults, "%g", b->ndims);
    PRINTAR("dim_labels:", b->dim_labels, "%s", b->ndims);
    PRINTAR("dim_units:", b->dim_units, "%s", b->ndims);
    PRINT("geometry:", sdf_geometry_c[b->geometry], "%s");
    PRINTAR("extents:", b->extents, "%g", 2*b->ndims);
    PRINTAR("dims:", b->dims, "%" PRIi64, b->ndims);
}


static void print_metadata_point_mesh(sdf_block_t *b)
{
    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - npoints   INTEGER(i8)
    // - speciesid CHARACTER(id_length)

    SET_WIDTH("dim_labels:");
    PRINTAR("dim_mults", b->dim_mults, "%g", b->ndims);
    PRINTAR("dim_labels", b->dim_labels, "%s", b->ndims);
    PRINTAR("dim_units", b->dim_units, "%s", b->ndims);
    PRINT("geometry:", sdf_geometry_c[b->geometry], "%s");
    PRINTAR("extents", b->extents, "%g", 2*b->ndims);
    //PRINTAR("dims", b->dims, "%" PRIi64, b->ndims);
    PRINT("nelements:", b->nelements, "%" PRIi64);
    if (b->material_id)
        PRINT("species id:", b->material_id, "%s");
}


static void print_metadata_plain_variable(sdf_block_t *b)
{
    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - dims      INTEGER(i4), DIMENSION(ndims)
    // - stagger   INTEGER(i4)

    SET_WIDTH("mesh id:");
    PRINT("mult:", b->mult, "%g");
    PRINT("units:", b->units, "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINTAR("dims:", b->dims, "%" PRIi64, b->ndims);
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
}


static void print_metadata_point_variable(sdf_block_t *b)
{
    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - npoints   INTEGER(i8)
    // - speciesid CHARACTER(id_length)

    SET_WIDTH("species id:");
    PRINT("mult:", b->mult, "%g");
    PRINT("units:", b->units, "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINT("nelements:", b->nelements, "%" PRIi64);
    if (b->material_id)
        PRINT("species id:", b->material_id, "%s");
}


static void print_metadata_constant(sdf_block_t *b)
{
    int32_t i4;
    int64_t i8;
    float r4;
    double r8;

    // Metadata is
    // - value     TYPE_SIZE

    printf("%svalue: ", indent);

    switch (b->datatype) {
    case SDF_DATATYPE_INTEGER4:
        memcpy(&i4, b->const_value, sizeof(i4));
        printf("%i", i4);
        break;
    case SDF_DATATYPE_INTEGER8:
        memcpy(&i8, b->const_value, sizeof(i8));
        printf("%" PRIi64, i8);
        break;
    case SDF_DATATYPE_REAL4:
        memcpy(&r4, b->const_value, sizeof(r4));
        printf("%g", r4);
        break;
    case SDF_DATATYPE_REAL8:
        memcpy(&r8, b->const_value, sizeof(r8));
        printf("%g", r8);
        break;
    //case SDF_DATATYPE_REAL16:
    //    printf("%g", (double)b->const_value);
    //    break;
    case SDF_DATATYPE_CHARACTER:
        printf("%c", *b->const_value);
        break;
    case SDF_DATATYPE_LOGICAL:
        if (*b->const_value)
            printf("True");
        else
            printf("False");
        break;
    }

    printf("\n");
}


static void print_metadata_array(sdf_block_t *b)
{
    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SET_WIDTH("dims:");
    PRINTAR("dims:", b->dims, "%" PRIi64, b->ndims);
}


static void print_metadata_cpu_split(sdf_block_t *b)
{
    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SET_WIDTH("geometry:");
    PRINT("geometry:", sdf_geometry_c[b->geometry], "%s");
    PRINTAR("dims:", b->dims, "%" PRIi64, b->ndims);
}


static void print_metadata_run(sdf_block_t *b)
{
    struct run_info *run = b->data;
    time_t time;
    char *stime;
    char version[32];

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

    SET_WIDTH("compile_machine:");
    snprintf(version, 32, "%i.%i.%i",
             run->version, run->revision, run->minor_rev);
    stime = version;
    PRINT("version:", stime, "%s");
    PRINT("commit id:", run->commit_id, "%s");
    PRINT("sha1sum:", run->sha1sum, "%s");
    PRINT("compile_machine:", run->compile_machine, "%s");
    PRINT("compile_flags:", run->compile_flags, "%s");
    PRINT("defines:", run->defines, "%" PRIi64);
    time = run->compile_date; stime = ctime(&time); stime[strlen(stime)-1]='\0';
    PRINT("compile_date:", stime, "%s");
    time = run->run_date; stime = ctime(&time); stime[strlen(stime)-1] = '\0';
    PRINT("run_date:", stime, "%s");
    time = run->io_date; stime = ctime(&time); stime[strlen(stime)-1] = '\0';
    PRINT("io_date:", stime, "%s");
}


static void print_metadata_stitched(sdf_block_t *b)
{
    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SET_WIDTH("variable ids:");
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINTAR("variable ids:", b->variable_ids, "%s", b->ndims);
}


static void print_metadata_stitched_material(sdf_block_t *b)
{
    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matnames  ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SET_WIDTH("material names:");
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINTAR("material names:", b->material_names, "%s", b->ndims);
    PRINTAR("variable ids:", b->variable_ids, "%s", b->ndims);
}


static void print_metadata_stitched_matvar(sdf_block_t *b)
{
    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SET_WIDTH("variable ids:");
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINT("material id:", b->material_id, "%s");
    PRINTAR("variable ids:", b->variable_ids, "%s", b->ndims);
}


static void print_metadata_stitched_species(sdf_block_t *b)
{
    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - matname   CHARACTER(string_length)
    // - specnames ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SET_WIDTH("species names:");
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
    PRINT("mesh id:", b->mesh_id, "%s");
    PRINT("material id:", b->material_id, "%s");
    PRINT("material name:", b->material_name, "%s");
    PRINTAR("species names:", b->material_names, "%s", b->ndims);
    PRINTAR("variable ids:", b->variable_ids, "%s", b->ndims);
}


static void print_metadata_stitched_obstacle_group(sdf_block_t *b)
{
    // Metadata is
    // - stagger         INTEGER(i4)
    // - obstacle_id     CHARACTER(id_length)
    // - vfm_id          CHARACTER(id_length)
    // - obstacle_names  ndims*CHARACTER(string_length)

    SET_WIDTH("volume fraction id:");
    PRINT("stagger:", sdf_stagger_c[b->stagger], "%s");
    PRINT("obstacle id:", b->obstacle_id, "%s");
    PRINT("volume fraction id:", b->vfm_id, "%s");
    PRINTAR("obstacle names:", b->material_names, "%s", b->ndims);
}


static void print_metadata_station(sdf_block_t *b)
{
    // Metadata is
    // - nelements INTEGER(i8)
    // - entry_len INTEGER(i4)
    // - nstations INTEGER(i4)
    // - nvars     INTEGER(i4)
    // - step0     INTEGER(i4)
    // - step_inc  INTEGER(i4)
    // - time0     REAL(r8)
    // - time_inc  REAL(r8)
    // - use_mult  CHARACTER(1)
    // - padding   CHARACTER(3)
    // - statids   CHARACTER(id_length), DIMENSION(nstations)
    // - statnames CHARACTER(string_length), DIMENSION(nstations)
    // - statnvars INTEGER(i4), DIMENSION(nstations)
    // - statmove  INTEGER(i4), DIMENSION(nstations)
    // - statx0    REAL(r8), DIMENSION(nstations*ndims)
    // - varids    CHARACTER(id_length), DIMENSION(nvars)
    // - varnames  CHARACTER(string_length), DIMENSION(nvars)
    // - vartypes  INTEGER(i4), DIMENSION(nvars)
    // - varunits  CHARACTER(id_length), DIMENSION(nvars)
    // - varmults  REAL(r8), DIMENSION(use_mult*nvars)

    SET_WIDTH("time_increment:");
    PRINT("nelements:", b->nelements, "%" PRIi64);
    PRINT("entry_len:", b->type_size, "%i");
    PRINT("nstations:", b->nstations, "%i");
    PRINT("nvariables:", b->nvariables, "%i");
    PRINT("step0:", b->step, "%i");
    PRINT("step_increment:", b->step_increment, "%i");
    PRINT("time0:", b->time, "%g");
    PRINT("time_increment:", b->time_increment, "%g");
    PRINTAR("station_ids:", b->station_ids, "%s", b->nstations);
    PRINTAR("station_names:", b->station_names, "%s", b->nstations);
    PRINTAR("station_nvars:", b->station_nvars, "%i", b->nstations);
    PRINTAR("station_move:", b->station_move, "%i", b->nstations);
    PRINTAR("station_x:", b->station_x, "%g", b->nstations);
    if (b->ndims > 1)
        PRINTAR("station_y:", b->station_y, "%g", b->nstations);
    if (b->ndims > 2)
        PRINTAR("station_z:", b->station_z, "%g", b->nstations);
    PRINTAR("variable_ids:", b->variable_ids, "%s", b->nvariables);
    PRINTAR("variable_names:", b->material_names, "%s", b->nvariables);
    PRINTDAR("variable_types:", b->variable_types, "%s", b->nvariables);
    PRINTAR("variable_units:", b->dim_units, "%s", b->nvariables);
    if (b->dim_mults)
        PRINTAR("variable_mults:", b->dim_mults, "%g", b->nvariables);
}


static void print_metadata(sdf_block_t *b, int inum, int nblocks)
{
    int digit = 0;
    static const int fmtlen = 32;
    char fmt[fmtlen];

    while (nblocks) {
        nblocks /= 10;
        digit++;
    }

    snprintf(fmt, fmtlen, "Block %%%ii, ID: %%s", digit);
    printf(fmt, inum, b->id);
    if (!b->in_file)
        printf("  (derived)");
    printf("\n");
    if (just_id) return;

    sprintf(indent, default_indent, 1);

    if (verbose_metadata)
        SET_WIDTH("block_location:");
    else
        SET_WIDTH("blocktype:");

    PRINT("name:", b->name, "%s");
    PRINT("blocktype:", sdf_blocktype_c[b->blocktype], "%s");
    PRINT("datatype:", sdf_datatype_c[b->datatype], "%s");

    if (verbose_metadata) {
        PRINT("ndims:", b->ndims, "%i");
        PRINT("data_length:", b->data_length, "%" PRIi64);
        PRINT("info_length:", b->info_length, "%i");
        PRINT("data_location:", b->data_location, "%" PRIi64);
        PRINT("block_location:", b->block_start, "%" PRIi64);
        PRINT("next_block:", b->next_block_location, "%" PRIi64);
    }

    sprintf(indent+strlen(indent), default_indent, 1);

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_PLAIN_MESH:
    case SDF_BLOCKTYPE_LAGRANGIAN_MESH:
        print_metadata_plain_mesh(b);
        break;
    case SDF_BLOCKTYPE_POINT_MESH:
        print_metadata_point_mesh(b);
        break;
    case SDF_BLOCKTYPE_PLAIN_VARIABLE:
    case SDF_BLOCKTYPE_PLAIN_DERIVED:
        print_metadata_plain_variable(b);
        break;
    case SDF_BLOCKTYPE_POINT_VARIABLE:
    case SDF_BLOCKTYPE_POINT_DERIVED:
        print_metadata_point_variable(b);
        break;
    case SDF_BLOCKTYPE_CONSTANT:
        print_metadata_constant(b);
        break;
    case SDF_BLOCKTYPE_ARRAY:
        print_metadata_array(b);
        break;
    case SDF_BLOCKTYPE_CPU_SPLIT:
        print_metadata_cpu_split(b);
        break;
    case SDF_BLOCKTYPE_RUN_INFO:
        print_metadata_run(b);
        break;
    case SDF_BLOCKTYPE_STITCHED:
    case SDF_BLOCKTYPE_CONTIGUOUS:
    case SDF_BLOCKTYPE_STITCHED_TENSOR:
    case SDF_BLOCKTYPE_CONTIGUOUS_TENSOR:
        print_metadata_stitched(b);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL:
        print_metadata_stitched_material(b);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATVAR:
        print_metadata_stitched_matvar(b);
        break;
    case SDF_BLOCKTYPE_STITCHED_SPECIES:
    case SDF_BLOCKTYPE_CONTIGUOUS_SPECIES:
        print_metadata_stitched_species(b);
        break;
    case SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP:
        print_metadata_stitched_obstacle_group(b);
        break;
    case SDF_BLOCKTYPE_STATION:
    case SDF_BLOCKTYPE_STATION_DERIVED:
        print_metadata_station(b);
        break;
    }

    printf("\n");
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
    stack_init();

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

        if (metadata && slice_direction == -1)
            print_metadata(b, idx, h->nblocks);

        if (!contents) continue;

        switch (b->blocktype) {
        case SDF_BLOCKTYPE_PLAIN_DERIVED:
            set_array_section(b);
            sdf_helper_read_data(h, b);
            if (b->station_id) {
                mesh = sdf_find_block_by_id(h, b->mesh_id);
                if (!mesh) continue;
                if (mesh->nelements > nelements_max) {
                    nelements_max = mesh->nelements;
                    mesh0 = mesh;
                }

                if (!mesh->done_data)
                    sdf_helper_read_data(h, mesh);

                list_append(station_blocks, b);
            } else
                pretty_print(h, b, idx);
            break;
        case SDF_BLOCKTYPE_PLAIN_VARIABLE:
        case SDF_BLOCKTYPE_PLAIN_MESH:
        case SDF_BLOCKTYPE_POINT_VARIABLE:
        case SDF_BLOCKTYPE_POINT_MESH:
            set_array_section(b);
            sdf_helper_read_data(h, b);
            pretty_print(h, b, idx);
            break;
        default:
            printf("Unsupported blocktype %s\n",
                   sdf_blocktype_c[b->blocktype]);
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
