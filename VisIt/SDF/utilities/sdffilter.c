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

int metadata, contents, debug, single, use_mmap, ignore_summary;
int exclude_variables;
char *output_file;
struct id_list {
  char *id;
  struct id_list *next;
} *variable_ids, *variable_last_id;

int nrange;

struct range_type {
  unsigned int start, end;
} *range_list;


int close_files(sdf_file_t *h);

void usage(int err)
{
    fprintf(stderr, "usage: sdffilter [options] <sdf_filename>\n");
    fprintf(stderr, "\noptions:\n\
  -h --help           Show this usage message\n\
  -n --no-metadata    Don't show metadata blocks (shown by default)\n\
  -c --contents       Show block's data content\n\
  -s --single         Convert block data to single precision\n\
  -v --variable=id    Find the block with id matching 'id'\n\
  -x --exclude=id     Exclude the block with id matching 'id'\n\
  -m --mmap           Use mmap'ed file I/O\n\
  -i --no-summary     Ignore the metadata summary\n\
");
/*
  -o --output         Output filename\n\
  -d --debug          Show the contents of the debug buffer\n\
*/

    exit(err);
}


int range_sort(const void *v1, const void *v2)
{
    struct range_type *a = (struct range_type *)v1;
    struct range_type *b = (struct range_type *)v2;

    return (a->start - b->start);
}


char *parse_args(int *argc, char ***argv)
{
    char *ptr, *file = NULL;
    int c, i, err, range, sz, nrange_max, got_include, got_exclude;
    struct range_type *range_tmp;
    struct stat statbuf;
    static struct option longopts[] = {
        { "no-metadata", no_argument,    NULL,      'n' },
        { "contents", no_argument,       NULL,      'c' },
        //{ "debug",    no_argument,       NULL,      'd' },
        { "help",     no_argument,       NULL,      'h' },
        { "variable", required_argument, NULL,      'v' },
        { "exclude",  required_argument, NULL,      'x' },
        { "mmap",     no_argument,       NULL,      'm' },
        { "no-summary", no_argument,     NULL,      'i' },
        //{ "output", required_argument,   NULL,      'o' },
        { NULL,       0,                 NULL,       0  }
    };

    metadata = debug = 1;
    contents = single = use_mmap = ignore_summary = exclude_variables = 0;
    variable_ids = NULL;
    variable_last_id = NULL;
    output_file = NULL;
    nrange_max = nrange = 0;
    sz = sizeof(struct range_type);

    got_include = got_exclude = 0;

    while ((c = getopt_long(*argc, *argv,
            "hncsmiv:x:", longopts, NULL)) != -1) {
        switch (c) {
        case 'h':
            usage(0);
            break;
        case 'n':
            metadata = 0;
            break;
        case 'c':
            contents = 1;
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
            }
            if (c == 'x') {
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
        case 'o':
            if (output_file) free(output_file);
            output_file = malloc(strlen(optarg)+1);
            memcpy(output_file, optarg, strlen(optarg)+1);
            break;
        case 'm':
            use_mmap = 1;
            break;
        case 'i':
            ignore_summary = 1;
            break;
        default:
            usage(1);
        }
    }

    if ((optind+1) == *argc) {
        file = (*argv)[optind];
        err = stat(file, &statbuf);
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

    sdf_read_blocklist_all(h);

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

        if (metadata) {
            printf("%4i id: %s\n", i+1, b->id);
        }

        if (contents) {
            sdf_read_data(h);
            if (b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED && b->station_id) {
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
        }
    }

    if (mesh0 && (variable_ids || nrange > 0)) {
        printf("# Stations Time History File\n");

        nelements_max = 0;
        b = list_start(station_blocks);
        for (i = 0; i < station_blocks->count; i++) {
            len = strlen(b->station_id);
            printf("# %s\t%s\t(%s)\n", b->id, &b->name[len+1], b->units);
            idx = b->opt + b->nelements - mesh0->opt;
            if (idx > nelements_max)
                nelements_max = idx;
            b = list_next(station_blocks);
        }

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
