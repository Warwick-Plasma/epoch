#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include "sdf.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

int metadata, contents, debug, single, mmap;
char *variable_id;


void usage(int err)
{
    fprintf(stderr, "usage: sdf2ascii [options] <sdf_filename>\n");
    fprintf(stderr, "\noptions:\n\
  -h --help           Show this usage message\n\
  -n --no-metadata    Don't show metadata blocks (shown by default)\n\
  -c --contents       Show block's data content\n\
  -s --single         Convert block data to single precision\n\
  -v --variable=id    Find the block with id matching 'id'\n\
  -m --mmap           Use mmap'ed file I/O\n\
");
/*
  -d --debug          Show the contents of the debug buffer\n\
*/

    exit(err);
}


char *parse_args(int *argc, char ***argv)
{
    char *file = NULL;
    char c;
    int err;
    struct stat statbuf;
    static struct option longopts[] = {
        { "no-metadata", no_argument,    NULL,      'n' },
        { "contents", no_argument,       NULL,      'c' },
        //{ "debug",    no_argument,       NULL,      'd' },
        { "help",     no_argument,       NULL,      'h' },
        { "variable", required_argument, NULL,      'v' },
        { "mmap",     no_argument,       NULL,      'm' },
        { NULL,       0,                 NULL,       0  }
    };

    metadata = debug = 1;
    contents = single = mmap = 0;
    variable_id = NULL;

    while ((c = getopt_long(*argc, *argv, "hncsmv:", longopts, NULL)) != -1) {
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
            if (variable_id) free(variable_id);
            variable_id = malloc(strlen(optarg)+1);
            memcpy(variable_id, optarg, strlen(optarg)+1);
            break;
        case 'm':
            mmap = 1;
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


    return file;
}


int main(int argc, char **argv)
{
    char *file = NULL;
    int i, buflen, len, block, err;
    sdf_file_t *h;
    sdf_block_t *b;
    int rank = 0, size = 1;
    comm_t comm;

    file = parse_args(&argc, &argv);

#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
#endif

    h = sdf_open(file, rank, comm, mmap);
    if (!h) {
        fprintf(stderr, "Error opening file %s\n", file);
        return 1;
    }
    h->use_float = single;
    h->print = debug;
    sdf_set_ncpus(h, size);

    //sdf_read_blocklist(h);
    // Read blocklist -- just like sdf_read_blocklist() but with the ability
    // to parse as we go
    sdf_read_header(h);
    h->current_block = NULL;

    // If nblocks is negative then the file is corrupt
    if (h->nblocks < 0) {
        block = (-h->nblocks) / 64;
        err = -h->nblocks - 64 * block;
        fprintf(stderr, "Error code %s found at block %i\n",
                sdf_error_codes_c[err], block);
        return 1;
    }

    // Read the whole summary block into a temporary buffer on rank 0
    buflen = h->summary_size;
    h->current_location = h->start_location = h->summary_location;
    h->buffer = malloc(buflen);

    if (h->rank == h->rank_master) {
        sdf_seek(h);
        sdf_read_bytes(h, h->buffer, buflen);
    }

    // Send the temporary buffer to all processors
    sdf_broadcast(h, h->buffer, buflen);

    // Construct the metadata blocklist using the contents of the buffer
    len = 0;
    if (variable_id && metadata) len = strlen(variable_id) + 1;
    for (i = 0; i < h->nblocks; i++) {
        if (len)
            h->dbg = h->dbg_buf; *h->dbg = '\0';
        sdf_read_block_info(h);
        if (len && !memcmp(h->current_block->id, variable_id, len)) {
            printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
            break;
        }
    }

    free(h->buffer);
    h->buffer = NULL;
    h->current_block = h->blocklist;
#ifdef SDF_DEBUG
    if (metadata)
        printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif

    // Done reading metadata

    if (!contents) return 0;

    len = 0;
    if (variable_id) len = strlen(variable_id) + 1;
printf("cont %s %i\n", variable_id, len);
    b = h->current_block = h->blocklist;
    for (i=0; i<h->nblocks; i++) {
        b = h->current_block;
        if (len) {
            h->print = 0;
            if (!memcmp(b->id, variable_id, len)) {
                h->print = 1;
                sdf_read_data(h);
                printf("c %s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
                break;
            }
        } else
            sdf_read_data(h);
        h->current_block = b->next;
#ifdef SDF_DEBUG_ALL
#ifdef SDF_DEBUG
        if (debug)
            printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
#endif
    }
#ifdef SDF_DEBUG_ALL
#ifdef SDF_DEBUG
    if (debug)
        printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
#endif
    sdf_close(h);
    printf("\n");

#ifdef PARALLEL
    MPI_Finalize();
#endif

    return 0;
}
