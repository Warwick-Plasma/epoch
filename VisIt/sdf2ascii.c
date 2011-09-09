#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "sdf_common.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
    struct stat statbuf;
    int i, j, err;
    sdf_file_t *h;
    sdf_block_t *b;
    int rank = 0, size = 1;
    comm_t comm;

    err = 1; 
    if (argc == 2) err = lstat(argv[1], &statbuf);

    if (err) {
        fprintf(stderr, "usage: sdf2ascii <sdf_filename>\n");
        return 1;
    }
#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
#endif

    h = sdf_open(argv[1], rank, comm, 1);
    if (!h) {
        fprintf(stderr, "Error opening file %s\n", argv[1]);
        return 1;
    }
    h->use_float = 0;
    h->print = 1;
    sdf_set_ncpus(h, size);

    sdf_read_blocklist(h);
#ifdef SDF_DEBUG
    printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    b = h->current_block = h->blocklist;
    for (i=0; i<h->nblocks; i++) {
        b = h->current_block;
        sdf_read_data(h);
        h->current_block = b->next;
#ifdef SDF_DEBUG
        printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    }
#ifdef SDF_DEBUG
    printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0';
#endif
    sdf_close(h);
    printf("\n");

#ifdef PARALLEL
    MPI_Finalize();
#endif

    return 0;
}
