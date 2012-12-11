#include "sdf.h"
#include <stdlib.h>

//#define SDF_COMMON_MESH_LENGTH (4 + 8 + SDF_ID_LENGTH + 4 * b->ndims)

#define SDF_COMMON_MESH_INFO() do { \
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


int sdf_point_factor(sdf_file_t *h, int *local_npoints);


int sdf_read_point_mesh_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;

    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - npoints   INTEGER(i8)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_ARRAY_REAL8(b->dim_mults, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_labels, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_units, b->ndims);

    SDF_READ_ENTRY_INT4(b->geometry);

    SDF_READ_ENTRY_ARRAY_REAL8(b->extents, 2*b->ndims);

    SDF_READ_ENTRY_INT8(b->npoints);
    for (i = 0; i < b->ndims; i++) b->dims[i] = b->npoints;

    b->stagger = SDF_STAGGER_VERTEX;

    return 0;
}



int sdf_read_point_variable_info(sdf_file_t *h)
{
    sdf_block_t *b;

    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - npoints   INTEGER(i8)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_REAL8(b->mult);

    SDF_READ_ENTRY_ID(b->units);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_INT8(b->npoints);
    b->dims[0] = b->npoints;

    return 0;
}



static int sdf_create_1d_distribution(sdf_file_t *h, int64_t global, int local,
        int start)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;
    int lengths[3];
    MPI_Aint disp[3];
    MPI_Datatype types[3];

    lengths[0] = 1;
    lengths[1] = local;
    lengths[2] = 1;
    disp[0] = 0;
    disp[1] = start * SDF_TYPE_SIZES[b->datatype];
    disp[2] = (int)global * SDF_TYPE_SIZES[b->datatype];
    types[0] = MPI_LB;
    types[1] = b->mpitype;
    types[2] = MPI_UB;

    MPI_Type_create_struct(3, lengths, disp, types, &b->distribution);
    MPI_Type_commit(&b->distribution);
#endif
    return 0;
}



static int sdf_free_distribution(sdf_file_t *h)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;

    MPI_Type_free(&b->distribution);
#endif
    return 0;
}



static int sdf_helper_read_array(sdf_file_t *h, void **var_in, int count)
{
    sdf_block_t *b = h->current_block;
    char **var = (char **)var_in;

    if (h->mmap) {
        *var = h->mmap + h->current_location;
        return 0;
    }

    if (*var) free(*var);
    *var = malloc(count * SDF_TYPE_SIZES[b->datatype]);
#ifdef PARALLEL
    MPI_File_set_view(h->filehandle, h->current_location, b->mpitype,
            b->distribution, "native", MPI_INFO_NULL);
    MPI_File_read_all(h->filehandle, *var, count, b->mpitype,
            MPI_STATUS_IGNORE);
    MPI_File_set_view(h->filehandle, 0, MPI_BYTE, MPI_BYTE, "native",
            MPI_INFO_NULL);
#else
    fseeko(h->filehandle, h->current_location, SEEK_SET);
    fread(*var, SDF_TYPE_SIZES[b->datatype], count, h->filehandle);
#endif

    return 0;
}



int sdf_read_point_mesh(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_point_mesh_info(h);

    sdf_point_factor(h, &b->nlocal);

    h->current_location = b->data_location;

    if (!b->grids) b->grids = calloc(b->ndims, sizeof(float*));

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n<b->ndims; n++) SDF_DPRNT("%i ",b->nlocal);
        SDF_DPRNT("\n");
        h->indent = 2;
    }
    for (n = 0; n < 3; n++) {
        if (b->ndims > n) {
#ifdef PARALLEL
            sdf_create_1d_distribution(h, b->npoints, b->nlocal,
                    b->starts[0]);
#endif
            sdf_helper_read_array(h, &b->grids[n], b->nlocal);
            sdf_free_distribution(h);
            sdf_convert_array_to_float(h, &b->grids[n], b->nlocal);
            if (h->use_random) sdf_randomize_array(h, &b->grids[n], b->nlocal);
            if (h->print) {
                SDF_DPRNT("%s: ", b->dim_labels[n]);
                SDF_DPRNTar(b->grids[n], b->nlocal);
            }
            h->current_location = h->current_location
                    + SDF_TYPE_SIZES[b->datatype] * b->npoints;
        }
    }

    b->done_data = 1;

    return 0;
}



int sdf_point_factor(sdf_file_t *h, int *local_npoints)
{
    sdf_block_t *b = h->current_block;
#ifdef PARALLEL
    int npoint_min, split_big;

    npoint_min = (int)b->npoints / h->ncpus;
    split_big = (int)b->npoints - h->ncpus * npoint_min;

    if (h->rank >= split_big) {
        b->starts[0] = split_big * (npoint_min + 1)
            + (h->rank - split_big) * npoint_min;
        *local_npoints = npoint_min;
    } else {
        b->starts[0] = h->rank * (npoint_min + 1);
        *local_npoints = npoint_min + 1;
    }
#else
    *local_npoints = (int)b->npoints;
#endif

    return 0;
}



int sdf_read_point_variable(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_blocklist(h);

    sdf_point_factor(h, &b->nlocal);

    h->current_location = b->data_location;

#ifdef PARALLEL
    sdf_create_1d_distribution(h, b->npoints, b->nlocal, b->starts[0]);
#endif
    sdf_helper_read_array(h, &b->data, b->nlocal);
    sdf_free_distribution(h);
    sdf_convert_array_to_float(h, &b->data, b->nlocal);
    if (h->use_random) sdf_randomize_array(h, &b->data, b->nlocal);
    h->current_location = h->current_location
            + SDF_TYPE_SIZES[b->datatype] * b->npoints;

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n<b->ndims; n++) SDF_DPRNT("%i ",b->nlocal);
        SDF_DPRNT("\n  ");
        SDF_DPRNTar(b->data, b->nlocal);
    }

    b->done_data = 1;

    return 0;
}
