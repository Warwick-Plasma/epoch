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


int sdf_read_plain_mesh_info(sdf_file_t *h)
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
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_ARRAY_REAL8(b->dim_mults, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_labels, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_units, b->ndims);

    SDF_READ_ENTRY_INT4(b->geometry);

    SDF_READ_ENTRY_ARRAY_REAL8(b->extents, 2*b->ndims);

    SDF_READ_ENTRY_ARRAY_INT4(b->dims_in, b->ndims);
    for (i = 0; i < b->ndims; i++) b->dims[i] = b->dims_in[i];

    b->stagger = SDF_STAGGER_VERTEX;
    for (i = 0; i < b->ndims; i++) b->const_value[i] = 1;

    return 0;
}



int sdf_read_plain_variable_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
#ifdef PARALLEL
    int local_start[SDF_MAXDIMS];
#endif

    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - dims      INTEGER(i4), DIMENSION(ndims)
    // - stagger   INTEGER(i4)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_REAL8(b->mult);

    SDF_READ_ENTRY_ID(b->units);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_INT4(b->dims_in, b->ndims);
    for (i = 0; i < b->ndims; i++) b->dims[i] = b->dims_in[i];

    SDF_READ_ENTRY_INT4(b->stagger);
    for (i = 0; i < b->ndims; i++) b->const_value[i] = (b->stagger & 1<<i);

#ifdef PARALLEL
    // Calculate per block parallel factorisation
    // This will be fixed up later once we have the whole block list.
    sdf_factor(h, local_start);
#else
    b->nlocal = 1;
    for (i=0; i < b->ndims; i++) {
        b->local_dims[i] = b->dims[i];
        b->nlocal *= b->dims[i];
    }
    for (i=b->ndims; i < 3; i++) b->local_dims[i] = 1;
#endif

    return 0;
}



static int sdf_create_1d_distribution(sdf_file_t *h, int global, int local,
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
    disp[1] = start * b->type_size;
    disp[2] = global * b->type_size;
    types[0] = MPI_LB;
    types[1] = b->mpitype;
    types[2] = MPI_UB;

    MPI_Type_create_struct(3, lengths, disp, types, &b->distribution);
    MPI_Type_commit(&b->distribution);
#endif
    return 0;
}



static int sdf_plain_mesh_distribution(sdf_file_t *h)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;
    int n;
    int local_start[SDF_MAXDIMS], sizes[SDF_MAXDIMS];

    if (b->ng) return 0;

    for (n=0; n < b->ndims; n++) sizes[n] = b->dims[n];

    // Get local_start for creating subarray
    sdf_factor(h, local_start);

    MPI_Type_create_subarray(b->ndims, sizes, b->local_dims, local_start,
        MPI_ORDER_FORTRAN, b->mpitype, &b->distribution);
    MPI_Type_commit(&b->distribution);
#endif

    return 0;
}



static int sdf_free_distribution(sdf_file_t *h)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;

    if (b->ng) return 0;

    MPI_Type_free(&b->distribution);
#endif
    return 0;
}



static int sdf_helper_read_array_halo(sdf_file_t *h, void **var_in)
{
    sdf_block_t *b = h->current_block;
    char **var_ptr = (char**)var_in;
    char *var = *var_ptr;
    char *vptr;
    char convert;
    int count, subsizes[SDF_MAXDIMS], starts[SDF_MAXDIMS];
    uint64_t offset;
#ifdef PARALLEL
    MPI_Datatype distribution, facetype;
    int face[SDF_MAXDIMS];
    int i, tag;
    char *p1, *p2;
#else
    int i, n;
    int idx, rem, orem, dim, npt;
#endif

    b->nlocal = 1;
    for (i=0; i < b->ndims; i++) {
        subsizes[i] = b->local_dims[i];
        starts[i] = b->ng;
        b->local_dims[i] += 2 * b->ng;
        b->dims[i] += 2 * b->ng;
        b->nlocal *= b->local_dims[i];
    }
    for (i=b->ndims; i < SDF_MAXDIMS; i++)
        subsizes[i] = 1;
    count = b->nlocal;

    if (*var_ptr) free(*var_ptr);
    *var_ptr = var = malloc(count * b->type_size);

#ifdef PARALLEL
    MPI_Type_create_subarray(b->ndims, b->local_dims, subsizes, starts,
        MPI_ORDER_FORTRAN, b->mpitype, &distribution);

    MPI_Type_commit(&distribution);

    MPI_File_set_view(h->filehandle, h->current_location, b->mpitype,
            b->distribution, "native", MPI_INFO_NULL);
    MPI_File_read_all(h->filehandle, var, 1, distribution,
            MPI_STATUS_IGNORE);
    MPI_File_set_view(h->filehandle, 0, MPI_BYTE, MPI_BYTE, "native",
            MPI_INFO_NULL);

    MPI_Type_free(&distribution);

    // Swap ghostcell faces
    for (i=0; i < b->ndims; i++) {
        face[i] = b->local_dims[i] - 2 * b->ng;
        starts[i] = b->ng;
    }

    tag = 1;
    offset = b->type_size;
    for (i=0; i < b->ndims; i++) {
        face[i] = b->ng;
        starts[i] = 0;

        MPI_Type_create_subarray(b->ndims, b->local_dims, face, starts,
            MPI_ORDER_FORTRAN, b->mpitype, &facetype);
        MPI_Type_commit(&facetype);

        p1 = b->data + b->ng * offset;
        p2 = b->data + (b->local_dims[i] - b->ng) * offset;
        MPI_Sendrecv(p1, 1, facetype, b->proc_min[i], tag, p2, 1, facetype,
            b->proc_max[i], tag, h->comm, MPI_STATUS_IGNORE);
        tag++;

        p1 = b->data + (b->local_dims[i] - 2 * b->ng) * offset;
        p2 = b->data;
        MPI_Sendrecv(p1, 1, facetype, b->proc_max[i], tag, p2, 1, facetype,
            b->proc_min[i], tag, h->comm, MPI_STATUS_IGNORE);
        tag++;

        MPI_Type_free(&facetype);

        face[i] = b->local_dims[i];
        offset *= b->local_dims[i];
    }
#else
    idx = 0;
    offset = b->ng;
    for (i=0; i < b->ndims; i++) {
        idx += offset;
        offset *= b->local_dims[i];
    }

    offset = b->local_dims[0];
    npt = b->local_dims[0] - 2 * b->ng;
    for (i=1; i < b->ndims; i++) {
        npt += (b->local_dims[i] - 2 * b->ng) * offset;
        offset *= b->local_dims[i];
    }

    vptr = var + idx * b->type_size;

    while (idx < npt) {
        fseeko(h->filehandle, h->current_location, SEEK_SET);
        fread(vptr, b->type_size, subsizes[0], h->filehandle);
        h->current_location += subsizes[0] * b->type_size;
        vptr += b->local_dims[0] * b->type_size;
        idx += b->local_dims[0];

        dim = b->local_dims[0];
        orem = idx / dim;
        offset = 2 * b->ng * dim;
        for (i=1; i < b->ndims; i++) {
            rem = orem / dim;
            n = orem - rem * dim;
            dim = b->local_dims[i];
            if (n == dim - b->ng) {
                vptr += offset * b->type_size;
                idx += offset;
            }
            offset *= dim;
            orem = rem;
        }
    }
#endif

    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8)
        convert = 1;
    else
        convert = 0;

    if (convert) {
        int i;
        float *r4;
        double *old_var, *r8;
        r8 = old_var = (double*)var;
        r4 = (float*)(*var_ptr);
        if (!r4) {
            *var_ptr = malloc(count * sizeof(float));
            r4 = (float*)(*var_ptr);
        }
        for (i=0; i < count; i++)
            *r4++ = (float)(*r8++);
        if (!h->mmap) free(old_var);
        b->datatype_out = SDF_DATATYPE_REAL4;
        b->type_size_out = 4;
#ifdef PARALLEL
        b->mpitype_out = MPI_FLOAT;
#endif
    }

    return 0;
}



static int sdf_helper_read_array(sdf_file_t *h, void **var_in, int count)
{
    sdf_block_t *b = h->current_block;
    char **var_ptr = (char**)var_in;
    char *var = *var_ptr;
    char convert;

    if (b->ng) return sdf_helper_read_array_halo(h, var_in);

    if (h->mmap) {
        *var_ptr = h->mmap + h->current_location;
        return 0;
    }

    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8) {
        convert = 1;
        *var_ptr = var = malloc(count * b->type_size);
    } else {
        convert = 0;
        if (!var) *var_ptr = var = malloc(count * b->type_size);
    }

#ifdef PARALLEL
    MPI_File_set_view(h->filehandle, h->current_location, b->mpitype,
            b->distribution, "native", MPI_INFO_NULL);
    MPI_File_read_all(h->filehandle, var, count, b->mpitype,
            MPI_STATUS_IGNORE);
    MPI_File_set_view(h->filehandle, 0, MPI_BYTE, MPI_BYTE, "native",
            MPI_INFO_NULL);
#else
    fseeko(h->filehandle, h->current_location, SEEK_SET);
    fread(var, b->type_size, count, h->filehandle);
#endif

    if (convert) {
        int i;
        float *r4;
        double *old_var, *r8;
        r8 = old_var = (double*)var;
        r4 = (float*)(*var_ptr);
        if (!r4) {
            *var_ptr = malloc(count * sizeof(float));
            r4 = (float*)(*var_ptr);
        }
        for (i=0; i < count; i++)
            *r4++ = (float)(*r8++);
        if (!h->mmap) free(old_var);
        b->datatype_out = SDF_DATATYPE_REAL4;
        b->type_size_out = 4;
#ifdef PARALLEL
        b->mpitype_out = MPI_FLOAT;
#endif
    }

    return 0;
}



int sdf_read_plain_mesh(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int local_start[SDF_MAXDIMS];
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_blocklist(h);

    sdf_factor(h, local_start);

    h->current_location = b->data_location;

    if (!b->grids) b->grids = calloc(3, sizeof(float*));

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n<b->ndims; n++) SDF_DPRNT("%i ",b->local_dims[n]);
        SDF_DPRNT("\n");
        h->indent = 2;
    }
    for (n = 0; n < 3; n++) {
        if (b->ndims > n) {
            sdf_create_1d_distribution(h, b->dims[n], b->local_dims[n],
                    local_start[n]);
            sdf_helper_read_array(h, &b->grids[n], b->local_dims[n]);
            sdf_free_distribution(h);
            if (h->print) {
                SDF_DPRNT("%s: ", b->dim_labels[n]);
                SDF_DPRNTar(b->grids[n], b->local_dims[n]);
            }
            h->current_location = h->current_location
                    + b->type_size * b->dims[n];
        } else {
            b->grids[n] = calloc(1, b->type_size);
        }
    }

    b->done_data = 1;

    return 0;
}


int sdf_read_plain_variable(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_plain_variable_info(h);

    h->current_location = b->data_location;

    sdf_plain_mesh_distribution(h);

    sdf_helper_read_array(h, &b->data, b->nlocal);

    sdf_free_distribution(h);

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n<b->ndims; n++) SDF_DPRNT("%i ",b->local_dims[n]);
        SDF_DPRNT("\n  ");
        SDF_DPRNTar(b->data, b->nlocal);
    }

    b->done_data = 1;

    return 0;
}
