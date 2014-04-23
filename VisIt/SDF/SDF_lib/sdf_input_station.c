#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sdf.h>

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


int sdf_read_station_info(sdf_file_t *h)
{
    sdf_block_t *b;
    char use_mult[5];
    char *use_mult_ptr = use_mult;

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

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_INT8(b->nelements);
    SDF_READ_ENTRY_INT4(b->type_size);
    SDF_READ_ENTRY_INT4(b->nstations);
    SDF_READ_ENTRY_INT4(b->nvariables);
    SDF_READ_ENTRY_INT4(b->step);
    SDF_READ_ENTRY_INT4(b->step_increment);
    SDF_READ_ENTRY_REAL8(b->time);
    SDF_READ_ENTRY_REAL8(b->time_increment);
    SDF_READ_ENTRY_STRINGLEN(use_mult_ptr,4);

    SDF_READ_ENTRY_ARRAY_ID(b->station_ids, b->nstations);
    SDF_READ_ENTRY_ARRAY_STRING(b->station_names, b->nstations);
    SDF_READ_ENTRY_ARRAY_INT4(b->station_nvars, b->nstations);
    SDF_READ_ENTRY_ARRAY_INT4(b->station_move, b->nstations);
    SDF_READ_ENTRY_ARRAY_REAL8(b->station_x, b->nstations);
    if (b->ndims > 1)
        SDF_READ_ENTRY_ARRAY_REAL8(b->station_y, b->nstations);
    if (b->ndims > 2)
        SDF_READ_ENTRY_ARRAY_REAL8(b->station_z, b->nstations);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->nvariables);
    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->nvariables);
    SDF_READ_ENTRY_ARRAY_INT4(b->variable_types, b->nvariables);
    SDF_READ_ENTRY_ARRAY_ID(b->dim_units, b->nvariables);
    if (use_mult[0])
        SDF_READ_ENTRY_ARRAY_REAL8(b->dim_mults, b->nvariables);

    b->stagger = SDF_STAGGER_VERTEX;

    return 0;
}
