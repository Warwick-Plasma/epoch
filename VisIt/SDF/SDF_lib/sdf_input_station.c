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


int sdf_read_station_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int nstations, nvariables, step0, step_inc;
    double time0, time_inc;
    char use_mult[4];
    char *use_mult_ptr = use_mult;

    // Metadata is
    // - nentries  INTEGER(i4)
    // - type_size INTEGER(i4)
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
    SDF_READ_ENTRY_INT4(nstations);
    SDF_READ_ENTRY_INT4(nvariables);
    SDF_READ_ENTRY_INT4(step0);
    SDF_READ_ENTRY_INT4(step_inc);
    SDF_READ_ENTRY_REAL8(time0);
    SDF_READ_ENTRY_REAL8(time_inc);
    SDF_READ_ENTRY_STRINGLEN(use_mult_ptr,4);

    SDF_READ_ENTRY_ARRAY_ID(b->station_ids, nstations);
    SDF_READ_ENTRY_ARRAY_STRING(b->station_names, nstations);
    SDF_READ_ENTRY_ARRAY_INT4(b->station_nvars, nstations);
    SDF_READ_ENTRY_ARRAY_INT4(b->station_move, nstations);
    SDF_READ_ENTRY_ARRAY_REAL8(b->station_x, nstations);
    if (b->ndims > 1)
        SDF_READ_ENTRY_ARRAY_REAL8(b->station_y, nstations);
    if (b->ndims > 2)
        SDF_READ_ENTRY_ARRAY_REAL8(b->station_z, nstations);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, nvariables);
    SDF_READ_ENTRY_ARRAY_STRING(b->variable_names, nvariables);
    SDF_READ_ENTRY_ARRAY_INT4(b->variable_types, nvariables);
    SDF_READ_ENTRY_ARRAY_ID(b->variable_units, nvariables);
    if (use_mult[0])
        SDF_READ_ENTRY_ARRAY_REAL8(b->variable_mults, nvariables);

    b->stagger = SDF_STAGGER_VERTEX;

    return 0;
}
