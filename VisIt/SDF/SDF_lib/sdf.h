#ifndef _SDF_COMMON_H_
#define _SDF_COMMON_H_

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <stdio.h>
#include <inttypes.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef SDF_DEBUG_ALL
#define SDF_DEBUG
#endif

#define SOI4  4
#define SOI8  8
#define SOF4  4
#define SOF8  8

#define SDF_MAXDIMS 4
#define SDF_ID_LENGTH 32
#define SDF_HEADER_LENGTH (11 * SOI4 + 2 * SOI8 + SOF8 + 12 + SDF_ID_LENGTH)
#define SDF_BLOCK_HEADER_LENGTH \
    (4 + 3 * SOI4 + 3 * SOI8 + SDF_ID_LENGTH + h->string_length)
#define SDF_SUMMARY_OFFSET (4 + 3 * SOI4 + SDF_ID_LENGTH + SOI8)
#define SDF_ENDIANNESS 16911887

#define SDF_VERSION  1
#define SDF_REVISION 2
#define SDF_LIB_VERSION  4
#define SDF_LIB_REVISION 0

#define SDF_MAGIC "SDF1"

#ifdef __cplusplus
extern "C" {
#endif

enum sdf_blocktype {
    SDF_BLOCKTYPE_SCRUBBED = -1,
    SDF_BLOCKTYPE_NULL = 0,
    SDF_BLOCKTYPE_PLAIN_MESH,
    SDF_BLOCKTYPE_POINT_MESH,
    SDF_BLOCKTYPE_PLAIN_VARIABLE,
    SDF_BLOCKTYPE_POINT_VARIABLE,
    SDF_BLOCKTYPE_CONSTANT,
    SDF_BLOCKTYPE_ARRAY,
    SDF_BLOCKTYPE_RUN_INFO,
    SDF_BLOCKTYPE_SOURCE,
    SDF_BLOCKTYPE_STITCHED_TENSOR,
    SDF_BLOCKTYPE_STITCHED_MATERIAL,
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
};

enum sdf_geometry {
    SDF_GEOMETRY_NULL = 0,
    SDF_GEOMETRY_CARTESIAN,
    SDF_GEOMETRY_CYLINDRICAL,
    SDF_GEOMETRY_SPHERICAL,
};

enum sdf_stagger {
    SDF_STAGGER_CELL_CENTRE = 0,
    SDF_STAGGER_FACE_X,
    SDF_STAGGER_FACE_Y,
    SDF_STAGGER_EDGE_Z,
    SDF_STAGGER_FACE_Z,
    SDF_STAGGER_EDGE_Y,
    SDF_STAGGER_EDGE_X,
    SDF_STAGGER_VERTEX,
};

enum sdf_datatype {
    SDF_DATATYPE_NULL = 0,
    SDF_DATATYPE_INTEGER4,
    SDF_DATATYPE_INTEGER8,
    SDF_DATATYPE_REAL4,
    SDF_DATATYPE_REAL8,
    SDF_DATATYPE_REAL16,
    SDF_DATATYPE_CHARACTER,
    SDF_DATATYPE_LOGICAL,
    SDF_DATATYPE_OTHER,
};

static const int SDF_TYPE_SIZES[] = {
    0,  // SDF_DATATYPE_NULL = 0,
    4,  // SDF_DATATYPE_INTEGER4,
    8,  // SDF_DATATYPE_INTEGER8,
    4,  // SDF_DATATYPE_REAL4,
    8,  // SDF_DATATYPE_REAL8,
    16, // SDF_DATATYPE_REAL16,
    1,  // SDF_DATATYPE_CHARACTER,
    1,  // SDF_DATATYPE_LOGICAL,
    0,  // SDF_DATATYPE_OTHER,
};

enum sdf_dimension {
    SDF_DIMENSION_IRRELEVANT = 0,
    SDF_DIMENSION_1D,
    SDF_DIMENSION_2D,
    SDF_DIMENSION_3D,
};

enum sdf_error_codes {
    SDF_ERR_SUCCESS = 0,
    SDF_ERR_ACCESS,
    SDF_ERR_AMODE,
    SDF_ERR_BAD_FILE,
    SDF_ERR_CONVERSION,
    SDF_ERR_DUP_DATAREP,
    SDF_ERR_FILE,
    SDF_ERR_FILE_EXISTS,
    SDF_ERR_FILE_IN_USE,
    SDF_ERR_INFO,
    SDF_ERR_INFO_KEY,
    SDF_ERR_INFO_NOKEY,
    SDF_ERR_INFO_VALUE,
    SDF_ERR_IO,
    SDF_ERR_NOT_SAME,
    SDF_ERR_NO_SPACE,
    SDF_ERR_NO_SUCH_FILE,
    SDF_ERR_QUOTA,
    SDF_ERR_READ_ONLY,
    SDF_ERR_UNSUPPORTED_DATAREP,
    SDF_ERR_UNSUPPORTED_OPERATION,
    SDF_ERR_UNKNOWN,
};


#define SDF_READ  1
#define SDF_WRITE 2

extern const char *sdf_blocktype_c[];
extern const char *sdf_geometry_c[];
extern const char *sdf_stagger_c[];
extern const char *sdf_datatype_c[];
extern const char *sdf_error_codes_c[];

extern const int sdf_blocktype_len;
extern const int sdf_geometry_len;
extern const int sdf_stagger_len;
extern const int sdf_datatype_len;
extern const int sdf_error_codes_len;

#ifdef PARALLEL
    typedef MPI_Comm comm_t;
#else
    typedef int comm_t;
#endif

typedef struct sdf_block sdf_block_t;
typedef struct sdf_file sdf_file_t;

struct sdf_block {
    // This struct must be changed with care and the SDF_LIB_VERSION bumped
    // if the resulting struct is not aligned the same.
    double *extents, *dim_mults;
    double *station_x, *station_y, *station_z;
    double mult, time, time_increment;
    uint64_t dims[3];
    uint64_t block_start, next_block_location, data_location;
    uint64_t nelements, data_length, *nelements_blocks, *data_length_blocks;
    uint32_t ndims, geometry, datatype, blocktype, info_length;
    uint32_t type_size, stagger, datatype_out, type_size_out;
    uint32_t nstations, nvariables, step, step_increment;
    uint32_t *dims_in, *station_nvars, *variable_types, *station_index;
    int32_t *station_move;
    int local_dims[3], nm, nelements_local, n_ids, opt, ng, nfaces;
    char const_value[16];
    char *id, *units, *mesh_id, *material_id;
    char *vfm_id, *obstacle_id, *station_id;
    char *name, *material_name, *must_read;
    char **dim_labels, **dim_units;
    char **station_ids, **variable_ids;
    char **station_names, **material_names;
    int *node_list, *boundary_cells;
    void **grids, *data;
    char done_header, done_info, done_data, dont_allocate, dont_display;
    char dont_own_data, use_mult;
    sdf_block_t *next, *prev;
    sdf_block_t *subblock, *subblock2;
    sdf_block_t *(*populate_data)(sdf_file_t *, sdf_block_t *);
#ifdef PARALLEL
    MPI_Datatype mpitype, distribution, mpitype_out;
    int cpu_split[SDF_MAXDIMS], starts[SDF_MAXDIMS];
    int proc_min[3], proc_max[3];
#endif
};

struct sdf_file {
    uint64_t dbg_count;
    uint32_t sdf_lib_version, sdf_lib_revision;
    uint32_t sdf_extension_version, sdf_extension_revision;
    uint32_t file_version, file_revision;
    char *dbg, *dbg_buf, **extension_names;
    // Lines above should never be changed.
    // Lines below must be changed with care and the SDF_LIB_VERSION bumped
    // if the resulting struct is not aligned the same.
    double time;
    uint64_t first_block_location, summary_location, start_location, soi, sof;
    uint64_t current_location;
    uint32_t jobid1, jobid2, endianness, summary_size;
    uint32_t block_header_length, string_length;
    uint32_t code_io_version, step;
    int32_t nblocks, error_code;
    int rank, ncpus, ndomains, rank_master, indent, print;
    char *buffer, *filename;
    char done_header, restart_flag, other_domains, use_float, use_summary;
    char use_random, station_file, swap;
    char *code_name, *error_message;
    sdf_block_t *blocklist, *tail, *current_block;
    char *mmap;
    void *ext_data;
#ifdef PARALLEL
    MPI_File filehandle;
#else
    FILE *filehandle;
#endif
    comm_t comm;
};

struct run_info {
    uint64_t defines;
    uint32_t version, revision, compile_date, run_date, io_date, minor_rev;
    char *commit_id, *sha1sum, *compile_machine, *compile_flags;
};


sdf_file_t *sdf_open(const char *filename, comm_t comm, int mode, int use_mmap);
int sdf_close(sdf_file_t *h);
int sdf_free_blocklist_data(sdf_file_t *h);
sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id);
sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name);
int sdf_read_header(sdf_file_t *h);
int sdf_read_summary(sdf_file_t *h);
int sdf_read_blocklist(sdf_file_t *h);
int sdf_read_blocklist_all(sdf_file_t *h);
int sdf_read_block_info(sdf_file_t *h);
int sdf_read_data(sdf_file_t *h);
int sdf_get_domain_extents(sdf_file_t *h, int rank, int *start, int *local);
int sdf_seek(sdf_file_t *h);
int sdf_seek_set(sdf_file_t *h, off_t offset);
int sdf_read_bytes(sdf_file_t *h, char *buf, int buflen);
int sdf_broadcast(sdf_file_t *h, void *buf, int size);
int sdf_write(sdf_file_t *h);

// internal routines

int sdf_factor(sdf_file_t *h);
int sdf_convert_array_to_float(sdf_file_t *h, void **var_in, int count);
int sdf_randomize_array(sdf_file_t *h, void **var_in, int count);
int sdf_set_rank_master(sdf_file_t *h, int rank);
int sdf_read_nblocks(sdf_file_t *h);

int sdf_abort(sdf_file_t *h);
int sdf_read_next_block_header(sdf_file_t *h);
int sdf_read_stitched_material(sdf_file_t *h);
int sdf_read_stitched_matvar(sdf_file_t *h);
int sdf_read_stitched_species(sdf_file_t *h);
int sdf_read_stitched_obstacle_group(sdf_file_t *h);
int sdf_read_stitched(sdf_file_t *h);
int sdf_read_constant(sdf_file_t *h);

int sdf_read_plain_mesh(sdf_file_t *h);
int sdf_read_plain_mesh_info(sdf_file_t *h);
int sdf_read_plain_variable(sdf_file_t *h);
int sdf_read_plain_variable_info(sdf_file_t *h);

int sdf_read_point_mesh(sdf_file_t *h);
int sdf_read_point_mesh_info(sdf_file_t *h);
int sdf_read_point_variable(sdf_file_t *h);
int sdf_read_point_variable_info(sdf_file_t *h);
int sdf_read_lagran_mesh(sdf_file_t *h);
int sdf_read_station_info(sdf_file_t *h);

void sdf_trim(char *str);
int sdf_header_copy(const sdf_file_t *h_in, sdf_file_t *h_out);


#ifdef SDF_DEBUG
  #define DBG_CHUNK 256

  #define SDF_RANK0 if (h->rank == h->rank_master)

  #define SDF_PRNT(...) SDF_RANK0 do { \
    sprintf(h->dbg, __VA_ARGS__); \
    h->dbg += strlen(h->dbg); \
    if (h->dbg_count - (h->dbg - h->dbg_buf) < DBG_CHUNK) { \
        char *old = h->dbg_buf; \
        h->dbg_buf = malloc(2 * h->dbg_count); \
        memcpy(h->dbg_buf, old, h->dbg_count); \
        h->dbg_count *= 2; \
        h->dbg = h->dbg_buf + (h->dbg - old); \
        free(old); \
    }} while (0)

  #define SDF_DPRNT(...) do { \
        int _a; for (_a=0; _a<h->indent; _a++) SDF_PRNT(" "); \
        SDF_PRNT(__VA_ARGS__); \
    } while (0)

  #define SDF_DPRNTa(a,f,len) SDF_RANK0 { \
            int _b, _c; \
            for (_b=0; _b<len; _b++) { \
                for (_c=0; _c<h->indent; _c++) SDF_PRNT(" "); \
                SDF_PRNT(#a "[%i]: %" #f "\n", _b, a[_b]); \
            } \
        }

 #ifdef SDF_DEBUG_ALL
  #define SDF_DPRNTar(a,len) SDF_RANK0 { \
        int _d, _i; \
        if (b->datatype_out == SDF_DATATYPE_REAL4) { \
            float *arr = (a); \
            SDF_PRNT("r4 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %g", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_REAL8) { \
            double *arr = (a); \
            SDF_PRNT("r8 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %g", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_CHARACTER) { \
            char *arr = (a); \
            SDF_PRNT("c1 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT("%c", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_LOGICAL) { \
            char *arr = (a); \
            SDF_PRNT("l1 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT("%x ", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_INTEGER4) { \
            int *arr = (a); \
            SDF_PRNT("i4 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %i", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_INTEGER8) { \
            uint64_t *arr = (a); \
            SDF_PRNT("i8 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < 10; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %" PRIu64, arr[_d]); \
                } \
            } \
        } \
        SDF_PRNT("\n"); \
    }
 #else
  #define SDF_DPRNTar(a,len) do {} while(0)
 #endif
#else
  #define SDF_DPRNT(...) do {} while(0)
  #define SDF_DPRNTa(a,f,len) do {} while(0)
  #define SDF_DPRNTar(a,len) do {} while(0)
#endif

#define _SDF_CBYTE_SWAP32(val) do { \
        char _c; \
        _c = (val)[0]; (val)[0] = (val)[3]; (val)[3] = _c; \
        _c = (val)[1]; (val)[1] = (val)[2]; (val)[2] = _c; \
    } while(0)

#define _SDF_CBYTE_SWAP64(val) do { \
        char _c; \
        _c = (val)[0]; (val)[0] = (val)[7]; (val)[7] = _c; \
        _c = (val)[1]; (val)[1] = (val)[6]; (val)[6] = _c; \
        _c = (val)[2]; (val)[2] = (val)[5]; (val)[5] = _c; \
        _c = (val)[3]; (val)[3] = (val)[4]; (val)[4] = _c; \
    } while(0)

#define _SDF_BYTE_SWAP64(val) \
    (val) = ((uint64_t)(val)&0xff00000000000000)>>56 | \
            ((uint64_t)(val)&0xff000000000000)>>40 | \
            ((uint64_t)(val)&0xff0000000000)>>24 | \
            ((uint64_t)(val)&0xff00000000)>>8 | \
            ((uint64_t)(val)&0xff000000)<<8 | \
            ((uint64_t)(val)&0xff0000)<<24 | \
            ((uint64_t)(val)&0xff00)<<40 | \
            ((uint64_t)(val)&0xff)<<56

#define _SDF_BYTE_SWAP32(val) \
    (val) = ((uint32_t)(val)&0xff000000)>>24 | \
            ((uint32_t)(val)&0xff0000)>>8 | \
            ((uint32_t)(val)&0xff00)<<8 | \
            ((uint32_t)(val)&0xff)<<24

#define SDF_READ_ENTRY_INT4(value) do { \
        (value) = *((uint32_t *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_BYTE_SWAP32(value); \
        h->current_location += 4; \
        SDF_DPRNT(#value ": %i\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_INT8(value) do { \
        (value) = *((uint64_t *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_BYTE_SWAP64(value); \
        h->current_location += 8; \
        SDF_DPRNT(#value ": %lli\n", (long long int)(value)); \
    } while (0)

#define SDF_READ_ENTRY_REAL4(value) do { \
        (value) = *((float *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_CBYTE_SWAP32(((char*)&value)); \
        h->current_location += 4; \
        SDF_DPRNT(#value ": %g\n", (float)(value)); \
    } while (0)

#define SDF_READ_ENTRY_REAL8(value) do { \
        (value) = *((double *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_CBYTE_SWAP64(((char*)&value)); \
        h->current_location += 8; \
        SDF_DPRNT(#value ": %g\n", (double)(value)); \
    } while (0)

#define SDF_READ_ENTRY_LOGICAL(value) do { \
        (value) = *((char *) \
            (h->buffer + h->current_location - h->start_location)); \
        h->current_location += 1; \
        if ((value)) { \
            SDF_DPRNT(#value ": true\n"); \
        } else { \
            SDF_DPRNT(#value ": false\n"); \
        } \
    } while (0)

#define SDF_READ_ENTRY_CONST(value) do { \
        char _c; \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            SDF_TYPE_SIZES[b->datatype]); \
        h->current_location += SDF_TYPE_SIZES[b->datatype]; \
        switch (b->datatype) { \
        case(SDF_DATATYPE_REAL4): \
          if (h->swap) { \
            _c = (value)[0]; (value)[0] = (value)[3]; (value)[3] = _c; \
            _c = (value)[1]; (value)[1] = (value)[2]; (value)[2] = _c; \
          } \
          SDF_DPRNT(#value ": %g\n", *((float*)(value))); \
          break; \
        case(SDF_DATATYPE_REAL8): \
          if (h->swap) { \
            _c = (value)[0]; (value)[0] = (value)[7]; (value)[7] = _c; \
            _c = (value)[1]; (value)[1] = (value)[6]; (value)[6] = _c; \
            _c = (value)[2]; (value)[2] = (value)[5]; (value)[5] = _c; \
            _c = (value)[3]; (value)[3] = (value)[4]; (value)[4] = _c; \
          } \
          SDF_DPRNT(#value ": %g\n", *((double*)(value))); \
          break; \
        case(SDF_DATATYPE_INTEGER4): \
          if (h->swap) { \
            _c = (value)[0]; (value)[0] = (value)[3]; (value)[3] = _c; \
            _c = (value)[1]; (value)[1] = (value)[2]; (value)[2] = _c; \
          } \
          SDF_DPRNT(#value ": %i\n", *((int32_t*)(value))); \
          break; \
        case(SDF_DATATYPE_INTEGER8): \
          if (h->swap) { \
            _c = (value)[0]; (value)[0] = (value)[7]; (value)[7] = _c; \
            _c = (value)[1]; (value)[1] = (value)[6]; (value)[6] = _c; \
            _c = (value)[2]; (value)[2] = (value)[5]; (value)[5] = _c; \
            _c = (value)[3]; (value)[3] = (value)[4]; (value)[4] = _c; \
          } \
          SDF_DPRNT(#value ": %lli\n", *((long long int*)(value))); \
          break; \
        } \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_INT4(value, length) do { \
        uint32_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(int32_t)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            4 * (length)); \
        if (h->swap) { \
            val = (uint32_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP32(*val); \
                val++; \
            } \
        } \
        h->current_location += 4 * (length); \
        SDF_DPRNTa(value, i, (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_INT8(value, length) do { \
        uint64_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(int64_t)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            8 * (length)); \
        if (h->swap) { \
            val = (uint64_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP64(*val); \
                val++; \
            } \
        } \
        h->current_location += 8 * (length); \
        SDF_DPRNTa(value, lli, (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_REAL4(value, length) do { \
        uint32_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(float)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            4 * (length)); \
        if (h->swap) { \
            val = (uint32_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP32(*val); \
                val++; \
            } \
        } \
        h->current_location += 4 * (length); \
        SDF_DPRNTa(value, g, (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_REAL8(value, length) do { \
        uint64_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(double)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            8 * (length)); \
        if (h->swap) { \
            val = (uint64_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP64(*val); \
                val++; \
            } \
        } \
        h->current_location += 8 * (length); \
        SDF_DPRNTa(value, g, (length)); \
    } while (0)

#define SDF_READ_ENTRY_TYPE(value) do { \
        (b->value) = *((uint32_t *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_BYTE_SWAP32(b->value); \
        h->current_location += 4; \
        if (b->value < sdf_ ## value ## _len) \
            SDF_DPRNT("b->" #value ": %s\n", sdf_ ## value ## _c[b->value]); \
        else \
            SDF_DPRNT("b->" #value ": %i (UNKNOWN)\n", b->value); \
    } while (0)

#define SDF_READ_ENTRY_STRINGLEN(value, length) do { \
        if (!(value)) value = calloc(length+1, sizeof(char)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            (length)); \
        value[length] = '\0'; \
        sdf_trim(value); \
        h->current_location += (length); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_STRINGLEN(value, length, clen) do { \
        int _e; \
        if (!(value)) value = calloc((length+1), sizeof(char*)); \
        for (_e=0; _e<(length); _e++) { \
            SDF_READ_ENTRY_STRINGLEN(value[_e], clen); \
            SDF_DPRNT(#value "[%i]: %s\n", _e, (value[_e])); \
        } \
    } while (0)

#define SDF_READ_ENTRY_ID(value) do { \
        SDF_READ_ENTRY_STRINGLEN(value, SDF_ID_LENGTH); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_ID(value, length) \
        SDF_READ_ENTRY_ARRAY_STRINGLEN(value, length, SDF_ID_LENGTH)

#define SDF_READ_ENTRY_STRING(value) do { \
        SDF_READ_ENTRY_STRINGLEN(value, h->string_length); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_STRING(value, length) \
        SDF_READ_ENTRY_ARRAY_STRINGLEN(value, length, h->string_length)

#define SDF_SET_ENTRY_STRINGLEN(value, strvalue, length) do { \
        if (!(value)) value = malloc(h->string_length+1); \
        strncpy((value), (strvalue), (length)); \
    } while (0)

#define SDF_SET_ENTRY_ID(value, strvalue) do { \
        SDF_SET_ENTRY_STRINGLEN(value, strvalue, SDF_ID_LENGTH); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#define SDF_SET_ENTRY_STRING(value, strvalue) do { \
        SDF_SET_ENTRY_STRINGLEN(value, strvalue, h->string_length); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif
