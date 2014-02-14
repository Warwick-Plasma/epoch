#ifndef _SDF_COMMON_H_
#define _SDF_COMMON_H_

#include <stdio.h>
#include <inttypes.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#define SDF_VERSION  1
#define SDF_REVISION 2
#define SDF_LIB_VERSION  5
#define SDF_LIB_REVISION 0

#define SDF_MAGIC "SDF1"

#define SDF_MAXDIMS 4

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
    uint64_t inline_block_start, inline_next_block_location;
    uint64_t summary_block_start, summary_next_block_location;
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
    char dont_own_data, use_mult, next_block_modified, rewrite_metadata;
    char in_file;
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
    uint32_t block_header_length, string_length, id_length;
    uint32_t code_io_version, step;
    int32_t nblocks, nblocks_file, error_code;
    int rank, ncpus, ndomains, rank_master, indent, print;
    char *buffer, *filename;
    char done_header, restart_flag, other_domains, use_float, use_summary;
    char use_random, station_file, swap;
    char inline_metadata_read, summary_metadata_read;
    char inline_metadata_invalid, summary_metadata_invalid, tmp_flag;
    char metadata_modified, can_truncate, first_block_modified;
    char *code_name, *error_message;
    sdf_block_t *blocklist, *tail, *current_block, *last_block_in_file;
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
int sdf_get_domain_bounds(sdf_file_t *h, int rank,
                          int *starts, int *local_dims);

int sdf_modify_array(sdf_file_t *h, sdf_block_t *b, void *data);
int sdf_modify_array_section(sdf_file_t *h, sdf_block_t *b, void *data,
                             uint64_t *starts, uint64_t *ends);
int sdf_modify_array_element(sdf_file_t *h, sdf_block_t *b, void *data,
                             uint64_t *index);
int sdf_modify_add_block(sdf_file_t *h, sdf_block_t *block);
int sdf_modify_add_block_copy(sdf_file_t *h, sdf_block_t *block);
int sdf_modify_remove_block(sdf_file_t *h, sdf_block_t *block);
int sdf_modify_remove_block_id(sdf_file_t *h, const char *id);
int sdf_modify_remove_block_name(sdf_file_t *h, const char *name);
int sdf_modify_add_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material);
int sdf_modify_remove_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material);
int sdf_modify_rewrite_metadata(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
