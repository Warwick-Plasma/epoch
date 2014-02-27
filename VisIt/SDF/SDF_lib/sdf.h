/**
   @file sdf.h

   @brief Declarations for the SDF C-library.
   @details Routines for reading and writing SDF files.
   @author Dr Keith Bennett
   @date 15/02/2014
   @version 5.0
*/

#ifndef _SDF_COMMON_H_
#define _SDF_COMMON_H_

#include <stdio.h>
#include <inttypes.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#define SDF_VERSION  1
#define SDF_REVISION 2
#define SDF_LIB_VERSION  6
#define SDF_LIB_REVISION 0

#define SDF_MAGIC "SDF1"

#define SDF_MAXDIMS 4

#ifdef __cplusplus
extern "C" {
#endif

/**
  Constants used for identifying the block's type
 */
enum sdf_blocktype {
    /** Deleted block. Should be ignored. */
    SDF_BLOCKTYPE_SCRUBBED = -1,
    /** Unknown block type. This is an error. */
    SDF_BLOCKTYPE_NULL = 0,
    /** Block describing a plain mesh or grid. */
    SDF_BLOCKTYPE_PLAIN_MESH,
    /** Block describing a point mesh or grid. */
    SDF_BLOCKTYPE_POINT_MESH,
    /** Block describing a variable on a plain mesh. */
    SDF_BLOCKTYPE_PLAIN_VARIABLE,
    /** Block describing a variable on a point mesh. */
    SDF_BLOCKTYPE_POINT_VARIABLE,
    /** A simple constant not associated with a grid. */
    SDF_BLOCKTYPE_CONSTANT,
    /** A simple array not associated with a grid. */
    SDF_BLOCKTYPE_ARRAY,
    /** Information about the simulation. */
    SDF_BLOCKTYPE_RUN_INFO,
    /** Embedded source code block. */
    SDF_BLOCKTYPE_SOURCE,
    /** List of blocks to combine as a tensor or vector. */
    SDF_BLOCKTYPE_STITCHED_TENSOR,
    /** List of blocks to combine as a multi-material mesh. */
    SDF_BLOCKTYPE_STITCHED_MATERIAL,
    /** List of blocks to combine as a multi-material variable. */
    SDF_BLOCKTYPE_STITCHED_MATVAR,
    /** List of blocks to combine as a species mesh. This is similar to a
      * multi-material mesh except there is no interface in a mixed cell. */
    SDF_BLOCKTYPE_STITCHED_SPECIES,
    /** Information about a particle species. */
    SDF_BLOCKTYPE_SPECIES,
    /** This blocktype is never actually written to an SDF file. It is used
      * within the C-library and VisIt to represent a plain variable whose
      * content is generated dynamically based on other data in the file. */
    SDF_BLOCKTYPE_PLAIN_DERIVED,
    /** As above, this blocktype is never actually written to an SDF file. It
      * serves the same purpose as SDF_BLOCKTYPE_PLAIN_DERIVED, except the
      * variable is defined on a point mesh. */
    SDF_BLOCKTYPE_POINT_DERIVED,
    /** This is the same as SDF_BLOCKTYPE_STITCHED_TENSOR, except that all the
      * data for the stitched variables is contained in the data section of
      * this block rather than the blocks which are referenced. */
    SDF_BLOCKTYPE_CONTIGUOUS_TENSOR,
    /** Same as above, for SDF_BLOCKTYPE_STITCHED_MATERIAL */
    SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL,
    /** Same as above, for SDF_BLOCKTYPE_STITCHED_MATVAR */
    SDF_BLOCKTYPE_CONTIGUOUS_MATVAR,
    /** Same as above, for SDF_BLOCKTYPE_STITCHED_SPECIES */
    SDF_BLOCKTYPE_CONTIGUOUS_SPECIES,
    /** Information about the parallel domain decomposition. */
    SDF_BLOCKTYPE_CPU_SPLIT,
    /** List of blocks to combine as an obstacle mesh. */
    SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP,
    /** Block describing an unstructured mesh or grid. */
    SDF_BLOCKTYPE_UNSTRUCTURED_MESH,
    /** Block describing a stitched variable.
      * This allows any arbitrary set of variables to be grouped together. */
    SDF_BLOCKTYPE_STITCHED,
    /** This is the same as SDF_BLOCKTYPE_STITCHED, except that all the
      * data for the stitched variables is contained in the data section of
      * this block rather than the blocks which are referenced. */
    SDF_BLOCKTYPE_CONTIGUOUS,
    /** Block describing a Lagrangian mesh or grid. */
    SDF_BLOCKTYPE_LAGRANGIAN_MESH,
    /** Block describing a station point. */
    SDF_BLOCKTYPE_STATION,
    /** As with SDF_BLOCKTYPE_PLAIN_DERIVED, this blocktype is never actually
      * written to an SDF file. It serves the same purpose as
      * SDF_BLOCKTYPE_PLAIN_DERIVED, except the variable is defined as a
      * station variable. */
    SDF_BLOCKTYPE_STATION_DERIVED,
};


/**
The "geometry_type" specifies the geometry of the current block and it can
take one of the following values:
*/
enum sdf_geometry {
    SDF_GEOMETRY_NULL = 0,    /**< Unspecified geometry. This is an error. */
    SDF_GEOMETRY_CARTESIAN,   /**< Cartesian geometry. */
    SDF_GEOMETRY_CYLINDRICAL, /**< Cylindrical geometry. */
    SDF_GEOMETRY_SPHERICAL,   /**< Spherical geometry. */
};


/**
The mesh associated with a variable is always node-centred, ie. the values
written as mesh data specify the nodal values of a grid. Variables may be
defined at points which are offset from this grid due to grid staggering in
the code. The "stagger" entry specifies where the variable is defined
relative to the mesh. Since we have already defined the number of points
that the associated mesh contains, this determines how many points are required
to display the variable. The entry is represented by a bit-mask where each
bit corresponds to a shift in coordinates by half a grid cell in the direction
given by the bit position. Therefore the value "1" (or "0001" in binary)
is a shift by \e dx/2 in the \e x direction, "2" (or "0010" in binary) is
a shift by \e dy/2 in the \e y direction and "4" (or "0100" in binary) is
a shift by \e dz/2 in the \e z direction. These can be combined to give shifts
in more than one direction. The system can also be extended to account for
more than three directions (eg. "8" for direction 4).

For convenience, a list of pre-defined constants are defined for the typical
cases.
 */
enum sdf_stagger {
    /** Cell centred. At the midpoint between nodes.
      * Implies an <em>(nx,ny,nz)</em> grid. */
    SDF_STAGGER_CELL_CENTRE = 0,
    /** Face centred in X. Located at the midpoint between nodes on the Y-Z
      * plane.
      * Implies an <em>(nx+1,ny,nz)</em> grid. */
    SDF_STAGGER_FACE_X,
    /** Face centred in Y. Located at the midpoint between nodes on the X-Z
      * plane.
      * Implies an <em>(nx,ny+1,nz)</em> grid. */
    SDF_STAGGER_FACE_Y,
    /** Face centred in Z. Located at the midpoint between nodes on the X-Y
      * plane.
      * Implies an <em>(nx,ny,nz+1)</em> grid. */
    SDF_STAGGER_EDGE_Z,
    /** Edge centred along X. Located at the midpoint between nodes along the
      * X-axis.
      * Implies an <em>(nx,ny+1,nz+1)</em> grid. */
    SDF_STAGGER_FACE_Z,
    /** Edge centred along Y. Located at the midpoint between nodes along the
      * Y-axis.
      * Implies an <em>(nx+1,ny,nz+1)</em> grid. */
    SDF_STAGGER_EDGE_Y,
    /** Edge centred along Z. Located at the midpoint between nodes along the
      * Z-axis.
      * Implies an <em>(nx+1,ny+1,nz)</em> grid. */
    SDF_STAGGER_EDGE_X,
    /** Node centred. At the same place as the mesh.
      * Implies an <em>(nx+1,ny+1,nz+1)</em> grid. */
    SDF_STAGGER_VERTEX,
};


/**
 * The datatype specifies the numberical representation of the block's data
 * array.
 */
enum sdf_datatype {
    SDF_DATATYPE_NULL = 0, /**< No datatype specified. This is an error. */
    SDF_DATATYPE_INTEGER4, /**< 4-byte integers. */
    SDF_DATATYPE_INTEGER8, /**< 8-byte integers. */
    SDF_DATATYPE_REAL4,    /**< 4-byte floating point (ie. single precision). */
    SDF_DATATYPE_REAL8,    /**< 8-byte floating point (ie. double precision). */
    SDF_DATATYPE_REAL16,   /**< 16-byte floating point (ie. quad precision). */
    SDF_DATATYPE_CHARACTER,/**< 1-byte characters. */
    SDF_DATATYPE_LOGICAL,  /**< Logical variables. (Represented as 1-byte
                               characters */
    SDF_DATATYPE_OTHER,    /**< Unspecified datatype. The type of data in the
                                block must be inferred from the block type. */
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
    int64_t dims[3], local_dims[3];
    int64_t block_start, next_block_location, data_location;
    int64_t inline_block_start, inline_next_block_location;
    int64_t summary_block_start, summary_next_block_location;
    int64_t nelements, nelements_local, data_length;
    int64_t *nelements_blocks, *data_length_blocks;
    int64_t *array_starts, *array_ends, *array_strides;
    int64_t *global_array_starts, *global_array_ends, *global_array_strides;
    int32_t ndims, geometry, datatype, blocktype, info_length;
    int32_t type_size, stagger, datatype_out, type_size_out;
    int32_t nstations, nvariables, step, step_increment;
    int32_t *dims_in, *station_nvars, *variable_types, *station_index;
    int32_t *station_move;
    int nm, n_ids, opt, ng, nfaces, ngrids;
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
    int cpu_split[SDF_MAXDIMS], starts[SDF_MAXDIMS];
    int proc_min[3], proc_max[3];
#ifdef PARALLEL
    MPI_Datatype mpitype, distribution, mpitype_out;
#endif
};

struct sdf_file {
    int64_t dbg_count;
    int32_t sdf_lib_version, sdf_lib_revision;
    int32_t sdf_extension_version, sdf_extension_revision;
    int32_t file_version, file_revision;
    char *dbg, *dbg_buf, **extension_names;
    // Lines above should never be changed.
    // Lines below must be changed with care and the SDF_LIB_VERSION bumped
    // if the resulting struct is not aligned the same.
    double time;
    int64_t first_block_location, summary_location, start_location, soi, sof;
    int64_t current_location;
    int32_t jobid1, jobid2, endianness, summary_size;
    int32_t block_header_length, string_length, id_length;
    int32_t code_io_version, step;
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
    int array_count;
#ifdef PARALLEL
    MPI_File filehandle;
#else
    FILE *filehandle;
#endif
    comm_t comm;
};

struct run_info {
    int64_t defines;
    int32_t version, revision, compile_date, run_date, io_date, minor_rev;
    char *commit_id, *sha1sum, *compile_machine, *compile_flags;
};


/**
 @brief Open an SDF file and return a file handle.

 This routine attempts to open the given filename as an SDF file.

 @param[in] filename Name of the SDF file to open
 @param[in] comm     MPI communicator to use
 @param[in] mode     File mode
 @param[in] use_mmap Flag which specifies wether mmap should be used

 @return SDF filehandle on success, NULL on error

 Example usage:
 @code
    sdf_file_t *h = sdf_open("myfile.sdf", MPI_COMM_WORLD, SDF_READ, 0);
 @endcode
 */
sdf_file_t *sdf_open(const char *filename, comm_t comm, int mode, int use_mmap);


/**
 @brief Close an SDF file and free the filehandle

 This routine closes the SDF file associated with this file handle
 and frees the file handle along with all associated data.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_close(sdf_file_t *h);


/**
 @brief Free all data blocks in the blocklist.

 This routine cycles through the blocklist and frees any data blocks that
 have been allocated, whilst leaving the metadata intact.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_free_blocklist_data(sdf_file_t *h);


/**
 @brief Find a block with the given ID

 This function finds a block in file's blocklist whose ID field matches the
 one given.

 @param[in] h        SDF file handle
 @param[in] id       ID of the block to find

 @return Pointer to SDF block on success, NULL on error
 */
sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id);


/**
 @brief Find a block with the given name

 This function finds a block in file's blocklist whose name matches the
 one given.

 @param[in] h        SDF file handle
 @param[in] name     Name of the block to find

 @return Pointer to SDF block on success, NULL on error
 */
sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name);


/**
 @ingroup input
 @{
 @brief Read the header of an SDF file

 The sdf_open reads just enough information to check that the file is a
 valid SDF file. This routine populates information about the file such as
 the number of blocks, etc.
 This information is required in order to read the rest of the file's contents.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_header(sdf_file_t *h);


/**
 @brief Read the summary section of an SDF file

 This routine reads the entire summary section of an SDF file into a buffer.
 No parsing is done. If the file does not contain a summary section then the
 buffer is populated by reading the inline metadata blocks.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_summary(sdf_file_t *h);


/**
 @brief Reads the metadata contents of the SDF file and populates a blocklist

 This routine reads the summary from an SDF file and then parses each of
 the metadata blocks. It builds a linked list of SDF block structures which
 contains the metadata information.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_blocklist(sdf_file_t *h);


/**
 @brief Reads the metadata contents of the SDF file and populates a blocklist
        containing both file contents and any derived blocks.

 This routine reads the summary from an SDF file and then parses each of
 the metadata blocks. It builds a linked list of SDF block structures which
 contains the metadata information.
 In addition, any derived variables which can be calculated based on the files
 contents are added to the blocklist.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_blocklist_all(sdf_file_t *h);


/**
 @brief Reads the metadata for a single block in the SDF file

 This routine reads and parses the metadata for the current block pointer.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_block_info(sdf_file_t *h);


/**
 @brief Reads the data for a single block in the SDF file

 This routine reads the data block for the current block pointer.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_read_data(sdf_file_t *h);
/**@}*/


/**
 @brief Returns the bounds for the current parallel domain decomposition.

 @param[in]  h          SDF file handle
 @param[in]  rank       The processor rank to return the extents for
 @param[out] starts     The starting index within the global domain for each
                        dimension.
 @param[out] local_dims The length of the local domain in each dimension.

 @return 0 on success, 1 on error
 */
int sdf_get_domain_bounds(sdf_file_t *h, int rank,
                          int *starts, int *local_dims);


/**
 @brief Modify an array in-place.

 This function modifies an entire array in-place.

 @param[in] h        SDF file handle
 @param[in] b        SDF block on which to act
 @param[in] data     Supplied data array

 @return 0 on success, 1 on error
 */
int sdf_modify_array(sdf_file_t *h, sdf_block_t *b, void *data);


/**
 @brief Modify an array slice in-place.

 This function allows an array to be modified in-place.

 @param[in] h        SDF file handle
 @param[in] b        SDF block on which to act
 @param[in] data     Supplied data array
 @param[in] starts   Array of starts
 @param[in] ends     Array of ends

 @return 0 on success, 1 on error

 Example usage:
 @code
    memcpy(starts, b->dims, sizeof(*starts));
    memcpy(ends, b->dims, sizeof(*ends));
    starts[2] = 3;
    ends[2] = 4;
    sdf_modify_array_section(h, b, data, starts, ends);
 @endcode
 */
int sdf_modify_array_section(sdf_file_t *h, sdf_block_t *b, void *data,
                             int64_t *starts, int64_t *ends);


/**
 @brief Modify an array element in-place.

 This function modifies a single element of an array in-place.

 @param[in] h        SDF file handle
 @param[in] b        SDF block on which to act
 @param[in] data     Supplied data array
 @param[in] index    Array containing the element index

 @return 0 on success, 1 on error

 Example usage:
 @code
    double value = 2;
    int64_t index[] = {1,2,3};
    sdf_modify_array_element(h, b, &value, index);
 @endcode
 */
int sdf_modify_array_element(sdf_file_t *h, sdf_block_t *b, void *data,
                             int64_t *index);


/**
 @brief Add a block in-place.

 This function adds the specified block to the SDF file.
 The block is appended to the end of both the file and blocklist and the
 SDF file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] block    The block to add to the file

 @return 0 on success, 1 on error
 */
int sdf_modify_add_block(sdf_file_t *h, sdf_block_t *block);


/**
 @brief Add the copy of a block in-place.

 This function copies the specified block and adds it to the SDF file.
 The block copy is appended to the end of both the file and blocklist and the
 SDF file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] block    The block to copy and add to the file

 @return 0 on success, 1 on error
 */
int sdf_modify_add_block_copy(sdf_file_t *h, sdf_block_t *block);


/**
 @brief Remove a block in-place.

 This function removes a block from an SDF file.
 The block's associated pointers are removed and deallocated and the SDF
 file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] block    The block to remove from the file

 @return 0 on success, 1 on error
 */
int sdf_modify_remove_block(sdf_file_t *h, sdf_block_t *block);


/**
 @brief Remove a block in-place by ID.

 This function removes a block from an SDF file.
 The block's associated pointers are removed and deallocated and the SDF
 file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] id       The ID string for the block to remove

 @return 0 on success, 1 on error
 */
int sdf_modify_remove_block_id(sdf_file_t *h, const char *id);


/**
 @brief Remove a block in-place by name.

 This function removes a block from an SDF file.
 The block's associated pointers are removed and deallocated and the SDF
 file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] name     The name of the block to remove

 @return 0 on success, 1 on error
 */
int sdf_modify_remove_block_name(sdf_file_t *h, const char *name);


/**
 @brief Add a single material in-place.

 This function adds a single material to a stitched block.
 The block pointer for the material is added to the blocklist
 and the SDF file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] stitched The stitched material block to remove the material from
 @param[in] material The material to remove

 @return 0 on success, 1 on error
 */
int sdf_modify_add_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material);


/**
 @brief Remove a single material in-place.

 This function removes a single material from a stitched block.
 The material's associated block pointers are removed and deallocated
 and the SDF file headers are updated in-place.

 @param[in] h        SDF file handle
 @param[in] stitched The stitched material block to remove the material from
 @param[in] material The material to remove

 @return 0 on success, 1 on error
 */
int sdf_modify_remove_material(sdf_file_t *h, sdf_block_t *stitched,
        sdf_block_t *material);


/**
 @brief Rewrites the metadata for the file.

 The sdf_modify_* routines merely update the blocklist for the file and do
 not actually perform any changes on the file itself. This routine runs through
 the blocklist and updates the metadata in the file where necessary.

 @param[in] h        SDF file handle

 @return 0 on success, 1 on error
 */
int sdf_modify_rewrite_metadata(sdf_file_t *h);


/**
 @brief Set the array section to use for a block.

 This routine sets up parameters for a block so that only a subsection of
 the total array is read into the data block.

 @param[in] b        SDF block on which to act
 @param[in] ndims    Size of starts, ends and strides arrays
 @param[in] starts   Array of start indices
 @param[in] ends     Array of end indices
 @param[in] strides  Array of strides

 @return 0 on success, 1 on error
 */
int sdf_block_set_array_section(sdf_block_t *b, const int ndims,
                                const int64_t *starts, const int64_t *ends,
                                const int64_t *strides);

#ifdef __cplusplus
}
#endif

#endif
