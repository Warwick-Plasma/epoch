/**
   @file sdf_helper.h

   @brief Helper routines for the SDF C-library.
   @details High-level routines for simplifying the use of the SDF library.
   @author Dr Keith Bennett
   @date 19/02/2014
*/

#ifndef _SDF_HELPER_H_
#define _SDF_HELPER_H_
#include <sdf.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 @brief Populate data for and SDF block.

 This routine does all the work necessary in order to populate
 the data in an SDF block. This includes allocating memory, if
 necessary, reading derived block dependencies and populating
 derived data using callback functions.

 @param[in] h        SDF file handle
 @param[in] b        The block whose data will be populated

 @return 0 on success, 1 on error
 */
int sdf_helper_read_data(sdf_file_t *h, sdf_block_t *b);

#ifdef __cplusplus
}
#endif

#endif
