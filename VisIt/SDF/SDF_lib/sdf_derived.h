#ifndef _SDF_DERIVED_H_
#define _SDF_DERIVED_H_
#include <sdf.h>

#ifdef __cplusplus
extern "C" {
#endif

int sdf_add_derived_blocks(sdf_file_t *h);
int sdf_add_derived_blocks_final(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
