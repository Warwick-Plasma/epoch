#ifndef _SDF_MODIFY_H_
#define _SDF_MODIFY_H_

#ifdef __cplusplus
extern "C" {
#endif

sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id);
sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name);

#ifdef __cplusplus
}
#endif

#endif
