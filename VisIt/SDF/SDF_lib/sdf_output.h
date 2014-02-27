#ifndef _SDF_OUTPUT_H_
#define _SDF_OUTPUT_H_

#ifdef __cplusplus
extern "C" {
#endif

int sdf_write_bytes(sdf_file_t *h, void *buf, int buflen);
int sdf_write_at(sdf_file_t *h, off_t offset, void *buf, int buflen);
int sdf_flush(sdf_file_t *h);
int64_t sdf_write_new_summary(sdf_file_t *h);
int sdf_write_meta(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
