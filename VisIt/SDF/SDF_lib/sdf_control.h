#ifndef _SDF_CONTROL_H_
#define _SDF_CONTROL_H_

#include "sdf_util.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SDF_ID_LENGTH 32
#define SDF_HEADER_LENGTH (11 * SOI4 + 2 * SOI8 + SOF8 + 12 + h->id_length)
#define SDF_BLOCK_HEADER_LENGTH \
    (4 + 3 * SOI4 + 3 * SOI8 + h->id_length + h->string_length)
#define SDF_ENDIANNESS 16911887

int sdf_seek_set(sdf_file_t *h, off_t offset);
int sdf_seek(sdf_file_t *h);
int sdf_free_block(sdf_file_t *h, sdf_block_t *b);
int sdf_broadcast(sdf_file_t *h, void *buf, int size);
int sdf_convert_array_to_float(sdf_file_t *h, void **var_in, int count);
int sdf_factor(sdf_file_t *h);
int sdf_randomize_array(sdf_file_t *h, void **var_in, int count);

#ifdef __cplusplus
}
#endif

#endif
