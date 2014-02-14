#ifndef _SDF_UTIL_H_
#define _SDF_UTIL_H_

#define SOI4  4
#define SOI8  8
#define SOF4  4
#define SOF8  8

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

#ifdef __cplusplus
extern "C" {
#endif

sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id);
sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name);

#ifdef __cplusplus
}
#endif

#endif
