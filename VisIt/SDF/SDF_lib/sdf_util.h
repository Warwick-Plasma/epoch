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

#ifdef SDF_DEBUG_ALL
#define SDF_DEBUG
#endif

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
                SDF_PRNT(#a "[%i]: %" f "\n", _b, a[_b]); \
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
                for (_i=0; _i < h->array_count; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %g", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_REAL8) { \
            double *arr = (a); \
            SDF_PRNT("r8 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < h->array_count; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %g", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_CHARACTER) { \
            char *arr = (a); \
            SDF_PRNT("c1 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < h->array_count; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT("%c", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_LOGICAL) { \
            char *arr = (a); \
            SDF_PRNT("l1 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < h->array_count; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT("%x ", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_INTEGER4) { \
            int *arr = (a); \
            SDF_PRNT("i4 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < h->array_count; _i++, _d++) { \
                    if (_d == (len)) break; \
                    SDF_PRNT(" %i", arr[_d]); \
                } \
            } \
        } else if (b->datatype_out == SDF_DATATYPE_INTEGER8) { \
            int64_t *arr = (a); \
            SDF_PRNT("i8 "); \
            _d=0; while (_d<(len)) { \
                SDF_PRNT("\n%i ",_d); \
                for (_i=0; _i < h->array_count; _i++, _d++) { \
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

#ifdef __cplusplus
extern "C" {
#endif

sdf_block_t *sdf_find_block_by_id(sdf_file_t *h, const char *id);
sdf_block_t *sdf_find_block_by_name(sdf_file_t *h, const char *name);

#ifdef __cplusplus
}
#endif

#endif
