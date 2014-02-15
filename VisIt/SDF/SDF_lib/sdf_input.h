/**
   @file sdf_input.h

   @brief Declarations for the SDF C-library.
   @details Routines for reading and writing SDF files.
   @author Dr Keith Bennett
   @date 15/02/2014
   @version 5.0
*/

#ifndef _SDF_INPUT_H_
#define _SDF_INPUT_H_

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include "sdf_util.h"

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
        SDF_READ_ENTRY_STRINGLEN(value, h->id_length); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_ID(value, length) \
        SDF_READ_ENTRY_ARRAY_STRINGLEN(value, length, h->id_length)

#define SDF_READ_ENTRY_STRING(value) do { \
        SDF_READ_ENTRY_STRINGLEN(value, h->string_length); \
        SDF_DPRNT(#value ": %s\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_STRING(value, length) \
        SDF_READ_ENTRY_ARRAY_STRINGLEN(value, length, h->string_length)

#ifdef __cplusplus
extern "C" {
#endif

void sdf_trim(char *str);
int sdf_read_bytes(sdf_file_t *h, char *buf, int buflen);

#ifdef __cplusplus
}
#endif

#endif
