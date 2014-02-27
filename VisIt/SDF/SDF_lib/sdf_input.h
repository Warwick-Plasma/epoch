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

#define SDF_READ_ENTRY_INT4(value) do { \
        (value) = *((int32_t *) \
            (h->buffer + h->current_location - h->start_location)); \
        if (h->swap) _SDF_BYTE_SWAP32(value); \
        h->current_location += 4; \
        SDF_DPRNT(#value ": %i\n", (value)); \
    } while (0)

#define SDF_READ_ENTRY_INT8(value) do { \
        (value) = *((int64_t *) \
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
        int32_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(int32_t)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            4 * (length)); \
        if (h->swap) { \
            val = (int32_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP32(*val); \
                val++; \
            } \
        } \
        h->current_location += 4 * (length); \
        SDF_DPRNTa(value, PRIi32, (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_INT8(value, length) do { \
        int64_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(int64_t)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            8 * (length)); \
        if (h->swap) { \
            val = (int64_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP64(*val); \
                val++; \
            } \
        } \
        h->current_location += 8 * (length); \
        SDF_DPRNTa(value, PRIi64, (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_REAL4(value, length) do { \
        int32_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(float)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            4 * (length)); \
        if (h->swap) { \
            val = (int32_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP32(*val); \
                val++; \
            } \
        } \
        h->current_location += 4 * (length); \
        SDF_DPRNTa(value, "g", (length)); \
    } while (0)

#define SDF_READ_ENTRY_ARRAY_REAL8(value, length) do { \
        int64_t *val; int _i; \
        if (!(value)) value = calloc((length), sizeof(double)); \
        memcpy((value), (h->buffer + h->current_location - h->start_location), \
            8 * (length)); \
        if (h->swap) { \
            val = (int64_t*)(value); \
            for (_i=0; _i<(length); _i++) { \
                _SDF_BYTE_SWAP64(*val); \
                val++; \
            } \
        } \
        h->current_location += 8 * (length); \
        SDF_DPRNTa(value, "g", (length)); \
    } while (0)

#define SDF_READ_ENTRY_TYPE(value) do { \
        (b->value) = *((int32_t *) \
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
