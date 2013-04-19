#ifndef _SDF_VECTOR_TYPE_H_
#define _SDF_VECTOR_TYPE_H_

/* Growable array datatype and helper routines */

typedef struct vector_type vector_t;

struct vector_type {
    int *data;
    int allocated, size;
};

static vector_t *vector_new(void)
{
    vector_t *vector;

    vector = (vector_t*)malloc(sizeof(vector_t));
    vector->allocated = 32;
    vector->size = 0;
    vector->data = (int*)malloc(vector->allocated * sizeof(int));

    return vector;
}

static void vector_push_back(vector_t *vector, int val)
{
    int *data;

    // Grow vector if necessary
    if (vector->size == vector->allocated) {
        vector->allocated = vector->allocated << 1;
        data = (int*)malloc(vector->allocated * sizeof(int));
        memcpy(data, vector->data, vector->size * sizeof(int));
        free(vector->data);
        vector->data = data;
    }

    vector->data[vector->size++] = val;
}

static void vector_truncate(vector_t *vector)
{
    int *data;

    vector->allocated = vector->size;
    data = (int*)malloc(vector->allocated * sizeof(int));
    memcpy(data, vector->data, vector->size * sizeof(int));
    free(vector->data);
    vector->data = data;
}

static void vector_free(vector_t *vector)
{
    free(vector->data);
    free(vector);
}

#endif
