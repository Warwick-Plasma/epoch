#ifndef _SDF_LIST_TYPE_H_
#define _SDF_LIST_TYPE_H_

/* Linked list datatype and helper routines */

typedef struct list_type list_t;


struct list_type {
    struct list_entry {
        void *value;
        struct list_entry *next;
    } *head, *tail;
    int count;
};

static void list_init(list_t **lt_ptr)
{
    list_t *lt = *lt_ptr;

    *lt_ptr = lt = calloc(1, sizeof(*lt));
    lt->head = lt->tail = calloc(1, sizeof(*(lt->head)));
}

static void list_append(list_t *lt, void *value)
{
    lt->tail = lt->tail->next = calloc(1, sizeof(*(lt->head)));
    lt->tail->value = value;
    lt->count++;
}

static void *list_start(list_t *lt)
{
    lt->tail = lt->head->next;
    return (lt->tail?lt->tail->value:NULL);
}

static void *list_next(list_t *lt)
{
    lt->tail = lt->tail->next;
    return (lt->tail?lt->tail->value:NULL);
}

static void list_destroy(list_t **lt_ptr)
{
    list_t *lt = *lt_ptr;
    int i;

    lt->tail = lt->head;
    for (i = 0; i < lt->count+1; i++) {
        lt->head = lt->tail->next;
        free(lt->tail);
        lt->tail = lt->head;
    }
    free(lt);
    *lt_ptr = NULL;
}

#endif
