#include <stdlib.h>
#include "SDFStack.h"
#include "sdf.h"

struct stack;

struct stack {
    sdf_block_t *block;
    struct stack *next;
};

static struct stack *stack_head = NULL;
static struct stack *stack_tail = NULL;

#define MAX_MEMORY 2147483648 // 2GB
static uint64_t memory_size = 0;


static inline void stack_alloc(sdf_block_t *b)
{
    struct stack *tail;
    b->data = calloc(b->nlocal, b->type_size_out);
    memory_size += b->nlocal * b->type_size_out;
    stack_tail->next = tail = (struct stack*)malloc(sizeof(struct stack));
    tail->block = b;
    tail->next = NULL;
    stack_tail = tail;
}


static inline void stack_free(void)
{
    sdf_block_t *b;
    struct stack *head;

    if (memory_size < MAX_MEMORY) return;

    while (stack_head->next) {
        head = stack_head;
        stack_head = stack_head->next;
        free(head);
        b = stack_head->block;
        stack_head->block = NULL;
        free(b->data);
        b->data = NULL;
        b->done_data = 0;
        memory_size -= b->nlocal * b->type_size_out;
        if (memory_size < MAX_MEMORY) break;
    }
}


static inline void stack_init(void)
{
    if (!stack_head) stack_head =
        stack_tail = (struct stack*)calloc(1, sizeof(struct stack));
}
