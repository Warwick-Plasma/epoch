#include <stdlib.h>
#include <sdf.h>
#include "stack_allocator.h"

// ****************************************************************************
//  Memory management stack
// ****************************************************************************

struct stack;

struct stack {
    sdf_block_t *block;
    struct stack *next;
};

static struct stack *stack_head = NULL;
static struct stack *stack_tail = NULL;

#define MAX_MEMORY 2147483648 // 2GB
static int64_t memory_size = 0;


void stack_alloc(sdf_block_t *b)
{
    int i;
    uint64_t sz;
    struct stack *tail;
    if (b->done_data || b->dont_own_data) return;
    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH ||
            b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
        b->ngrids = b->ndims;
        sz = b->ngrids * sizeof(*b->grids);
        b->grids = calloc(1, sz);
        memory_size += sz;
        for (i = 0; i < b->ngrids; i++) {
            sz = b->local_dims[i] * SDF_TYPE_SIZES[b->datatype_out];
            b->grids[i] = calloc(1, sz);
            memory_size += sz;
        }
    } else {
        sz = b->nelements_local * SDF_TYPE_SIZES[b->datatype_out];
        b->data = calloc(1, sz);
        memory_size += sz;
    }
    stack_tail->next = tail = (struct stack*)malloc(sizeof(struct stack));
    tail->block = b;
    tail->next = NULL;
    stack_tail = tail;
}


static void stack_free_data_or_grid(sdf_block_t *b)
{
    int i;

    if (b->grids) {
        for (i = 0; i < b->ngrids; i++) {
            free(b->grids[i]);
            memory_size -= b->local_dims[i] * SDF_TYPE_SIZES[b->datatype_out];
        }
        memory_size -= b->ngrids * sizeof(*b->grids);
    } else {
        free(b->data);
        memory_size -= b->nelements_local * SDF_TYPE_SIZES[b->datatype_out];
    }
    b->grids = NULL;
    b->data = NULL;
    b->done_data = 0;
}


void stack_free_block(sdf_block_t *b)
{
    struct stack *old_stack_entry = stack_head;
    struct stack *stack_entry = stack_head;

    while (stack_entry) {
        if (stack_entry->block == b) {
            stack_free_data_or_grid(b);
            old_stack_entry->next = stack_entry->next;
            if (stack_entry == stack_tail) stack_tail = old_stack_entry;
            free(stack_entry);
            return;
        }
        old_stack_entry = stack_entry;
        stack_entry = stack_entry->next;
    }
}


void stack_push_to_bottom(sdf_block_t *b)
{
    struct stack *old_stack_entry = stack_head;
    struct stack *stack_entry = stack_head;

    while (stack_entry) {
        if (stack_entry->block == b) {
            old_stack_entry->next = stack_entry->next;
            stack_tail->next = stack_entry;
            stack_tail = stack_entry;
            stack_tail->next = NULL;
            return;
        }
        old_stack_entry = stack_entry;
        stack_entry = stack_entry->next;
    }
}


void stack_freeup_memory(void)
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
        stack_free_data_or_grid(b);
        if (memory_size < MAX_MEMORY) break;
    }
}


void stack_free(void)
{
    sdf_block_t *b;
    struct stack *head;

    while (stack_head->next) {
        head = stack_head;
        stack_head = stack_head->next;
        free(head);
        b = stack_head->block;
        stack_head->block = NULL;
        stack_free_data_or_grid(b);
    }
    memory_size = 0;
}


void stack_init(void)
{
    if (!stack_head) stack_head =
        stack_tail = (struct stack*)calloc(1, sizeof(struct stack));
}
