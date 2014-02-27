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
    struct stack *tail;
    if (b->done_data || b->dont_own_data) return;
    b->data = calloc(b->nelements_local, SDF_TYPE_SIZES[b->datatype_out]);
    memory_size += b->nelements_local * SDF_TYPE_SIZES[b->datatype_out];
    stack_tail->next = tail = (struct stack*)malloc(sizeof(struct stack));
    tail->block = b;
    tail->next = NULL;
    stack_tail = tail;
}


void stack_free_block(sdf_block_t *b)
{
    struct stack *old_stack_entry = stack_head;
    struct stack *stack_entry = stack_head;

    while (stack_entry) {
        if (stack_entry->block == b) {
            free(b->data);
            b->data = NULL;
            b->done_data = 0;
            memory_size -= b->nelements_local * SDF_TYPE_SIZES[b->datatype_out];
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
        free(b->data);
        b->data = NULL;
        b->done_data = 0;
        memory_size -= b->nelements_local * SDF_TYPE_SIZES[b->datatype_out];
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
        free(b->data);
        b->data = NULL;
        b->done_data = 0;
    }
    memory_size = 0;
}


void stack_init(void)
{
    if (!stack_head) stack_head =
        stack_tail = (struct stack*)calloc(1, sizeof(struct stack));
}
