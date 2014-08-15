/**
   @internal
   @file stack_allocator.h

   @brief Declarations for the SDF stack allocator helper routines.
   @details Routines for stack-based memory management.
   @author Dr Keith Bennett
   @date 15/02/2014
*/

#ifndef _STACK_ALLOCATOR_H_
#define _STACK_ALLOCATOR_H_
#include <sdf.h>

#ifdef __cplusplus
extern "C" {
#endif

void stack_alloc(sdf_block_t *b);
void stack_free_block(sdf_block_t *b);
void stack_push_to_bottom(sdf_block_t *b);
void stack_freeup_memory(void);
void stack_free(void);
void stack_init(void);

#ifdef __cplusplus
}
#endif

#endif
