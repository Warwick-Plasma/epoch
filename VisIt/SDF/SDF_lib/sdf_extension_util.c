#include <dlfcn.h>
#include "sdf.h"
#include "sdf_derived.h"
#include "sdf_extension.h"


static void *sdf_global_extension = NULL;
static void *sdf_global_extension_dlhandle = NULL;


void *sdf_extension_load(sdf_file_t *h)
{
    sdf_extension_create_t *sdf_extension_create;

    if (sdf_global_extension) return sdf_global_extension;

    sdf_global_extension_dlhandle = dlopen("sdf_extension.so", RTLD_LAZY);

    if (!sdf_global_extension_dlhandle) {
        h->error_message = dlerror();
        return NULL;
    }

    sdf_extension_create =
            (sdf_extension_create_t *)dlsym(sdf_global_extension_dlhandle,
            "sdf_extension_create");

    sdf_global_extension = sdf_extension_create(h);

    return sdf_global_extension;
}


void sdf_extension_unload(void)
{
    sdf_extension_destroy_t *sdf_extension_destroy;

    if (!sdf_global_extension_dlhandle) return;

    if (sdf_global_extension) {
        sdf_extension_destroy =
                (sdf_extension_destroy_t *)dlsym(sdf_global_extension_dlhandle,
                "sdf_extension_destroy");

        sdf_extension_destroy(sdf_global_extension);
    }

    dlclose(sdf_global_extension_dlhandle);

    sdf_extension_destroy = NULL;
    sdf_global_extension = NULL;
    sdf_global_extension_dlhandle = NULL;

    return;
}


int sdf_read_blocklist_all(sdf_file_t *h)
{
    sdf_extension_t *ext;
    char **preload;
    sdf_block_t *b, *cur;

    // Retrieve the extended interface library from the plugin manager
    sdf_extension_load(h);

    sdf_read_blocklist(h);

    // Append derived data to the blocklist using built-in library.
    sdf_add_derived_blocks(h);

    ext = sdf_global_extension;
    if (ext) {
        preload = ext->preload(ext, h);
        // For each entry in the preload array, try to find the block
        // and populate its data.
        if (preload) {
            int n = 0;
            cur = h->current_block;
            while(preload[n]) {
                b = sdf_find_block_by_id(h, preload[n]);
                if (b && !b->data) {
                    h->current_block = b;
                    if (!b->done_data && !b->dont_own_data) {
                        if (b->data) free(b->data);
                        b->data = calloc(b->nelements_local,
                                SDF_TYPE_SIZES[b->datatype_out]);
                    }
                    sdf_read_data(h);
                }
                free(preload[n]);
                n++;
            }
            free(preload);
            h->current_block = cur;
        }

        // Append derived data to the blocklist using the extension library.
        ext->read_blocklist(ext, h);
    }

    // Append additional derived data for blocks added by the extension.
    sdf_add_derived_blocks_final(h);

    return 0;
}
