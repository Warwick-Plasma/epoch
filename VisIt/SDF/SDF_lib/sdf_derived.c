#include <stdlib.h>
#include "sdf.h"



sdf_block_t *sdf_callback_grid_component(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);

    if (!b->grids) b->grids = calloc(1, sizeof(float*));
    b->data = b->grids[0] = mesh->grids[b->nm];
    b->datatype_out = mesh->datatype_out;
    if (b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED)
        b->local_dims[0] = b->nlocal = mesh->nlocal;
    else
        b->local_dims[0] = mesh->local_dims[b->nm];
    return b;
}



sdf_block_t *sdf_callback_cpu_split(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *split_block = sdf_find_block_by_id(h, b->variable_ids[0]);
    int *var, *split, *xsplit, *ysplit, *zsplit;
    int nproc, nprocx, nprocy, nprocz, np, idx, i, j, k;
    int x0, x1, y0, y1, z0, z1;
    int nx, ny, nz, npx, npy, npz;

    nx = b->dims[0];
    ny = b->dims[1];
    nz = b->dims[2];

    if (split_block->geometry == 1) {
        var = split_block->data;
        np = 0;

        nprocx = split_block->dims[0] + 1;
        xsplit = malloc(nprocx * sizeof(int));
        for (i = 0; i < nprocx-1; i++) xsplit[i] = var[np++];
        xsplit[nprocx-1] = nx = b->dims[0];

        if (b->ndims > 1) {
            nprocy = split_block->dims[1] + 1;
            ysplit = malloc(nprocy * sizeof(int));
            for (i = 0; i < nprocy-1; i++) ysplit[i] = var[np++];
            ysplit[nprocy-1] = ny = b->dims[1];
        } else {
            nprocy = 1;
            ysplit = malloc(nprocy * sizeof(int));
            ysplit[nprocy-1] = ny = 1;
        }

        if (b->ndims > 2) {
            nprocz = split_block->dims[2] + 1;
            zsplit = malloc(nprocz * sizeof(int));
            for (i = 0; i < nprocz-1; i++) zsplit[i] = var[np++];
            zsplit[nprocz-1] = nz = b->dims[2];
        } else {
            nprocz = 1;
            zsplit = malloc(nprocz * sizeof(int));
            zsplit[nprocz-1] = nz = 1;
        }

        var = b->data;

        np = 0;
        z0 = 0;
        for (npz = 0; npz < nprocz; npz++) {
            z1 = zsplit[npz];
            y0 = 0;
            for (npy = 0; npy < nprocy; npy++) {
                y1 = ysplit[npy];
                x0 = 0;
                for (npx = 0; npx < nprocx; npx++) {
                    x1 = xsplit[npx];
                    for (k = z0; k < z1; k++) {
                    for (j = y0; j < y1; j++) {
                        idx = (j + k * ny) * nx + x0;
                        for (i = x0; i < x1; i++) {
                            var[idx++] = np;
                    }}}
                    np++;
                    x0 = x1;
                }
                y0 = y1;
            }
            z0 = z1;
        }
        free(xsplit);
        if (b->ndims > 1) free(ysplit);
        if (b->ndims > 2) free(zsplit);

    } else if (split_block->geometry == 2) {
        nprocx = split_block->dims[0];
        nprocy = split_block->dims[1];
        nprocz = split_block->dims[2];

        var = b->data;
        split = split_block->data;
        np = 0;

        if (b->ndims == 2) {
            x0 = 0;
            for (npx = 0; npx < nprocx; npx++) {
                x1 = split[npx];
                y0 = 0;
                for (npy = 0; npy < nprocy; npy++) {
                    y1 = split[nprocx + npy + npx * nprocy];
                    np = npx + npy * nprocx;
                    for (j = y0; j < y1; j++) {
                        idx = j * nx + x0;
                        for (i = x0; i < x1; i++) {
                            var[idx++] = np;
                    }}
                    y0 = y1;
                }
                x0 = x1;
            }
        } else if (b->ndims == 3) {
            x0 = 0;
            for (npx = 0; npx < nprocx; npx++) {
                x1 = split[npx];
                y0 = 0;
                for (npy = 0; npy < nprocy; npy++) {
                    y1 = split[nprocx + npy + npx * nprocy];
                    z0 = 0;
                    for (npz = 0; npz < nprocz; npz++) {
                        z1 = split[nprocx * (1 + nprocy) +
                                   npz + (npy + npx * nprocy) * nprocz];
                        np = npx + (npy + npz * nprocy) * nprocx;
                        for (k = z0; k < z1; k++) {
                        for (j = y0; j < y1; j++) {
                            idx = (j + k * ny) * nx + x0;
                            for (i = x0; i < x1; i++) {
                                var[idx++] = np;
                        }}}
                        z0 = z1;
                    }
                    y0 = y1;
                }
                x0 = x1;
            }
        }

    } else if (split_block->geometry == 3) {
        split = split_block->data;
        nproc = split_block->dims[1];

        var = b->data;

        if (b->ndims == 1) {
            for (np = 0; np < nproc; np++) {
                x0 = split[2*np]-1;
                x1 = split[2*np+1];

                idx = x0;
                for (i = x0; i < x1; i++) {
                    var[idx++] = np;
                }
            }
        } else if (b->ndims == 2) {
            for (np = 0; np < nproc; np++) {
                x0 = split[4*np]-1;
                y0 = split[4*np+1]-1;
                x1 = split[4*np+2];
                y1 = split[4*np+3];

                for (j = y0; j < y1; j++) {
                    idx = j * nx + x0;
                    for (i = x0; i < x1; i++) {
                        var[idx++] = np;
                }}
            }
        } else if (b->ndims == 3) {
            for (np = 0; np < nproc; np++) {
                x0 = split[6*np]-1;
                y0 = split[6*np+1]-1;
                z0 = split[6*np+2]-1;
                x1 = split[6*np+3];
                y1 = split[6*np+4];
                z1 = split[6*np+5];

                for (k = z0; k < z1; k++) {
                for (j = y0; j < y1; j++) {
                    idx = (j + k * ny) * nx + x0;
                    for (i = x0; i < x1; i++) {
                        var[idx++] = np;
                }}}
            }
        }
    }

    return b;
}



int sdf_add_derived_blocks(sdf_file_t *h)
{
    sdf_block_t *b, *cur, *append, *append_head, *append_tail;
    sdf_block_t *mesh;
    int i, len1, len2, nappend = 0;
    char *str, *name1, *name2;
    char *grid_ids[] = { "x", "y", "z" };

    cur = h->current_block;
    append = append_head = calloc(1, sizeof(sdf_block_t));

    b = h->blocklist;
    while (b) {
        if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH
            || b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
            for (i = 0; i < b->ndims; i++) {
                // Add 1d arrays for each coordinate dimension of the
                // particles. (ie. all the x's, all the y's, all the z's).
                // These are required in order to perform scatter plots.
                append->next = calloc(1, sizeof(sdf_block_t));

                nappend++;
                append = append_tail = append->next;
                append->next = NULL;

                if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
                    name1 = b->id;
                    name2 = grid_ids[i];
                    len1 = strlen(name1);
                    len2 = strlen(name2);
                    str = (char*)malloc(len1 + len2 + 2);
                    memcpy(str, name1, len1);
                    str[len1] = '/';
                    memcpy(str+len1+1, name2, len2);
                    str[len1+len2+1] = '\0';
                    append->id = str;

                    append->blocktype = SDF_BLOCKTYPE_POINT_DERIVED;
                    name2 = b->dim_labels[i];
                } else {
                    append->id = NULL;
                    SDF_SET_ENTRY_ID(append->id, grid_ids[i]);

                    append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
                    name2 = grid_ids[i];
                    append->dont_display = 1;
                }

                name1 = b->name;
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 2);
                memcpy(str, name1, len1);
                str[len1] = '/';
                memcpy(str+len1+1, name2, len2);
                str[len1+len2+1] = '\0';
                append->name = str;

                SDF_SET_ENTRY_ID(append->units, b->dim_units[i]);
                SDF_SET_ENTRY_ID(append->mesh_id, b->id);
                append->nm = i;
                append->ndims = 1;
                append->n_ids = 1;
                append->variable_ids = calloc(append->n_ids, sizeof(char*));
                SDF_SET_ENTRY_ID(append->variable_ids[0], b->id);
                append->must_read = calloc(append->n_ids, sizeof(char*));
                append->must_read[0] = 1;
                append->populate_data = sdf_callback_grid_component;
                append->done_header = 1;
                // Hack to prevent storage being allocated for this variable.
                append->dont_allocate = 1;
            }
        } else if (b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT) {
            append->next = calloc(1, sizeof(sdf_block_t));

            nappend++;
            append = append_tail = append->next;
            append->next = NULL;

            len1 = strlen(b->id);
            append->id = b->id;
            b->id = malloc(len1+7);
            memcpy(b->id, append->id, len1);
            memcpy(b->id+len1, "_orig", 6);

            len1 = strlen(b->name);
            append->name = b->name;
            b->name = malloc(len1+7);
            memcpy(b->name, append->name, len1);
            memcpy(b->name+len1, "_orig", 6);

            append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;

            // Find first grid mesh
            mesh = h->blocklist;
            while (mesh) {
                if (mesh->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) break;
                mesh = mesh->next;
            }
            SDF_SET_ENTRY_ID(append->units, "CPU");
            SDF_SET_ENTRY_ID(append->mesh_id, mesh->id);
            append->ndims = mesh->ndims;
            for (i=0; i<3; i++) append->dims[i] = mesh->dims[i] - 1;
            append->n_ids = 1;
            append->variable_ids = calloc(append->n_ids, sizeof(char*));
            SDF_SET_ENTRY_ID(append->variable_ids[0], b->id);
            append->must_read = calloc(append->n_ids, sizeof(char*));
            append->must_read[0] = 1;
            append->populate_data = sdf_callback_cpu_split;
            append->datatype = append->datatype_out = SDF_DATATYPE_INTEGER4;
        }
        b = b->next;
    }

    if (nappend) {
        h->tail->next = append_head->next;
        h->tail = append_tail;
        h->nblocks += nappend;
    }

    h->current_block = cur;

    return 0;
}
