#include <stdlib.h>
#include "sdf.h"

#define IJK(i,j,k) ((i) + nx * ((j) + ny * (k)))
#define IJK1(i,j,k) ((i) + (nx+1) * ((j) + (ny+1) * (k)))
#define IJK2(i,j,k) ((i)+1 + (nx+2) * ((j)+1 + (ny+2) * ((k)+1)))


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



sdf_block_t *sdf_callback_boundary_mesh(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *grid = sdf_find_block_by_id(h, b->subblock->mesh_id);
    sdf_block_t *current_block = h->current_block;

    int i, j, k, ii, jj, kk, i0, j0, k0, i1, j1, k1;
    int imin, imax, jmin, jmax, kmin, kmax;
    int nx, ny, nz;
    int gotmat, gotobst, npoints;
    int *vertex_index, *obdata, *ijk;
    float *vertex;

    vector_t *faces = vector_new();
    vector_t *boundary = vector_new();
    vector_t *vertijk = vector_new();

    h->current_block = grid;
    sdf_read_data(h);
    h->current_block = b->subblock;
    sdf_read_data(h);
    h->current_block = current_block;

    obdata = (int *)b->subblock->data;
    nx = b->subblock->local_dims[0] - 2;
    ny = b->subblock->local_dims[1] - 2;
    nz = b->subblock->local_dims[2] - 2;

#ifdef PARALLEL
    imin = (b->subblock->proc_min[0] < 0) ? 0  : nx;
    imax = (b->subblock->proc_max[0] < 0) ? nx : 0;
    jmin = (b->subblock->proc_min[1] < 0) ? 0  : ny;
    jmax = (b->subblock->proc_max[1] < 0) ? ny : 0;
    kmin = (b->subblock->proc_min[2] < 0) ? 0  : nz;
    kmax = (b->subblock->proc_max[2] < 0) ? nz : 0;
#else
    imin = 0; imax = nx;
    jmin = 0; jmax = ny;
    kmin = 0; kmax = nz;
#endif

    switch(b->nm) {
        case 1: // x_min
            imax = 0;
            jmin = ny; jmax = 0;
            kmin = nz; kmax = 0;
            break;
        case 2: // x_max
            imin = nx;
            jmin = ny; jmax = 0;
            kmin = nz; kmax = 0;
            break;
        case 3: // y_min
            imin = nx; imax = 0;
            jmax = 0;
            kmin = nz; kmax = 0;
            break;
        case 4: // y_max
            imin = nx; imax = 0;
            jmin = ny;
            kmin = nz; kmax = 0;
            break;
        case 5: // z_min
            imin = nx; imax = 0;
            jmin = ny; jmax = 0;
            kmax = 0;
            break;
        case 6: // z_max
            imin = nx; imax = 0;
            jmin = ny; jmax = 0;
            kmin = nz;
            break;
    }

    npoints = 0;
    vertex_index = (int*)malloc((nx+1)*(ny+1)*(nz+1)*sizeof(int));

    ii = 0;
    for (i = imin; i <= imax; i += nx) {
        k0 = 0;
        for (k = 0; k <= nz; k++) {
            k1 = (k == nz) ? k : nz-1;
            j0 = 0;
            for (j = 0; j <= ny; j++) {
                j1 = (j == ny) ? ny-1 : j;

                gotmat = gotobst = 0;

                for (kk = k0; kk <= k1; kk++) {
                for (jj = j0; jj <= j1; jj++) {
                    if (!obdata[IJK2(ii,jj,kk)])
                        gotmat = 1;
                }}

                if (gotmat) {
                    vertex_index[IJK1(i,j,k)] = npoints++;
                    vector_push_back(vertijk, i);
                    vector_push_back(vertijk, j);
                    vector_push_back(vertijk, k);
                }

                j0 = j;
            }
            k0 = k;
        }
        ii = nx-1;
    }

    jj = 0;
    for (j = jmin; j <= jmax; j += ny) {
        k0 = 0;
        for (k = 0; k <= nz; k++) {
            k1 = (k == nz) ? k : nz-1;
            i0 = 0;
            for (i = 0; i <= nx; i++) {
                i1 = (i == nx) ? nx-1 : i;

                gotmat = gotobst = 0;

                for (kk = k0; kk <= k1; kk++) {
                for (ii = i0; ii <= i1; ii++) {
                    if (!obdata[IJK2(ii,jj,kk)])
                        gotmat = 1;
                }}

                if (gotmat) {
                    vertex_index[IJK1(i,j,k)] = npoints++;
                    vector_push_back(vertijk, i);
                    vector_push_back(vertijk, j);
                    vector_push_back(vertijk, k);
                }

                i0 = i;
            }
            k0 = k;
        }
        jj = ny-1;
    }

    kk = 0;
    for (k = kmin; k <= kmax; k += nz) {
        j0 = 0;
        for (j = 0; j <= ny; j++) {
            j1 = (j == ny) ? j : ny-1;
            i0 = 0;
            for (i = 0; i <= nx; i++) {
                i1 = (i == nx) ? nx-1 : i;

                gotmat = gotobst = 0;

                for (jj = j0; jj <= j1; jj++) {
                for (ii = i0; ii <= i1; ii++) {
                    if (!obdata[IJK2(ii,jj,kk)])
                        gotmat = 1;
                }}

                if (gotmat) {
                    vertex_index[IJK1(i,j,k)] = npoints++;
                    vector_push_back(vertijk, i);
                    vector_push_back(vertijk, j);
                    vector_push_back(vertijk, k);
                }

                i0 = i;
            }
            j0 = j;
        }
        kk = nz-1;
    }

    b->npoints = npoints;
    b->data = malloc(b->ndims * npoints * sizeof(float));
    vertex = (float*)b->data;
    ijk = vertijk->data;

    if (grid->datatype_out == SDF_DATATYPE_REAL8) {
        double *x = grid->grids[0];
        double *y = grid->grids[1];
        double *z = grid->grids[2];
        for (i = 0; i < npoints; i++) {
            *vertex++ = x[*ijk++];
            *vertex++ = y[*ijk++];
            *vertex++ = z[*ijk++];
        }
    } else {
        float *x = grid->grids[0];
        float *y = grid->grids[1];
        float *z = grid->grids[2];
        for (i = 0; i < npoints; i++) {
            *vertex++ = x[*ijk++];
            *vertex++ = y[*ijk++];
            *vertex++ = z[*ijk++];
        }
    }

    vector_free(vertijk);

    // Scan faces in x-direction
    ii = 0;
    for (i = imin; i <= imax; i += nx) {
        for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
            if (!obdata[IJK2(ii,j,k)]) {
                vector_push_back(faces, vertex_index[IJK1(i,j  ,k  )]);
                vector_push_back(faces, vertex_index[IJK1(i,j+1,k  )]);
                vector_push_back(faces, vertex_index[IJK1(i,j+1,k+1)]);
                vector_push_back(faces, vertex_index[IJK1(i,j  ,k+1)]);
                vector_push_back(boundary, IJK(ii,j,k));
            }
        }}
        ii = nx - 1;
    }

    // Scan faces in y-direction
    jj = 0;
    for (j = jmin; j <= jmax; j += ny) {
        for (k = 0; k < nz; k++) {
        for (i = 0; i < nx; i++) {
            if (!obdata[IJK2(i,jj,k)]) {
                vector_push_back(faces, vertex_index[IJK1(i  ,j,k  )]);
                vector_push_back(faces, vertex_index[IJK1(i+1,j,k  )]);
                vector_push_back(faces, vertex_index[IJK1(i+1,j,k+1)]);
                vector_push_back(faces, vertex_index[IJK1(i  ,j,k+1)]);
                vector_push_back(boundary, IJK(i,jj,k));
            }
        }}
        jj = ny - 1;
    }

    // Scan faces in z-direction
    kk = 0;
    for (k = kmin; k <= kmax; k += nz) {
        for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            if (!obdata[IJK2(i,j,kk)]) {
                vector_push_back(faces, vertex_index[IJK1(i  ,j  ,k)]);
                vector_push_back(faces, vertex_index[IJK1(i+1,j  ,k)]);
                vector_push_back(faces, vertex_index[IJK1(i+1,j+1,k)]);
                vector_push_back(faces, vertex_index[IJK1(i  ,j+1,k)]);
                vector_push_back(boundary, IJK(i,j,kk));
            }
        }}
        kk = nz - 1;
    }

    free(vertex_index);

    if (b->boundary_cells) free(b->boundary_cells);

    b->nfaces = boundary->size;
    vector_truncate(boundary);
    b->boundary_cells = boundary->data;
    free(boundary);

    vector_truncate(faces);
    b->node_list = faces->data;
    free(faces);

    return b;
}



sdf_block_t *sdf_callback_surface_mesh(sdf_file_t *h, sdf_block_t *b)
{
    int grp = b->nm + 1;
    sdf_block_t *grid = sdf_find_block_by_id(h, b->subblock->mesh_id);
    sdf_block_t *current_block = h->current_block;

    int i, j, k, i1, j1, k1;
    int nx, ny, nz;
    int gotmat, gotobst, npoints, face, which, cell, ocell;
    int *vertex_index, *obdata;
    float *vertex;

    vector_t *faces = vector_new();
    vector_t *boundary = vector_new();

    h->current_block = grid;
    sdf_read_data(h);
    h->current_block = b->subblock;
    sdf_read_data(h);
    h->current_block = current_block;

    obdata = (int *)b->subblock->data;
    nx = b->subblock->local_dims[0] - 2;
    ny = b->subblock->local_dims[1] - 2;
    nz = b->subblock->local_dims[2] - 2;

    npoints = 0;
    vertex_index = (int*)malloc((nx+1)*(ny+1)*(nz+1)*sizeof(int));

    for (k = 0; k <= nz; k++) {
    for (j = 0; j <= ny; j++) {
    for (i = 0; i <= nx; i++) {
        vertex_index[IJK1(i,j,k)] = -1;
        gotmat = gotobst = 0;
        for (k1 = k-1; k1 <= k; k1++) {
        for (j1 = j-1; j1 <= j; j1++) {
        for (i1 = i-1; i1 <= i; i1++) {
            if (obdata[IJK2(i1,j1,k1)] == grp)
                gotobst = 1;
            else if (obdata[IJK2(i1,j1,k1)] == 0)
                gotmat = 1;
        }}}

        if (gotmat && gotobst)
            vertex_index[IJK1(i,j,k)] = npoints++;
    }}}

    b->npoints = npoints;
    b->data = malloc(b->ndims * npoints * sizeof(float));
    vertex = (float*)b->data;

    if (grid->datatype_out == SDF_DATATYPE_REAL8) {
        double *x = grid->grids[0];
        double *y = grid->grids[1];
        double *z = grid->grids[2];
        for (k = 0; k <= nz; k++) {
        for (j = 0; j <= ny; j++) {
        for (i = 0; i <= nx; i++) {
            if (vertex_index[IJK1(i,j,k)] != -1) {
                *vertex++ = x[i];
                *vertex++ = y[j];
                *vertex++ = z[k];
            }
        }}}
    } else {
        float *x = grid->grids[0];
        float *y = grid->grids[1];
        float *z = grid->grids[2];
        for (k = 0; k <= nz; k++) {
        for (j = 0; j <= ny; j++) {
        for (i = 0; i <= nx; i++) {
            if (vertex_index[IJK1(i,j,k)] != -1) {
                *vertex++ = x[i];
                *vertex++ = y[j];
                *vertex++ = z[k];
            }
        }}}
    }

    // Scan faces in x-direction
    for (k = 0; k <  nz; k++) {
    for (j = 0; j <  ny; j++) {
    for (i = 0; i <= nx; i++) {
        ocell = obdata[IJK2(i-1,j,k)];
        cell  = obdata[IJK2(i  ,j,k)];
        which = cell / grp;

        if (i == 0)
            face = (ocell == grp && cell == 0);
        else if (i == nx)
            face = (ocell == 0 && cell == grp);
        else
            face = (ocell == grp && cell == 0)
                    || (ocell == 0 && cell == grp);

        if (face) {
            vector_push_back(faces, vertex_index[IJK1(i,j  ,k  )]);
            vector_push_back(faces, vertex_index[IJK1(i,j+1,k  )]);
            vector_push_back(faces, vertex_index[IJK1(i,j+1,k+1)]);
            vector_push_back(faces, vertex_index[IJK1(i,j  ,k+1)]);
            vector_push_back(boundary, IJK(i-which,j,k));
        }
    }}}

    // Scan faces in y-direction
    for (k = 0; k <  nz; k++) {
    for (i = 0; i <  nx; i++) {
    for (j = 0; j <= ny; j++) {
        ocell = obdata[IJK2(i,j-1,k)];
        cell  = obdata[IJK2(i,j  ,k)];
        which = cell / grp;

        if (j == 0)
            face = (ocell == grp && cell == 0);
        else if (j == ny)
            face = (ocell == 0 && cell == grp);
        else
            face = (ocell == grp && cell == 0)
                    || (ocell == 0 && cell == grp);

        if (face) {
            vector_push_back(faces, vertex_index[IJK1(i  ,j,k  )]);
            vector_push_back(faces, vertex_index[IJK1(i+1,j,k  )]);
            vector_push_back(faces, vertex_index[IJK1(i+1,j,k+1)]);
            vector_push_back(faces, vertex_index[IJK1(i  ,j,k+1)]);
            vector_push_back(boundary, IJK(i,j-which,k));
        }
    }}}

    // Scan faces in z-direction
    for (j = 0; j <  ny; j++) {
    for (i = 0; i <  nx; i++) {
    for (k = 0; k <= nz; k++) {
        ocell = obdata[IJK2(i,j,k-1)];
        cell  = obdata[IJK2(i,j,k  )];
        which = cell / grp;

        if (k == 0)
            face = (ocell == grp && cell == 0);
        else if (k == nz)
            face = (ocell == 0 && cell == grp);
        else
            face = (ocell == grp && cell == 0)
                    || (ocell == 0 && cell == grp);

        if (face) {
            vector_push_back(faces, vertex_index[IJK1(i  ,j  ,k)]);
            vector_push_back(faces, vertex_index[IJK1(i+1,j  ,k)]);
            vector_push_back(faces, vertex_index[IJK1(i+1,j+1,k)]);
            vector_push_back(faces, vertex_index[IJK1(i  ,j+1,k)]);
            vector_push_back(boundary, IJK(i,j,k-which));
        }
    }}}

    free(vertex_index);

    if (b->boundary_cells) free(b->boundary_cells);

    b->nfaces = boundary->size;
    vector_truncate(boundary);
    b->boundary_cells = boundary->data;
    free(boundary);

    vector_truncate(faces);
    b->node_list = faces->data;
    free(faces);

    return b;
}



sdf_block_t *sdf_callback_surface(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
    sdf_block_t *current_block = h->current_block;
    int i, idx, sz;
    int *indexes = mesh->boundary_cells;
    char *ptr, *dptr;

    if (!b->subblock->data) {
        h->current_block = b->subblock;
        sdf_read_data(h);
        h->current_block = current_block;
    }

    b->nlocal = mesh->nfaces;
    b->datatype_out = b->subblock->datatype_out;
    sz = b->type_size_out = b->subblock->type_size_out;
    ptr = b->data = malloc(b->nlocal * b->type_size_out);
    dptr = b->subblock->data;
    for (i=0; i < b->nlocal; i++) {
        idx = *indexes++;
        memcpy(ptr, dptr + idx * sz, sz);
        ptr += sz;
    }

    return b;
}



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



static void add_surface_variables(sdf_file_t *h, sdf_block_t *surf,
        sdf_block_t **append_ptr, sdf_block_t **append_tail_ptr,
        int *nappend_ptr)
{
    sdf_block_t *b, *next, *append, *append_tail;
    int len1, len2, nappend = *nappend_ptr;
    char *str, *name1, *name2;

    append = *append_ptr;
    append_tail = *append_tail_ptr;

    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;

        if ((b->blocktype != SDF_BLOCKTYPE_PLAIN_VARIABLE &&
                b->blocktype != SDF_BLOCKTYPE_PLAIN_DERIVED) ||
                b->dont_display || b->stagger != SDF_STAGGER_CELL_CENTRE)
                    continue;

        append->next = calloc(1, sizeof(sdf_block_t));

        nappend++;
        append = append_tail = append->next;
        append->next = NULL;

        name1 = surf->id;
        name2 = b->id;
        len1 = strlen(name1);
        len2 = strlen(name2);
        str = (char*)malloc(len1 + len2 + 2);
        memcpy(str, name1, len1);
        str[len1] = '/';
        memcpy(str+len1+1, name2, len2);
        str[len1+len2+1] = '\0';
        append->id = str;

        name1 = surf->name;
        name2 = b->name;
        len1 = strlen(name1);
        len2 = strlen(name2);
        str = (char*)malloc(len1 + len2 + 2);
        memcpy(str, name1, len1);
        str[len1] = '/';
        memcpy(str+len1+1, name2, len2);
        str[len1+len2+1] = '\0';
        append->name = str;

        append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
        append->datatype = b->datatype;
        append->datatype_out = b->datatype_out;
        append->ndims = b->ndims;
        append->mult = b->mult;
        append->stagger = b->stagger;
        append->subblock = b;
        memcpy(append->dims, b->dims, 3 * sizeof(b->dims[0]));

        SDF_SET_ENTRY_ID(append->units, b->units);
        SDF_SET_ENTRY_ID(append->mesh_id, surf->id);
        append->done_header = 1;
        // Hack to prevent storage being allocated for this variable.
        append->dont_allocate = 1;
        append->populate_data = sdf_callback_surface;
    }

    *nappend_ptr = nappend;
    *append_tail_ptr = append_tail;
    *append_ptr = append;
}



int sdf_add_derived_blocks(sdf_file_t *h)
{
    sdf_block_t *b, *next, *append, *append_head, *append_tail;
    sdf_block_t *mesh;
    sdf_block_t *current_block = h->current_block;
    int i, len1, len2, nappend = 0;
    char *str, *name1, *name2;
    char *grid_ids[] = { "x", "y", "z" };

    append = append_head = calloc(1, sizeof(sdf_block_t));

    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;

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
    }

    if (nappend) {
        h->tail->next = append_head->next;
        h->tail = append_tail;
        h->nblocks += nappend;
    }

    h->current_block = current_block;

    return 0;
}



int sdf_add_derived_blocks_final(sdf_file_t *h)
{
    sdf_block_t *b, *next, *append, *append_head, *append_tail;
    sdf_block_t *mesh, *vfm, *obst;
    sdf_block_t *current_block = h->current_block;
    int i, len1, len2, nappend = 0;
    char *str, *name1, *name2;
    char *boundary_names[] =
        { "", "_x_min", "_x_max", "_y_min", "_y_max", "_z_min", "_z_max" };

    append = append_head = calloc(1, sizeof(sdf_block_t));

    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;

        if (b->blocktype == SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP) {
            obst = sdf_find_block_by_id(h, b->obstacle_id);
            if (!obst) continue;

            vfm = sdf_find_block_by_id(h, b->vfm_id);
            if (!vfm) continue;

            vfm->subblock = b;
            b->subblock = obst;
            obst->ng = 1;
            obst->dont_display = 1;
            mesh = sdf_find_block_by_id(h, obst->mesh_id);

            for (i = 0; i < b->ndims; i++) {
                append->next = calloc(1, sizeof(sdf_block_t));

                nappend++;
                append = append_tail = append->next;
                append->next = NULL;

                append->blocktype = SDF_BLOCKTYPE_UNSTRUCTURED_MESH;
                append->ndims = obst->ndims;
                append->subblock = obst;
                append->nm = i;
                append->populate_data = sdf_callback_surface_mesh;

                name1 = "obmsh";
                name2 = b->material_names[i];
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 1);
                memcpy(str, name1, len1);
                memcpy(str+len1, name2, len2);
                str[len1+len2] = '\0';
                append->id = str;

                name1 = "Surface Values/";
                name2 = b->material_names[i];
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 1);
                memcpy(str, name1, len1);
                memcpy(str+len1, name2, len2);
                str[len1+len2] = '\0';
                append->name = str;

                add_surface_variables(h, append, &append, &append_tail,
                    &nappend);
            }

            for (i = 0; i < 7; i++) {
                append->next = calloc(1, sizeof(sdf_block_t));

                nappend++;
                append = append_tail = append->next;
                append->next = NULL;

                append->blocktype = SDF_BLOCKTYPE_UNSTRUCTURED_MESH;
                append->ndims = obst->ndims;
                append->subblock = obst;
                append->nm = i;
                append->populate_data = sdf_callback_boundary_mesh;

                name1 = "boundary";
                name2 = boundary_names[i];
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 1);
                memcpy(str, name1, len1);
                memcpy(str+len1, name2, len2);
                str[len1+len2] = '\0';
                append->id = str;

                name1 = "Surface Values/boundary";
                name2 = boundary_names[i];
                len1 = strlen(name1);
                len2 = strlen(name2);
                str = (char*)malloc(len1 + len2 + 1);
                memcpy(str, name1, len1);
                memcpy(str+len1, name2, len2);
                str[len1+len2] = '\0';
                append->name = str;

                add_surface_variables(h, append, &append, &append_tail,
                    &nappend);
            }
        }
    }

    if (nappend) {
        h->tail->next = append_head->next;
        h->tail = append_tail;
        h->nblocks += nappend;
    }

    h->current_block = current_block;

    return 0;
}
