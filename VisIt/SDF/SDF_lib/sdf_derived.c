#include <stdlib.h>
#include "sdf.h"
#include "sdf_vector_type.h"
#include "sdf_list_type.h"

#define IJK(i,j,k) ((i) + nx * ((j) + ny * (k)))
#define IJK1(i,j,k) ((i) + (nx+1) * ((j) + (ny+1) * (k)))
#define IJK2(i,j,k) ((i)+1 + (nx+2) * ((j)+1 + (ny+2) * ((k)+1)))


static char *strcat_alloc(char *base, char *sfx)
{
    int len1 = strlen(base);
    int len2 = strlen(sfx);
    char *str = malloc(len1+len2+2);
    memcpy(str, base, len1);
    str[len1] = '/';
    memcpy(str+len1+1, sfx, len2); \
    str[len1+len2+1] = '\0';
    return str;
}



sdf_block_t *sdf_callback_boundary_mesh(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *grid = sdf_find_block_by_id(h, b->subblock->mesh_id);
    sdf_block_t *current_block = h->current_block;

    int i, j, k, ii, jj, kk, i0, j0, k0, i1, j1, k1;
    int imin, imax, jmin, jmax, kmin, kmax;
    int nx, ny, nz;
    int gotmat, gotobst, nelements;
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

    nelements = 0;
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
                    vertex_index[IJK1(i,j,k)] = nelements++;
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
                    vertex_index[IJK1(i,j,k)] = nelements++;
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
                    vertex_index[IJK1(i,j,k)] = nelements++;
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

    b->dims[0] = b->nelements = nelements;
    b->data = malloc(b->ndims * nelements * sizeof(float));
    vertex = (float*)b->data;
    ijk = vertijk->data;

    if (grid->datatype_out == SDF_DATATYPE_REAL8) {
        double *x = grid->grids[0];
        double *y = grid->grids[1];
        double *z = grid->grids[2];
        for (i = 0; i < nelements; i++) {
            *vertex++ = x[*ijk++];
            *vertex++ = y[*ijk++];
            *vertex++ = z[*ijk++];
        }
    } else {
        float *x = grid->grids[0];
        float *y = grid->grids[1];
        float *z = grid->grids[2];
        for (i = 0; i < nelements; i++) {
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
    int gotmat, gotobst, nelements, face, which, cell, ocell;
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

    nelements = 0;
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
            vertex_index[IJK1(i,j,k)] = nelements++;
    }}}

    b->dims[0] = b->nelements = nelements;
    b->data = malloc(b->ndims * nelements * sizeof(float));
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

    b->nelements_local = mesh->nfaces;
    b->datatype_out = b->subblock->datatype_out;
    sz = SDF_TYPE_SIZES[b->datatype_out];
    ptr = b->data = malloc(b->nelements_local * sz);
    dptr = b->subblock->data;
    for (i=0; i < b->nelements_local; i++) {
        idx = *indexes++;
        memcpy(ptr, dptr + idx * sz, sz);
        ptr += sz;
    }

    return b;
}



sdf_block_t *sdf_callback_grid_component(sdf_file_t *h, sdf_block_t *b)
{
    sdf_block_t *mesh = sdf_find_block_by_id(h, b->mesh_id);
    sdf_block_t *current_block = h->current_block;

    if (!b->grids) b->grids = calloc(1, sizeof(float*));

    if (!mesh->done_data) {
        h->current_block = mesh;
        sdf_read_data(h);
        h->current_block = current_block;
    }

    b->data = b->grids[0] = mesh->grids[b->nm];
    b->datatype_out = mesh->datatype_out;
    if (b->blocktype == SDF_BLOCKTYPE_POINT_DERIVED)
        b->local_dims[0] = b->nelements_local = mesh->nelements_local;
    else
        b->local_dims[0] = b->nelements_local = mesh->local_dims[b->nm];
    return b;
}



sdf_block_t *sdf_callback_face_grid(sdf_file_t *h, sdf_block_t *b)
{
    int i, n, sz;
    sdf_block_t *old = b->subblock;
    sdf_block_t *current_block = h->current_block;

    if (b->done_data) return b;

    if (!old->grids) {
        h->current_block = old;
        sdf_read_data(h);
        h->current_block = current_block;
    }

    memcpy(b->local_dims, old->local_dims, 3 * sizeof(int));

    b->grids = calloc(3, sizeof(float*));
    for (i = 0; i < 3; i++) {
        if (i != b->stagger) {
            sz = b->local_dims[i] * SDF_TYPE_SIZES[b->datatype_out];
            b->grids[i] = malloc(sz);
            memcpy(b->grids[i], old->grids[i], sz);
        }
    }

    n = b->stagger;
    b->local_dims[n] += 1;
    sz = b->local_dims[n] * SDF_TYPE_SIZES[b->datatype_out];
    b->grids[n] = malloc(sz);
    if (b->datatype_out == SDF_DATATYPE_REAL8) {
        double *oldx_ptr = (double*)old->grids[n];
        double *newx_ptr = (double*)b->grids[n];
        double x0, x1;
        x0 = *oldx_ptr++;
        x1 = *oldx_ptr;
        *newx_ptr++ = 0.5 * (3 * x0 - x1);

        for (i = 0; i < old->local_dims[n]-1; i++) {
            *newx_ptr++ = 0.5 * (x0 + *oldx_ptr);
            x0 = *oldx_ptr++;
        }

        x0 = *(oldx_ptr-1);
        x1 = *(oldx_ptr-2);
        *newx_ptr = 0.5 * (3 * x0 - x1);
    } else {
        float *oldx_ptr = (float*)old->grids[n];
        float *newx_ptr = (float*)b->grids[n];
        float x0, x1;
        x0 = *oldx_ptr++;
        x1 = *oldx_ptr;
        *newx_ptr++ = 0.5 * (3 * x0 - x1);

        for (i = 0; i < old->local_dims[n]-1; i++) {
            *newx_ptr++ = 0.5 * (x0 + *oldx_ptr);
            x0 = *oldx_ptr++;
        }

        x0 = *(oldx_ptr-1);
        x1 = *(oldx_ptr-2);
        *newx_ptr = 0.5 * (3 * x0 - x1);
    }
    b->done_data = 1;

    return b;
}



sdf_block_t *sdf_callback_cpu_mesh(sdf_file_t *h, sdf_block_t *b)
{
    int i, n, sz, np, nx;
    int i0, i1, idx;
    int *index;
    char *x, *xmesh;
    sdf_block_t *split = b->subblock;
    sdf_block_t *mesh = b->subblock2;
    sdf_block_t *current_block = h->current_block;
#ifdef PARALLEL
    void *buf;
#endif

    if (b->done_data) return b;

    if (!split->done_data) {
        h->current_block = split;
        sdf_read_data(h);
        h->current_block = current_block;
    }

    if (!mesh->done_data) {
        h->current_block = mesh;
        sdf_read_data(h);
        h->current_block = current_block;
    }
#ifdef PARALLEL
    memcpy(b->cpu_split, mesh->cpu_split, SDF_MAXDIMS * sizeof(int));
#endif

    b->datatype_out = mesh->datatype_out;
    b->grids = calloc(3, sizeof(float*));
    sz = SDF_TYPE_SIZES[b->datatype_out];

    if (split->geometry == 1) {
        index = split->data;
        np = 0;

        b->nelements_local = 1;
        for (n=0; n < b->ndims; n++) {
            xmesh = mesh->grids[n];

            nx = b->local_dims[n] = (int)b->dims[n];
            x = calloc(nx, sz);
#ifdef PARALLEL
            i0 = mesh->starts[n];
#else
            i0 = 0;
#endif
            i1 = mesh->local_dims[n] - 1;
            for (i = 1; i < nx-1; i++) {
                idx = index[np++] - i0;
                if (idx > i1) break;
                if (idx >= 0)
                    memcpy(x+sz*i, xmesh+sz*idx, sz);
            }
            if (i0 == 0)
                memcpy(x, xmesh, sz);
            idx = mesh->local_dims[n] - 1;
            if (mesh->dims[n] == mesh->local_dims[n] + i0)
                memcpy(x+sz*(nx-1), xmesh+sz*idx, sz);

#ifdef PARALLEL
            buf = malloc(nx*sz);
            MPI_Reduce(x, buf, nx*sz, MPI_BYTE, MPI_BOR, 0, h->comm);
            free(x);
            x = buf;

            if (h->rank)
                free(x);
            else
#endif
            b->grids[n] = x;
            b->nelements_local *= nx;
        }
        for (n=b->ndims; n < 3; n++) {
            b->local_dims[n] = 1;
            b->grids[n] = calloc(1, sz);
        }

#ifdef PARALLEL
        if (h->rank) {
            for (n=0; n < 3; n++) {
                b->local_dims[n] = 1;
                b->grids[n] = calloc(1, sz);
            }
            b->nelements_local = 1;
        }
#endif
    }

    b->done_data = 1;

    return b;
}



sdf_block_t *sdf_callback_current_cpu_mesh(sdf_file_t *h, sdf_block_t *b)
{
    int n, nx, sz, idx, i0 = 0, pmax = -2;
    char *x;
    sdf_block_t *mesh = b->subblock;
    sdf_block_t *current_block = h->current_block;
#ifdef PARALLEL
    int rem = h->rank;
    char *buf;
#endif

    if (b->done_data) return b;

    if (!mesh->done_data) {
        h->current_block = mesh;
        sdf_read_data(h);
        h->current_block = current_block;
    }
#ifdef PARALLEL
    memcpy(b->cpu_split, mesh->cpu_split, SDF_MAXDIMS * sizeof(int));
#endif

    b->datatype_out = mesh->datatype_out;
    b->grids = malloc(3 * sizeof(float*));
    sz = SDF_TYPE_SIZES[b->datatype_out];

    b->nelements_local = 1;
    for (n = 0; n < b->ndims; n++) {
#ifdef PARALLEL
        nx = mesh->cpu_split[n];
        i0 = rem%nx;
        rem = rem / nx;
        b->local_dims[n] = ++nx;
        pmax = mesh->proc_max[n];
#else
        nx = b->local_dims[n] = 2;
#endif
        idx = mesh->local_dims[n] - 1;
        x = calloc(nx, sz);
        memcpy(x+sz*i0, mesh->grids[n], sz);
        if (pmax < 0) memcpy(x+sz*(i0+1), (char*)mesh->grids[n]+sz*idx, sz);
        b->nelements_local *= nx;
#ifdef PARALLEL
        buf = malloc(nx*sz);
        MPI_Reduce(x, buf, nx*sz, MPI_BYTE, MPI_BOR, 0, h->comm);
        free(x);
        x = buf;

        if (h->rank)
            free(x);
        else
#endif
        b->grids[n] = x;
    }

    if (h->rank) {
        b->nelements_local = 1;
        nx = 0;
    } else
        nx = b->ndims;

    for (n = nx; n < 3; n++) {
        b->local_dims[n] = 1;
        b->grids[n] = calloc(1, sz);
    }

    b->subblock->nelements_local = b->nelements_local;
    b->done_data = 1;

    return b;
}



sdf_block_t *sdf_callback_cpu_data(sdf_file_t *h, sdf_block_t *b)
{
    int n, *var = b->data;

    if (b->done_data) return b;

    for (n=0; n < b->nelements_local; n++) var[n] = n;

    b->done_data = 1;

    return b;
}



sdf_block_t *sdf_callback_station_time(sdf_file_t *h, sdf_block_t *b)
{
    int i, n, idx, sz, data_offset0, varoffset;
    float *r4, dt4, time4;
    double *r8, dt8, time8;
    char *ptr;
    sdf_block_t *block, *current_block, *mesh = b->subblock;

    if (b->done_data) return b;

    if (mesh) {
        if (!mesh->done_data) {
            current_block = h->current_block;
            h->current_block = mesh;
            sdf_read_data(h);
            h->current_block = current_block;
        }

        if (b->datatype_out == SDF_DATATYPE_REAL4)
            b->data = (float*)mesh->data + b->opt;
        else
            b->data = (double*)mesh->data + b->opt;

        if (!b->grids) b->grids = malloc(sizeof(void**));
        b->grids[0] = b->data;

        b->done_data = 1;

        return b;
    }

    sz = SDF_TYPE_SIZES[b->datatype_out];
    if (!b->data) b->data = malloc(b->nelements_local * sz);

    // Find first station block
    block = h->blocklist;
    for (i = 0; i < h->nblocks; i++) {
        if (block->blocktype == SDF_BLOCKTYPE_STATION) break;
        block = block->next;
    }

    // If there is a fixed time increment then generate the data, otherwise
    // read it from file.
    if (block->time_increment > 0.0) {
        if (b->datatype_out == SDF_DATATYPE_REAL4) {
            r4 = (float*)b->data;
            dt4 = block->time_increment;
            time4 = block->time;
            for (n = 0; n < b->nelements; n++) {
                *r4++ = time4;
                time4 += dt4;
            }
        } else {
            r8 = (double*)b->data;
            dt8 = block->time_increment;
            time8 = block->time;
            for (n = 0; n < b->nelements; n++) {
                *r8++ = time8;
                time8 += dt8;
            }
        }
    } else {
        ptr = b->data;
        idx = i;
        data_offset0 = 0;
        varoffset = 0;
        if (block->step_increment <= 0) {
            data_offset0 += SDF_TYPE_SIZES[block->variable_types[varoffset]];
            varoffset++;
        }

        for (i = idx; i < h->nblocks; i++) {
            if (block->blocktype == SDF_BLOCKTYPE_STATION) {
                h->current_location = block->data_location + data_offset0;

                for (n = 0; n < block->nelements; n++) {
                    sdf_seek(h);
                    sdf_read_bytes(h, ptr, sz);
                    ptr += sz;
                    h->current_location += block->type_size;
                }
            }
            block = block->next;
        }
    }

    b->done_data = 1;

    return b;
}



sdf_block_t *sdf_callback_station(sdf_file_t *h, sdf_block_t *b)
{
    int i, j, k, sz, len, vidx, varoffset;
    int data_offset0, data_offset;
    list_t *station_blocks;
    char *varid, *ptr = b->data;
    sdf_block_t *block;

    if (b->done_data) return b;

    // Build list of all station blocks in the file
    list_init(&station_blocks);
    block = h->blocklist;
    for (i = 0; i < h->nblocks; i++) {
        if (block->blocktype == SDF_BLOCKTYPE_STATION)
            list_append(station_blocks, block);
        block = block->next;
    }

    len = strlen(b->station_id);
    varid = &b->id[len+1];
    sz = SDF_TYPE_SIZES[b->datatype_out];
    if (!b->data) b->data = malloc(b->nelements_local * sz);
    ptr = b->data;

    // Find station blocks containing this station id.
    block = list_start(station_blocks);
    varoffset = 0;
    data_offset0 = 0;
    if (block->step_increment <= 0) {
        data_offset0 += SDF_TYPE_SIZES[block->variable_types[varoffset]];
        varoffset++;
    }
    if (block->time_increment <= 0.0) {
        data_offset0 += SDF_TYPE_SIZES[block->variable_types[varoffset]];
        varoffset++;
    }

    for (i = 0; i < station_blocks->count; i++) {
        for (j = 0; j < block->nstations; j++) {
            if (strncmp(block->station_ids[j], b->station_id, SDF_ID_LENGTH))
                continue;
            if (block->station_move[j] < 0) continue;

            // Found station block. Now find the variable.
            vidx = varoffset;
            for (k = 0; k < j; k++) vidx += block->station_nvars[k];

            for (k = 0; k < block->station_nvars[j]; k++, vidx++) {
                if (strncmp(block->variable_ids[vidx], varid,
                    SDF_ID_LENGTH) == 0) break;
            }

            // Now read the data.
            data_offset = data_offset0;
            for (k = 0; k < vidx; k++)
                data_offset += SDF_TYPE_SIZES[block->variable_types[k]];

            h->current_location = block->data_location + data_offset;
            for (k = 0; k < block->nelements; k++) {
                sdf_seek(h);
                sdf_read_bytes(h, ptr, sz);
                ptr += sz;
                h->current_location += block->type_size;
            }
        }
        block = list_next(station_blocks);
    }

    list_destroy(&station_blocks);

    b->done_data = 1;

    return b;
}



static void add_surface_variables(sdf_file_t *h, sdf_block_t *surf,
        sdf_block_t **append_ptr, sdf_block_t **append_tail_ptr,
        int *nappend_ptr)
{
    sdf_block_t *b, *next, *append, *append_tail;
    int nappend = *nappend_ptr;
    size_t len1, len2;
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



static void add_station_variables(sdf_file_t *h, sdf_block_t **append,
        sdf_block_t **append_tail, int *nappend)
{
    sdf_block_t *b, *sb, *global_mesh, *mesh, *station = *append;
    char *meshid = "global_station/time";
    char *meshname = "Station/Time";
    int varoffset, mesh_datatype, var, i, j, n, is, nelements, nsofar;
    int *nelements_array;
    list_t *station_blocks;

    varoffset = 0;
    mesh_datatype = SDF_DATATYPE_REAL8;
    if (station->step_increment <= 0) varoffset++;

    if (station->time_increment <= 0.0) {
        varoffset++;
        mesh_datatype = station->variable_types[0];
    }

    // Build list of all station blocks in the file
    list_init(&station_blocks);
    b = h->blocklist;
    for (i = 0; i < h->nblocks; i++) {
        if (b->blocktype == SDF_BLOCKTYPE_STATION)
            list_append(station_blocks, b);
        b = b->next;
    }

    // Build list of nelements for each station block in the file
    nelements_array = calloc(station_blocks->count, sizeof(int));
    b = list_start(station_blocks);
    for (i = 0; i < station_blocks->count; i++) {
        nelements_array[i] = b->nelements;
        b = list_next(station_blocks);
    }

    // Reverse step over the nelements_array and fill with a cumulative
    // sum. nelements_array will then contain the number of time entries
    // for the current station block and all which come after it in the file.
    nelements = 0;
    for (i = station_blocks->count-1; i >= 0; i--) {
        nelements += nelements_array[i];
        nelements_array[i] = nelements;
    }

    /* Add global time mesh */
    global_mesh = mesh = (*append)->next = calloc(1, sizeof(sdf_block_t));
    (*nappend)++;
    (*append) = (*append_tail) = (*append)->next;

    SDF_SET_ENTRY_ID(mesh->id, meshid);
    SDF_SET_ENTRY_ID(mesh->units, "s");
    SDF_SET_ENTRY_STRING(mesh->name, meshname);
    mesh->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
    mesh->ndims = 1;
    mesh->dim_units = calloc(mesh->ndims, sizeof(char*));
    mesh->dim_labels = calloc(mesh->ndims, sizeof(char*));
    SDF_SET_ENTRY_ID(mesh->dim_units[0], "s");
    SDF_SET_ENTRY_ID(mesh->dim_labels[0], "Time");
    mesh->populate_data = sdf_callback_station_time;
    mesh->nelements_local = mesh->nelements = station->nelements;
    mesh->dims[0] = mesh->nelements;
    mesh->datatype = mesh->datatype_out = mesh_datatype;

    /* Add per-station-block time meshes */
    nsofar = 0;
    b = list_start(station_blocks);
    for (i = 0; i < station_blocks->count; i++) {
        mesh = (*append)->next = calloc(1, sizeof(sdf_block_t));
        (*nappend)++;
        (*append) = (*append_tail) = (*append)->next;

        mesh->id = strcat_alloc(meshid, b->id);
        mesh->name = strcat_alloc(meshname, b->name);
        SDF_SET_ENTRY_ID(mesh->units, "s");
        mesh->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
        mesh->ndims = 1;
        mesh->dim_units = calloc(mesh->ndims, sizeof(char*));
        mesh->dim_labels = calloc(mesh->ndims, sizeof(char*));
        SDF_SET_ENTRY_ID(mesh->dim_units[0], "s");
        SDF_SET_ENTRY_ID(mesh->dim_labels[0], "Time");
        mesh->populate_data = sdf_callback_station_time;
        mesh->nelements_local = mesh->nelements = nelements_array[i];
        mesh->dims[0] = mesh->nelements;
        mesh->datatype = mesh->datatype_out = mesh_datatype;
        mesh->subblock = global_mesh;
        mesh->opt = nsofar;
        mesh->dont_own_data = 1;

        nsofar += b->nelements;
        b = list_next(station_blocks);
    }

    free(nelements_array);

    var = varoffset;
    for (i = 0; i < station->nstations; i++) {
        // Sum nelements for all station blocks containing this station
        nelements = 0;
        sb = list_start(station_blocks);
        for (is = 0; is < station_blocks->count; is++) {
            for (n = 0; n < sb->nstations; n++) {
                if (strncmp(sb->station_ids[n], station->station_ids[i],
                        SDF_ID_LENGTH) == 0) {
                    if (sb->station_move[n] >= 0)
                        nelements += sb->nelements;
                    break;
                }
            }
            sb = list_next(station_blocks);
        }
        for (j = 0; j < station->station_nvars[i]; j++) {
            b = (*append)->next = calloc(1, sizeof(sdf_block_t));
            (*nappend)++;
            (*append) = (*append_tail) = (*append)->next;

            b->id = strcat_alloc(station->station_ids[i],
                    station->variable_ids[var]);
            b->name = strcat_alloc(station->station_names[i],
                    station->material_names[var]);
            b->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;
            SDF_SET_ENTRY_ID(b->units, station->dim_units[var]);
            SDF_SET_ENTRY_ID(b->station_id, station->station_ids[i]);
            b->ndims = 1;
            b->populate_data = sdf_callback_station;
            b->datatype = b->datatype_out = station->variable_types[var];
            b->nelements_local = b->nelements = nelements;
            b->dims[0] = b->nelements;
            var++;

            // Find the first station block which contains this station id
            nsofar = 0;
            sb = list_start(station_blocks);
            for (is = 0; is < station_blocks->count; is++) {
                for (n = 0; n < sb->nstations; n++) {
                    if (strncmp(sb->station_ids[n], b->station_id,
                            SDF_ID_LENGTH) != 0) continue;
                    if (sb->station_move[n] < 0) continue;
                    b->mesh_id = strcat_alloc(meshid, sb->id);
                    break;
                }
                if (b->mesh_id) break;
                nsofar += sb->nelements;
                sb = list_next(station_blocks);
            }
            b->opt = nsofar;
        }
    }

    list_destroy(&station_blocks);
/*
    b = (*append)->next = calloc(1, sizeof(sdf_block_t));
    (*nappend)++;
    (*append) = (*append_tail) = (*append)->next;

    b->blocktype == SDF_BLOCKTYPE_STATION_DERIVED;
    SDF_SET_ENTRY_ID(b->id, "stat");
*/
}



static void add_global_station(sdf_file_t *h, sdf_block_t **append,
        sdf_block_t **append_tail, int *nappend)
{
    sdf_block_t *b, *next, *new;
    int nextra, nstat_total, nstat_max = 2;
    int i, j, m, n, varoffset, nelements;
    int nidx; // variable index for new block
    int bidx; // variable index for current station block
    int *extra_id;
    char found, *found_id, *ctmp;
    list_t *station_blocks;

    // Create the new derived block and add it to the list
    new = (*append)->next = calloc(1, sizeof(sdf_block_t));
    (*nappend)++;
    (*append) = (*append_tail) = (*append)->next;

    new->blocktype = SDF_BLOCKTYPE_STATION_DERIVED;
    SDF_SET_ENTRY_ID(new->id, "global_stations");
    SDF_SET_ENTRY_STRING(new->name, "Global Stations");

    nelements = 0;
    extra_id = calloc(nstat_max, sizeof(int));
    found_id = malloc(sizeof(char));
    list_init(&station_blocks);

    // Build list of globally unique station_ids and for each
    // station block, assign the station indexes into this global array
    next = h->blocklist;
    for (i = 0; i < h->nblocks; i++) {
        h->current_block = b = next;
        next = b->next;
        if (b->blocktype != SDF_BLOCKTYPE_STATION) continue;

        nelements += b->nelements;
        list_append(station_blocks, b);

        if (b->nstations > nstat_max) {
            nstat_max = b->nstations * 11 / 10 + 2;
            free(extra_id);
            extra_id = calloc(nstat_max, sizeof(int));
        }

        b->station_index = calloc(b->nstations, sizeof(int));
        memset(found_id, 0, new->nstations);
        nextra = 0;
        for (n = 0; n < b->nstations; n++) {
            found = 0;
            for (m = 0; m < new->nstations; m++) {
                if (found_id[m]) continue;
                if (strncmp(b->station_ids[n], new->station_ids[m],
                        SDF_ID_LENGTH)) continue;
                b->station_index[n] = m;
                found = 1;
                found_id[m] = 1;
                break;
            }
            if (!found) {
                extra_id[nextra] = n;
                b->station_index[n] = new->nstations + n;
                nextra++;
            }
        }

        if (nextra == 0) continue;

        nstat_total = new->nstations + nextra;

        if (new->nstations == 0) {
            new->station_ids = calloc(nstat_total, sizeof(char*));
        } else {
            ctmp = malloc(new->nstations);
            memcpy(ctmp, new->station_ids, new->nstations * sizeof(char*));
            free(new->station_ids);
            new->station_ids = calloc(nstat_total, sizeof(char*));
            memcpy(new->station_ids, ctmp, new->nstations * sizeof(char*));
            free(ctmp);
        }

        free(found_id);
        found_id = malloc(nstat_total * sizeof(char));

        m = new->nstations;
        for (n = 0; n < nextra; n++) {
            SDF_SET_ENTRY_ID(new->station_ids[m], b->station_ids[extra_id[n]]);
            m++;
        }

        new->nstations = nstat_total;
    }


    // Populate each station_id with its associated data
    if (new->nstations > 0) {
        // Assign the station block variables
        b = list_start(station_blocks);
        new->step = b->step;
        new->step_increment = b->step_increment;
        new->time = b->time;
        new->time_increment = b->time_increment;
        new->use_mult = b->use_mult;

        varoffset = 0;
        if (new->step_increment <= 0) varoffset++;
        if (new->time_increment <= 0.0) varoffset++;
        new->nvariables = varoffset;
        new->nelements_local = new->nelements = nelements;

        // Allocate the stations
        new->station_names = calloc(new->nstations, sizeof(char*));
        new->station_nvars = calloc(new->nstations, sizeof(int*));
        new->station_move = calloc(new->nstations, sizeof(int*));
        new->station_index = calloc(new->nstations, sizeof(int*));
        new->station_x = calloc(new->nstations, sizeof(double*));
        if (new->ndims > 1)
            new->station_y = calloc(new->nstations, sizeof(double*));
        if (new->ndims > 2)
            new->station_z = calloc(new->nstations, sizeof(double*));

        // Assign the stations
        nstat_total = new->nstations;
        memset(found_id, 0, nstat_total);
        next = list_start(station_blocks);
        for (i = 0; i < station_blocks->count; i++) {
            b = next;
            next = list_next(station_blocks);
            new->station_index[i] = i;
            for (n = 0; n < b->nstations; n++) {
                m = b->station_index[n];
                if (found_id[m]) continue;
                nstat_total--;
                found_id[m] = 1;
                SDF_SET_ENTRY_STRING(new->station_names[m],
                        b->station_names[n]);
                new->station_nvars[m] = b->station_nvars[n];
                new->station_move[m] = b->station_move[n];
                new->station_x[m] = b->station_x[n];
                if (new->ndims > 1) new->station_y[m] = b->station_y[n];
                if (new->ndims > 2) new->station_z[m] = b->station_z[n];
                new->nvariables += new->station_nvars[m];
                if (!nstat_total) break;
            }
            if (!nstat_total) break;
        }

        // Allocate the variables
        new->variable_ids = calloc(new->nvariables, sizeof(char*));
        new->dim_units = calloc(new->nvariables, sizeof(char*));
        new->material_names = calloc(new->nvariables, sizeof(char*));
        new->variable_types = calloc(new->nvariables, sizeof(int*));
        if (new->use_mult)
            new->dim_mults = calloc(new->nvariables, sizeof(double*));

        // Assign starting variables (time, step)
        b = list_start(station_blocks);
        i = 0;
        nidx = 0;
        for (i = 0; i < varoffset; i++) {
            SDF_SET_ENTRY_ID(new->variable_ids[nidx], b->variable_ids[i]);
            SDF_SET_ENTRY_ID(new->dim_units[nidx], b->dim_units[i]);
            SDF_SET_ENTRY_STRING(new->material_names[nidx],
                b->material_names[i]);
            new->variable_types[nidx] = b->variable_types[i];
            if (new->use_mult)
                new->dim_mults[nidx] = b->dim_mults[i];
            nidx++;
        }

        // Assign the remaining variables
        nstat_total = new->nstations;
        memset(found_id, 0, nstat_total);
        next = list_start(station_blocks);
        for (i = 0; i < station_blocks->count; i++) {
            b = next;
            next = list_next(station_blocks);
            for (n = 0; n < b->nstations; n++) {
                m = b->station_index[n];
                if (found_id[m]) continue;
                nstat_total--;
                found_id[m] = 1;
                bidx = varoffset;
                for (j = 0; j < n; j++)
                    bidx += b->station_nvars[j];
                for (j = 0; j < b->station_nvars[n]; j++) {
                    SDF_SET_ENTRY_ID(new->variable_ids[nidx],
                            b->variable_ids[bidx]);
                    SDF_SET_ENTRY_ID(new->dim_units[nidx], b->dim_units[bidx]);
                    SDF_SET_ENTRY_STRING(new->material_names[nidx],
                            b->material_names[bidx]);
                    new->variable_types[nidx] = b->variable_types[bidx];
                    if (new->use_mult)
                        new->dim_mults[nidx] = b->dim_mults[bidx];
                    bidx++; nidx++;
                }
                if (!nstat_total) break;
            }
            if (!nstat_total) break;
        }
    }

    list_destroy(&station_blocks);
    free(extra_id);
    free(found_id);

    add_station_variables(h, append, append_tail, nappend);
}



int sdf_add_derived_blocks(sdf_file_t *h)
{
    sdf_block_t *b, *next, *append, *append_head, *append_tail;
    sdf_block_t *mesh, *first_mesh = NULL;
    sdf_block_t *current_block = h->current_block;
    int i, n, nd, nappend = 0;
    size_t len1, len2;
    char *str, *name1, *name2;
    char *grid_ids[] = { "x", "y", "z" };

    append = append_head = calloc(1, sizeof(sdf_block_t));

    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;

        if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH
            || b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
            if (!first_mesh && b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH)
                first_mesh = b;

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
            // Find first grid mesh
            mesh = h->blocklist;
            while (mesh) {
                if (mesh->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) break;
                mesh = mesh->next;
            }

            // First add CPU mesh
            append->next = calloc(1, sizeof(sdf_block_t));

            nappend++;
            append = append_tail = append->next;

            len1 = strlen(b->id);
            str = malloc(len1+6);
            memcpy(str, "grid_", 5);
            memcpy(str+5, b->id, len1+1);
            append->id = str;

            len1 = strlen(b->name);
            str = malloc(len1+6);
            memcpy(str, "Grid/", 5);
            memcpy(str+5, b->name, len1+1);
            append->name = str;

            append->subblock = b;
            append->subblock2 = mesh;
            append->populate_data = sdf_callback_cpu_mesh;
            append->blocktype = SDF_BLOCKTYPE_PLAIN_MESH;

            append->ndims = b->ndims;
            for (i = 0; i < b->ndims; i++) append->dims[i] = b->dims[i] + 2;
            for (i = b->ndims; i < 3; i++) append->dims[i] = 1;

            nd = mesh->ndims;
            len1 = 2 * nd * sizeof(double*);
            append->extents = malloc(len1);
            memcpy(append->extents, mesh->extents, len1);
            append->dim_labels = calloc(nd, sizeof(char*));
            append->dim_units = calloc(nd, sizeof(char*));
            for (n = 0; n < nd; n++) {
                SDF_SET_ENTRY_STRING(append->dim_labels[n],
                    mesh->dim_labels[n]);
                SDF_SET_ENTRY_STRING(append->dim_units[n],
                    mesh->dim_units[n]);
            }
            mesh = append;

            // Now add CPU data block
            append->next = calloc(1, sizeof(sdf_block_t));

            nappend++;
            append = append_tail = append->next;
            append->next = NULL;

            // Rename original block so that we can use the original name
            // for plotting
            str = b->id;
            append->id = str;
            len1 = strlen(str);
            str = malloc(len1+7);
            memcpy(str, append->id, len1);
            memcpy(str+len1, "_orig", 6);
            b->id = str;

            str = b->name;
            append->name = str;
            len1 = strlen(str);
            str = malloc(len1+7);
            memcpy(str, append->name, len1);
            memcpy(str+len1, "_orig", 6);
            b->name = str;

            append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;

            SDF_SET_ENTRY_ID(append->mesh_id, mesh->id);
            SDF_SET_ENTRY_ID(append->units, "CPU");
            append->ndims = b->ndims;
            append->nelements_local = 1;
            for (i=0; i<b->ndims; i++) {
                append->dims[i] = b->dims[i] + 1;
                append->local_dims[i] = (int)append->dims[i];
                append->nelements_local *= append->local_dims[i];
            }
            for (i=b->ndims; i<3; i++)
                append->local_dims[i] = append->dims[i] = 1;
            append->n_ids = 1;
            append->variable_ids = calloc(append->n_ids, sizeof(char*));
            SDF_SET_ENTRY_ID(append->variable_ids[0], b->id);
            append->must_read = calloc(append->n_ids, sizeof(char*));
            append->must_read[0] = 1;
            append->populate_data = sdf_callback_cpu_data;
            append->datatype = append->datatype_out = SDF_DATATYPE_INTEGER4;
        }
    }

    if (first_mesh) {
        // First add CPU mesh
        append->next = calloc(1, sizeof(sdf_block_t));

        nappend++;
        append = append_tail = append->next;

        SDF_SET_ENTRY_ID(append->id, "grid_cpus_current");
        SDF_SET_ENTRY_STRING(append->name, "Grid/CPUs/Current rank");

        append->subblock = first_mesh;
        append->populate_data = sdf_callback_current_cpu_mesh;
        append->blocktype = SDF_BLOCKTYPE_PLAIN_MESH;

        nd = append->ndims = first_mesh->ndims;

        len1 = 2 * nd * sizeof(double*);
        append->extents = malloc(len1);
        memcpy(append->extents, first_mesh->extents, len1);
        append->dim_labels = calloc(nd, sizeof(char*));
        append->dim_units = calloc(nd, sizeof(char*));
        for (n = 0; n < nd; n++) {
            SDF_SET_ENTRY_STRING(append->dim_labels[n],
                first_mesh->dim_labels[n]);
            SDF_SET_ENTRY_STRING(append->dim_units[n],
                first_mesh->dim_units[n]);
        }
        mesh = append;

        // Now add CPU data block
        append->next = calloc(1, sizeof(sdf_block_t));

        nappend++;
        append = append_tail = append->next;
        append->next = NULL;

        SDF_SET_ENTRY_ID(append->id, "cpus_current");
        SDF_SET_ENTRY_ID(append->name, "CPUs/Current rank");

        append->blocktype = SDF_BLOCKTYPE_PLAIN_DERIVED;

        SDF_SET_ENTRY_ID(append->units, "CPU");
        SDF_SET_ENTRY_ID(append->mesh_id, mesh->id);
        append->ndims = mesh->ndims;
        append->subblock = mesh;
        append->populate_data = sdf_callback_cpu_data;
        append->datatype = append->datatype_out = SDF_DATATYPE_INTEGER4;

        mesh->subblock2 = append;
    }

    if (h->station_file)
        add_global_station(h, &append, &append_tail, &nappend);

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
    sdf_block_t *mesh, *old_mesh, *vfm, *obst;
    sdf_block_t *current_block = h->current_block;
    int i, n, stagger, dont_add, nappend = 0;
    size_t nd, len1, len2;
    char *str, *name1, *name2;
    char *boundary_names[] =
        { "", "_x_min", "_x_max", "_y_min", "_y_max", "_z_min", "_z_max" };
    char *face_grid_ids[] = { "xface", "yface", "zface" };

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
            mesh = sdf_find_block_by_id(h, obst->mesh_id);

            if (obst->ng == 0) {
                obst->ng = 1;
                obst->dont_display = 1;
                obst->nelements_local = 1;
                for (i = 0; i < obst->ndims; i++) {
                    obst->dims[i] += 2 * obst->ng;
                    obst->local_dims[i] += 2 * obst->ng;
                    obst->nelements_local *= obst->local_dims[i];
                }
            }

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
        } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE
                || b->blocktype == SDF_BLOCKTYPE_PLAIN_DERIVED) {
            // Add grids for face-staggered variables
            for (i = 0; i < 3; i++) {
                stagger = 1<<i;
                if (b->stagger == stagger) {
                    old_mesh = sdf_find_block_by_id(h, b->mesh_id);
                    // For now, only add the mesh if variables are the correct
                    // dimensions
                    dont_add = 0;
                    for (n = 0; n < b->ndims; n++) {
                        nd = (size_t)old_mesh->dims[n];
                        if (n == i) nd += b->ng;
                        if (b->dims[n] != nd) dont_add = 1;
                    }
                    if (dont_add) continue;

                    name1 = b->mesh_id;
                    name2 = face_grid_ids[i];
                    len1 = strlen(name1);
                    len2 = strlen(name2);
                    str = (char*)malloc(len1 + len2 + 2);
                    memcpy(str, name1, len1);
                    str[len1] = '/';
                    memcpy(str+len1+1, name2, len2);
                    str[len1+len2+1] = '\0';
                    free(b->mesh_id);
                    b->mesh_id = str;

                    // Add new face grid if it doesn't exist
                    mesh = sdf_find_block_by_id(h, b->mesh_id);
                    if (!mesh) {
                        append->next = calloc(1, sizeof(sdf_block_t));

                        nappend++;
                        append = append_tail = append->next;

                        memcpy(append, old_mesh, sizeof(sdf_block_t));
                        append->next = NULL;

                        str = (char*)malloc(len1 + len2 + 2);
                        memcpy(str, b->mesh_id, len1+len2+2);
                        append->id = str;

                        name1 = old_mesh->name;
                        len1 = strlen(name1);
                        str = (char*)malloc(len1 + len2 + 2);
                        memcpy(str, name1, len1);
                        str[len1] = '/';
                        memcpy(str+len1+1, name2, len2);
                        str[len1+len2+1] = '\0';
                        append->name = str;

                        append->stagger = i;
                        append->subblock = old_mesh;
                        append->populate_data = sdf_callback_face_grid;

                        append->dims[i] += 1;
                        append->dim_mults = NULL;
                        append->extents = NULL;

                        nd = old_mesh->ndims;
                        len1 = 2 * nd * sizeof(double*);
                        append->extents = malloc(len1);
                        memcpy(append->extents, old_mesh->extents, len1);
                        append->dim_labels = calloc(nd, sizeof(char*));
                        append->dim_units = calloc(nd, sizeof(char*));
                        for (n = 0; n < nd; n++) {
                            SDF_SET_ENTRY_STRING(append->dim_labels[n],
                                old_mesh->dim_labels[n]);
                            SDF_SET_ENTRY_STRING(append->dim_units[n],
                                old_mesh->dim_units[n]);
                        }
                    }
                }
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
