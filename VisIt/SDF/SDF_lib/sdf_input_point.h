#ifndef _SDF_INPUT_POINT_H_
#define _SDF_INPUT_POINT_H_

#ifdef __cplusplus
extern "C" {
#endif

int sdf_read_point_mesh(sdf_file_t *h);
int sdf_read_point_mesh_info(sdf_file_t *h);
int sdf_read_point_variable(sdf_file_t *h);
int sdf_read_point_variable_info(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
