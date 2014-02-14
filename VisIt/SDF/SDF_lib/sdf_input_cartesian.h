#ifndef _SDF_INPUT_CARTESIAN_H_
#define _SDF_INPUT_CARTESIAN_H_

#ifdef __cplusplus
extern "C" {
#endif

int sdf_read_plain_mesh(sdf_file_t *h);
int sdf_read_plain_mesh_info(sdf_file_t *h);
int sdf_read_lagran_mesh(sdf_file_t *h);
int sdf_read_plain_variable(sdf_file_t *h);
int sdf_read_plain_variable_info(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
