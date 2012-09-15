#ifndef _SDF_EXTENSION_H_
#define _SDF_EXTENSION_H_

#include "sdf.h"

#define SDF_EXTENSION_VERSION  1
#define SDF_EXTENSION_REVISION 0

typedef struct derived_template_struct derived_template_t;
typedef struct sdf_extension_struct sdf_extension_t;

struct sdf_extension_struct {
    derived_template_t *derived_list;
    int (*read_blocklist)(sdf_extension_t *, sdf_file_t *h);
    int (*timestate_update)(sdf_extension_t *, sdf_file_t *h);
    void (*get_version)(sdf_extension_t *, int *major, int *minor);
    char *(*get_name)(sdf_extension_t *);
    char **(*preload)(sdf_extension_t *, sdf_file_t *);
};

typedef sdf_extension_t *sdf_extension_create_t(sdf_file_t *h);
typedef void sdf_extension_destroy_t(sdf_extension_t *);

#endif
