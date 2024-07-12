#ifndef DEPENDENCIES_WTI60KR7_6FEYKYPLISIF_WRAPPER_H
#define DEPENDENCIES_WTI60KR7_6FEYKYPLISIF_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef DEPENDENCIES_WTI60KR7_6FEYKYPLISIF_WRAPPER

void bind_c_assemble_matrix_wti60kr7(void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pydependencies_wti60kr7_6feykyplisif_API;


/*........................................*/
static int dependencies_wti60kr7_6feykyplisif_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/Users/patricklagarrigue/psydac_workspace/sympde/__gpyccel__/__epyccel__"));
    Pydependencies_wti60kr7_6feykyplisif_API = (void**)PyCapsule_Import("dependencies_wti60kr7_6feykyplisif._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pydependencies_wti60kr7_6feykyplisif_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // DEPENDENCIES_WTI60KR7_6FEYKYPLISIF_WRAPPER_H
