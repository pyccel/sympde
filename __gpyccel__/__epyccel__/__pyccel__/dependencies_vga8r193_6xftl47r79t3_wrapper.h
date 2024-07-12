#ifndef DEPENDENCIES_VGA8R193_6XFTL47R79T3_WRAPPER_H
#define DEPENDENCIES_VGA8R193_6XFTL47R79T3_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef DEPENDENCIES_VGA8R193_6XFTL47R79T3_WRAPPER

void bind_c_lo_dot_vga8r193(void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pydependencies_vga8r193_6xftl47r79t3_API;


/*........................................*/
static int dependencies_vga8r193_6xftl47r79t3_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/Users/patricklagarrigue/psydac_workspace/sympde/__gpyccel__/__epyccel__"));
    Pydependencies_vga8r193_6xftl47r79t3_API = (void**)PyCapsule_Import("dependencies_vga8r193_6xftl47r79t3._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pydependencies_vga8r193_6xftl47r79t3_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // DEPENDENCIES_VGA8R193_6XFTL47R79T3_WRAPPER_H
