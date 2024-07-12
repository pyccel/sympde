#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#define DEPENDENCIES_VGA8R193_6XFTL47R79T3_WRAPPER

#include "dependencies_vga8r193_6xftl47r79t3_wrapper.h"
#include <stdlib.h>
#include <stdint.h>
#include "ndarrays.h"


/*........................................*/


/*........................................*/

/*........................................*/
static PyObject* bind_c_lo_dot_vga8r193_wrapper(PyObject* self, PyObject* args, PyObject* kwargs)
{
    PyObject* mat00_obj;
    PyObject* mat01_obj;
    PyObject* mat10_obj;
    PyObject* mat11_obj;
    PyObject* x0_obj;
    PyObject* x1_obj;
    PyObject* out0_obj;
    PyObject* out1_obj;
    PyObject* s00_1_obj;
    PyObject* s00_2_obj;
    PyObject* s01_1_obj;
    PyObject* s01_2_obj;
    PyObject* s10_1_obj;
    PyObject* s10_2_obj;
    PyObject* s11_1_obj;
    PyObject* s11_2_obj;
    PyObject* n00_1_obj;
    PyObject* n00_2_obj;
    PyObject* n01_1_obj;
    PyObject* n01_2_obj;
    PyObject* n10_1_obj;
    PyObject* n10_2_obj;
    PyObject* n11_1_obj;
    PyObject* n11_2_obj;
    PyObject* ne00_1_obj;
    PyObject* ne00_2_obj;
    PyObject* ne01_1_obj;
    PyObject* ne01_2_obj;
    PyObject* ne10_1_obj;
    PyObject* ne10_2_obj;
    PyObject* ne11_1_obj;
    PyObject* ne11_2_obj;
    t_ndarray mat00 = {.shape = NULL};
    void* bound_mat00;
    int64_t bound_mat00_shape_1;
    int64_t bound_mat00_shape_2;
    int64_t bound_mat00_shape_3;
    int64_t bound_mat00_shape_4;
    int64_t bound_mat00_stride_1;
    int64_t bound_mat00_stride_2;
    int64_t bound_mat00_stride_3;
    int64_t bound_mat00_stride_4;
    t_ndarray mat01 = {.shape = NULL};
    void* bound_mat01;
    int64_t bound_mat01_shape_1;
    int64_t bound_mat01_shape_2;
    int64_t bound_mat01_shape_3;
    int64_t bound_mat01_shape_4;
    int64_t bound_mat01_stride_1;
    int64_t bound_mat01_stride_2;
    int64_t bound_mat01_stride_3;
    int64_t bound_mat01_stride_4;
    t_ndarray mat10 = {.shape = NULL};
    void* bound_mat10;
    int64_t bound_mat10_shape_1;
    int64_t bound_mat10_shape_2;
    int64_t bound_mat10_shape_3;
    int64_t bound_mat10_shape_4;
    int64_t bound_mat10_stride_1;
    int64_t bound_mat10_stride_2;
    int64_t bound_mat10_stride_3;
    int64_t bound_mat10_stride_4;
    t_ndarray mat11 = {.shape = NULL};
    void* bound_mat11;
    int64_t bound_mat11_shape_1;
    int64_t bound_mat11_shape_2;
    int64_t bound_mat11_shape_3;
    int64_t bound_mat11_shape_4;
    int64_t bound_mat11_stride_1;
    int64_t bound_mat11_stride_2;
    int64_t bound_mat11_stride_3;
    int64_t bound_mat11_stride_4;
    t_ndarray x0 = {.shape = NULL};
    void* bound_x0;
    int64_t bound_x0_shape_1;
    int64_t bound_x0_shape_2;
    int64_t bound_x0_stride_1;
    int64_t bound_x0_stride_2;
    t_ndarray x1 = {.shape = NULL};
    void* bound_x1;
    int64_t bound_x1_shape_1;
    int64_t bound_x1_shape_2;
    int64_t bound_x1_stride_1;
    int64_t bound_x1_stride_2;
    t_ndarray out0 = {.shape = NULL};
    void* bound_out0;
    int64_t bound_out0_shape_1;
    int64_t bound_out0_shape_2;
    int64_t bound_out0_stride_1;
    int64_t bound_out0_stride_2;
    t_ndarray out1 = {.shape = NULL};
    void* bound_out1;
    int64_t bound_out1_shape_1;
    int64_t bound_out1_shape_2;
    int64_t bound_out1_stride_1;
    int64_t bound_out1_stride_2;
    int64_t s00_1;
    int64_t s00_2;
    int64_t s01_1;
    int64_t s01_2;
    int64_t s10_1;
    int64_t s10_2;
    int64_t s11_1;
    int64_t s11_2;
    int64_t n00_1;
    int64_t n00_2;
    int64_t n01_1;
    int64_t n01_2;
    int64_t n10_1;
    int64_t n10_2;
    int64_t n11_1;
    int64_t n11_2;
    int64_t ne00_1;
    int64_t ne00_2;
    int64_t ne01_1;
    int64_t ne01_2;
    int64_t ne10_1;
    int64_t ne10_2;
    int64_t ne11_1;
    int64_t ne11_2;
    static char *kwlist[] = {
        "mat00",
        "mat01",
        "mat10",
        "mat11",
        "x0",
        "x1",
        "out0",
        "out1",
        "s00_1",
        "s00_2",
        "s01_1",
        "s01_2",
        "s10_1",
        "s10_2",
        "s11_1",
        "s11_2",
        "n00_1",
        "n00_2",
        "n01_1",
        "n01_2",
        "n10_1",
        "n10_2",
        "n11_1",
        "n11_2",
        "ne00_1",
        "ne00_2",
        "ne01_1",
        "ne01_2",
        "ne10_1",
        "ne10_2",
        "ne11_1",
        "ne11_2",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO", kwlist, &mat00_obj, &mat01_obj, &mat10_obj, &mat11_obj, &x0_obj, &x1_obj, &out0_obj, &out1_obj, &s00_1_obj, &s00_2_obj, &s01_1_obj, &s01_2_obj, &s10_1_obj, &s10_2_obj, &s11_1_obj, &s11_2_obj, &n00_1_obj, &n00_2_obj, &n01_1_obj, &n01_2_obj, &n10_1_obj, &n10_2_obj, &n11_1_obj, &n11_2_obj, &ne00_1_obj, &ne00_2_obj, &ne01_1_obj, &ne01_2_obj, &ne10_1_obj, &ne10_2_obj, &ne11_1_obj, &ne11_2_obj))
    {
        return NULL;
    }
    if (pyarray_check(mat00_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        mat00 = pyarray_to_ndarray(mat00_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument mat00");
        return NULL;
    }
    bound_mat00 = nd_data(&mat00);
    bound_mat00_shape_1 = nd_ndim(&mat00, INT64_C(0));
    bound_mat00_shape_2 = nd_ndim(&mat00, INT64_C(1));
    bound_mat00_shape_3 = nd_ndim(&mat00, INT64_C(2));
    bound_mat00_shape_4 = nd_ndim(&mat00, INT64_C(3));
    bound_mat00_stride_1 = nd_nstep_C(&mat00, INT64_C(0));
    bound_mat00_stride_2 = nd_nstep_C(&mat00, INT64_C(1));
    bound_mat00_stride_3 = nd_nstep_C(&mat00, INT64_C(2));
    bound_mat00_stride_4 = nd_nstep_C(&mat00, INT64_C(3));
    if (pyarray_check(mat01_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        mat01 = pyarray_to_ndarray(mat01_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument mat01");
        return NULL;
    }
    bound_mat01 = nd_data(&mat01);
    bound_mat01_shape_1 = nd_ndim(&mat01, INT64_C(0));
    bound_mat01_shape_2 = nd_ndim(&mat01, INT64_C(1));
    bound_mat01_shape_3 = nd_ndim(&mat01, INT64_C(2));
    bound_mat01_shape_4 = nd_ndim(&mat01, INT64_C(3));
    bound_mat01_stride_1 = nd_nstep_C(&mat01, INT64_C(0));
    bound_mat01_stride_2 = nd_nstep_C(&mat01, INT64_C(1));
    bound_mat01_stride_3 = nd_nstep_C(&mat01, INT64_C(2));
    bound_mat01_stride_4 = nd_nstep_C(&mat01, INT64_C(3));
    if (pyarray_check(mat10_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        mat10 = pyarray_to_ndarray(mat10_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument mat10");
        return NULL;
    }
    bound_mat10 = nd_data(&mat10);
    bound_mat10_shape_1 = nd_ndim(&mat10, INT64_C(0));
    bound_mat10_shape_2 = nd_ndim(&mat10, INT64_C(1));
    bound_mat10_shape_3 = nd_ndim(&mat10, INT64_C(2));
    bound_mat10_shape_4 = nd_ndim(&mat10, INT64_C(3));
    bound_mat10_stride_1 = nd_nstep_C(&mat10, INT64_C(0));
    bound_mat10_stride_2 = nd_nstep_C(&mat10, INT64_C(1));
    bound_mat10_stride_3 = nd_nstep_C(&mat10, INT64_C(2));
    bound_mat10_stride_4 = nd_nstep_C(&mat10, INT64_C(3));
    if (pyarray_check(mat11_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        mat11 = pyarray_to_ndarray(mat11_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument mat11");
        return NULL;
    }
    bound_mat11 = nd_data(&mat11);
    bound_mat11_shape_1 = nd_ndim(&mat11, INT64_C(0));
    bound_mat11_shape_2 = nd_ndim(&mat11, INT64_C(1));
    bound_mat11_shape_3 = nd_ndim(&mat11, INT64_C(2));
    bound_mat11_shape_4 = nd_ndim(&mat11, INT64_C(3));
    bound_mat11_stride_1 = nd_nstep_C(&mat11, INT64_C(0));
    bound_mat11_stride_2 = nd_nstep_C(&mat11, INT64_C(1));
    bound_mat11_stride_3 = nd_nstep_C(&mat11, INT64_C(2));
    bound_mat11_stride_4 = nd_nstep_C(&mat11, INT64_C(3));
    if (pyarray_check(x0_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        x0 = pyarray_to_ndarray(x0_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument x0");
        return NULL;
    }
    bound_x0 = nd_data(&x0);
    bound_x0_shape_1 = nd_ndim(&x0, INT64_C(0));
    bound_x0_shape_2 = nd_ndim(&x0, INT64_C(1));
    bound_x0_stride_1 = nd_nstep_C(&x0, INT64_C(0));
    bound_x0_stride_2 = nd_nstep_C(&x0, INT64_C(1));
    if (pyarray_check(x1_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        x1 = pyarray_to_ndarray(x1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument x1");
        return NULL;
    }
    bound_x1 = nd_data(&x1);
    bound_x1_shape_1 = nd_ndim(&x1, INT64_C(0));
    bound_x1_shape_2 = nd_ndim(&x1, INT64_C(1));
    bound_x1_stride_1 = nd_nstep_C(&x1, INT64_C(0));
    bound_x1_stride_2 = nd_nstep_C(&x1, INT64_C(1));
    if (pyarray_check(out0_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        out0 = pyarray_to_ndarray(out0_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument out0");
        return NULL;
    }
    bound_out0 = nd_data(&out0);
    bound_out0_shape_1 = nd_ndim(&out0, INT64_C(0));
    bound_out0_shape_2 = nd_ndim(&out0, INT64_C(1));
    bound_out0_stride_1 = nd_nstep_C(&out0, INT64_C(0));
    bound_out0_stride_2 = nd_nstep_C(&out0, INT64_C(1));
    if (pyarray_check(out1_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        out1 = pyarray_to_ndarray(out1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument out1");
        return NULL;
    }
    bound_out1 = nd_data(&out1);
    bound_out1_shape_1 = nd_ndim(&out1, INT64_C(0));
    bound_out1_shape_2 = nd_ndim(&out1, INT64_C(1));
    bound_out1_stride_1 = nd_nstep_C(&out1, INT64_C(0));
    bound_out1_stride_2 = nd_nstep_C(&out1, INT64_C(1));
    if (PyIs_Int64(s00_1_obj))
    {
        s00_1 = PyInt64_to_Int64(s00_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s00_1");
        return NULL;
    }
    if (PyIs_Int64(s00_2_obj))
    {
        s00_2 = PyInt64_to_Int64(s00_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s00_2");
        return NULL;
    }
    if (PyIs_Int64(s01_1_obj))
    {
        s01_1 = PyInt64_to_Int64(s01_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s01_1");
        return NULL;
    }
    if (PyIs_Int64(s01_2_obj))
    {
        s01_2 = PyInt64_to_Int64(s01_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s01_2");
        return NULL;
    }
    if (PyIs_Int64(s10_1_obj))
    {
        s10_1 = PyInt64_to_Int64(s10_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s10_1");
        return NULL;
    }
    if (PyIs_Int64(s10_2_obj))
    {
        s10_2 = PyInt64_to_Int64(s10_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s10_2");
        return NULL;
    }
    if (PyIs_Int64(s11_1_obj))
    {
        s11_1 = PyInt64_to_Int64(s11_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s11_1");
        return NULL;
    }
    if (PyIs_Int64(s11_2_obj))
    {
        s11_2 = PyInt64_to_Int64(s11_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument s11_2");
        return NULL;
    }
    if (PyIs_Int64(n00_1_obj))
    {
        n00_1 = PyInt64_to_Int64(n00_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n00_1");
        return NULL;
    }
    if (PyIs_Int64(n00_2_obj))
    {
        n00_2 = PyInt64_to_Int64(n00_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n00_2");
        return NULL;
    }
    if (PyIs_Int64(n01_1_obj))
    {
        n01_1 = PyInt64_to_Int64(n01_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n01_1");
        return NULL;
    }
    if (PyIs_Int64(n01_2_obj))
    {
        n01_2 = PyInt64_to_Int64(n01_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n01_2");
        return NULL;
    }
    if (PyIs_Int64(n10_1_obj))
    {
        n10_1 = PyInt64_to_Int64(n10_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n10_1");
        return NULL;
    }
    if (PyIs_Int64(n10_2_obj))
    {
        n10_2 = PyInt64_to_Int64(n10_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n10_2");
        return NULL;
    }
    if (PyIs_Int64(n11_1_obj))
    {
        n11_1 = PyInt64_to_Int64(n11_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n11_1");
        return NULL;
    }
    if (PyIs_Int64(n11_2_obj))
    {
        n11_2 = PyInt64_to_Int64(n11_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n11_2");
        return NULL;
    }
    if (PyIs_Int64(ne00_1_obj))
    {
        ne00_1 = PyInt64_to_Int64(ne00_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne00_1");
        return NULL;
    }
    if (PyIs_Int64(ne00_2_obj))
    {
        ne00_2 = PyInt64_to_Int64(ne00_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne00_2");
        return NULL;
    }
    if (PyIs_Int64(ne01_1_obj))
    {
        ne01_1 = PyInt64_to_Int64(ne01_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne01_1");
        return NULL;
    }
    if (PyIs_Int64(ne01_2_obj))
    {
        ne01_2 = PyInt64_to_Int64(ne01_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne01_2");
        return NULL;
    }
    if (PyIs_Int64(ne10_1_obj))
    {
        ne10_1 = PyInt64_to_Int64(ne10_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne10_1");
        return NULL;
    }
    if (PyIs_Int64(ne10_2_obj))
    {
        ne10_2 = PyInt64_to_Int64(ne10_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne10_2");
        return NULL;
    }
    if (PyIs_Int64(ne11_1_obj))
    {
        ne11_1 = PyInt64_to_Int64(ne11_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne11_1");
        return NULL;
    }
    if (PyIs_Int64(ne11_2_obj))
    {
        ne11_2 = PyInt64_to_Int64(ne11_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne11_2");
        return NULL;
    }
    bind_c_lo_dot_vga8r193(bound_mat00, bound_mat00_shape_1, bound_mat00_shape_2, bound_mat00_shape_3, bound_mat00_shape_4, bound_mat00_stride_1, bound_mat00_stride_2, bound_mat00_stride_3, bound_mat00_stride_4, bound_mat01, bound_mat01_shape_1, bound_mat01_shape_2, bound_mat01_shape_3, bound_mat01_shape_4, bound_mat01_stride_1, bound_mat01_stride_2, bound_mat01_stride_3, bound_mat01_stride_4, bound_mat10, bound_mat10_shape_1, bound_mat10_shape_2, bound_mat10_shape_3, bound_mat10_shape_4, bound_mat10_stride_1, bound_mat10_stride_2, bound_mat10_stride_3, bound_mat10_stride_4, bound_mat11, bound_mat11_shape_1, bound_mat11_shape_2, bound_mat11_shape_3, bound_mat11_shape_4, bound_mat11_stride_1, bound_mat11_stride_2, bound_mat11_stride_3, bound_mat11_stride_4, bound_x0, bound_x0_shape_1, bound_x0_shape_2, bound_x0_stride_1, bound_x0_stride_2, bound_x1, bound_x1_shape_1, bound_x1_shape_2, bound_x1_stride_1, bound_x1_stride_2, bound_out0, bound_out0_shape_1, bound_out0_shape_2, bound_out0_stride_1, bound_out0_stride_2, bound_out1, bound_out1_shape_1, bound_out1_shape_2, bound_out1_stride_1, bound_out1_stride_2, s00_1, s00_2, s01_1, s01_2, s10_1, s10_2, s11_1, s11_2, n00_1, n00_2, n01_1, n01_2, n10_1, n10_2, n11_1, n11_2, ne00_1, ne00_2, ne01_1, ne01_2, ne10_1, ne10_2, ne11_1, ne11_2);
    free_pointer(&mat00);
    free_pointer(&mat01);
    free_pointer(&mat10);
    free_pointer(&mat11);
    free_pointer(&x0);
    free_pointer(&x1);
    free_pointer(&out0);
    free_pointer(&out1);
    Py_INCREF(Py_None);
    return Py_None;
}
/*........................................*/

/*........................................*/

static PyMethodDef dependencies_vga8r193_6xftl47r79t3_methods[] = {
    {
        "lo_dot_vga8r193",
        (PyCFunction)bind_c_lo_dot_vga8r193_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef dependencies_vga8r193_6xftl47r79t3_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "dependencies_vga8r193_6xftl47r79t3",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    dependencies_vga8r193_6xftl47r79t3_methods,
};

/*........................................*/

PyMODINIT_FUNC PyInit_dependencies_vga8r193_6xftl47r79t3(void)
{
    PyObject* mod;
    static void* Pydependencies_vga8r193_6xftl47r79t3_API[0];
    PyObject* c_api_object_0001;
    mod = PyModule_Create(&dependencies_vga8r193_6xftl47r79t3_module);
    if (mod == NULL)
    {
        return NULL;
    }
    c_api_object_0001 = PyCapsule_New((void *)Pydependencies_vga8r193_6xftl47r79t3_API, "dependencies_vga8r193_6xftl47r79t3._C_API", NULL);
    if (PyModule_AddObject(mod, "_C_API", c_api_object_0001) < INT64_C(0))
    {
        Py_DECREF(mod);
        return NULL;
    }
    Py_INCREF(c_api_object_0001);
    import_array();
    return mod;
}
