#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#define DEPENDENCIES_E0EL77NQ_MX2LVWH9FED2_WRAPPER

#include "dependencies_e0el77nq_mx2lvwh9fed2_wrapper.h"
#include <stdlib.h>
#include <stdint.h>
#include "ndarrays.h"


/*........................................*/


/*........................................*/

/*........................................*/
static PyObject* bind_c_assemble_matrix_e0el77nq_wrapper(PyObject* self, PyObject* args, PyObject* kwargs)
{
    PyObject* global_test_basis_v1_0_1_obj;
    PyObject* global_test_basis_v1_0_2_obj;
    PyObject* global_test_basis_v1_1_1_obj;
    PyObject* global_test_basis_v1_1_2_obj;
    PyObject* global_trial_basis_u1_0_1_obj;
    PyObject* global_trial_basis_u1_0_2_obj;
    PyObject* global_trial_basis_u1_1_1_obj;
    PyObject* global_trial_basis_u1_1_2_obj;
    PyObject* global_span_v1_0_1_obj;
    PyObject* global_span_v1_0_2_obj;
    PyObject* global_span_v1_1_1_obj;
    PyObject* global_span_v1_1_2_obj;
    PyObject* global_x1_obj;
    PyObject* global_x2_obj;
    PyObject* test_v1_0_p1_obj;
    PyObject* test_v1_0_p2_obj;
    PyObject* test_v1_1_p1_obj;
    PyObject* test_v1_1_p2_obj;
    PyObject* trial_u1_0_p1_obj;
    PyObject* trial_u1_0_p2_obj;
    PyObject* trial_u1_1_p1_obj;
    PyObject* trial_u1_1_p2_obj;
    PyObject* n_element_1_obj;
    PyObject* n_element_2_obj;
    PyObject* k1_obj;
    PyObject* k2_obj;
    PyObject* pad1_obj;
    PyObject* pad2_obj;
    PyObject* g_mat_u1_0_v1_0_e0el77nq_obj;
    PyObject* g_mat_u1_1_v1_0_e0el77nq_obj;
    PyObject* g_mat_u1_0_v1_1_e0el77nq_obj;
    PyObject* g_mat_u1_1_v1_1_e0el77nq_obj;
    t_ndarray global_test_basis_v1_0_1 = {.shape = NULL};
    void* bound_global_test_basis_v1_0_1;
    int64_t bound_global_test_basis_v1_0_1_shape_1;
    int64_t bound_global_test_basis_v1_0_1_shape_2;
    int64_t bound_global_test_basis_v1_0_1_shape_3;
    int64_t bound_global_test_basis_v1_0_1_shape_4;
    int64_t bound_global_test_basis_v1_0_1_stride_1;
    int64_t bound_global_test_basis_v1_0_1_stride_2;
    int64_t bound_global_test_basis_v1_0_1_stride_3;
    int64_t bound_global_test_basis_v1_0_1_stride_4;
    t_ndarray global_test_basis_v1_0_2 = {.shape = NULL};
    void* bound_global_test_basis_v1_0_2;
    int64_t bound_global_test_basis_v1_0_2_shape_1;
    int64_t bound_global_test_basis_v1_0_2_shape_2;
    int64_t bound_global_test_basis_v1_0_2_shape_3;
    int64_t bound_global_test_basis_v1_0_2_shape_4;
    int64_t bound_global_test_basis_v1_0_2_stride_1;
    int64_t bound_global_test_basis_v1_0_2_stride_2;
    int64_t bound_global_test_basis_v1_0_2_stride_3;
    int64_t bound_global_test_basis_v1_0_2_stride_4;
    t_ndarray global_test_basis_v1_1_1 = {.shape = NULL};
    void* bound_global_test_basis_v1_1_1;
    int64_t bound_global_test_basis_v1_1_1_shape_1;
    int64_t bound_global_test_basis_v1_1_1_shape_2;
    int64_t bound_global_test_basis_v1_1_1_shape_3;
    int64_t bound_global_test_basis_v1_1_1_shape_4;
    int64_t bound_global_test_basis_v1_1_1_stride_1;
    int64_t bound_global_test_basis_v1_1_1_stride_2;
    int64_t bound_global_test_basis_v1_1_1_stride_3;
    int64_t bound_global_test_basis_v1_1_1_stride_4;
    t_ndarray global_test_basis_v1_1_2 = {.shape = NULL};
    void* bound_global_test_basis_v1_1_2;
    int64_t bound_global_test_basis_v1_1_2_shape_1;
    int64_t bound_global_test_basis_v1_1_2_shape_2;
    int64_t bound_global_test_basis_v1_1_2_shape_3;
    int64_t bound_global_test_basis_v1_1_2_shape_4;
    int64_t bound_global_test_basis_v1_1_2_stride_1;
    int64_t bound_global_test_basis_v1_1_2_stride_2;
    int64_t bound_global_test_basis_v1_1_2_stride_3;
    int64_t bound_global_test_basis_v1_1_2_stride_4;
    t_ndarray global_trial_basis_u1_0_1 = {.shape = NULL};
    void* bound_global_trial_basis_u1_0_1;
    int64_t bound_global_trial_basis_u1_0_1_shape_1;
    int64_t bound_global_trial_basis_u1_0_1_shape_2;
    int64_t bound_global_trial_basis_u1_0_1_shape_3;
    int64_t bound_global_trial_basis_u1_0_1_shape_4;
    int64_t bound_global_trial_basis_u1_0_1_stride_1;
    int64_t bound_global_trial_basis_u1_0_1_stride_2;
    int64_t bound_global_trial_basis_u1_0_1_stride_3;
    int64_t bound_global_trial_basis_u1_0_1_stride_4;
    t_ndarray global_trial_basis_u1_0_2 = {.shape = NULL};
    void* bound_global_trial_basis_u1_0_2;
    int64_t bound_global_trial_basis_u1_0_2_shape_1;
    int64_t bound_global_trial_basis_u1_0_2_shape_2;
    int64_t bound_global_trial_basis_u1_0_2_shape_3;
    int64_t bound_global_trial_basis_u1_0_2_shape_4;
    int64_t bound_global_trial_basis_u1_0_2_stride_1;
    int64_t bound_global_trial_basis_u1_0_2_stride_2;
    int64_t bound_global_trial_basis_u1_0_2_stride_3;
    int64_t bound_global_trial_basis_u1_0_2_stride_4;
    t_ndarray global_trial_basis_u1_1_1 = {.shape = NULL};
    void* bound_global_trial_basis_u1_1_1;
    int64_t bound_global_trial_basis_u1_1_1_shape_1;
    int64_t bound_global_trial_basis_u1_1_1_shape_2;
    int64_t bound_global_trial_basis_u1_1_1_shape_3;
    int64_t bound_global_trial_basis_u1_1_1_shape_4;
    int64_t bound_global_trial_basis_u1_1_1_stride_1;
    int64_t bound_global_trial_basis_u1_1_1_stride_2;
    int64_t bound_global_trial_basis_u1_1_1_stride_3;
    int64_t bound_global_trial_basis_u1_1_1_stride_4;
    t_ndarray global_trial_basis_u1_1_2 = {.shape = NULL};
    void* bound_global_trial_basis_u1_1_2;
    int64_t bound_global_trial_basis_u1_1_2_shape_1;
    int64_t bound_global_trial_basis_u1_1_2_shape_2;
    int64_t bound_global_trial_basis_u1_1_2_shape_3;
    int64_t bound_global_trial_basis_u1_1_2_shape_4;
    int64_t bound_global_trial_basis_u1_1_2_stride_1;
    int64_t bound_global_trial_basis_u1_1_2_stride_2;
    int64_t bound_global_trial_basis_u1_1_2_stride_3;
    int64_t bound_global_trial_basis_u1_1_2_stride_4;
    t_ndarray global_span_v1_0_1 = {.shape = NULL};
    void* bound_global_span_v1_0_1;
    int64_t bound_global_span_v1_0_1_shape_1;
    int64_t bound_global_span_v1_0_1_stride_1;
    t_ndarray global_span_v1_0_2 = {.shape = NULL};
    void* bound_global_span_v1_0_2;
    int64_t bound_global_span_v1_0_2_shape_1;
    int64_t bound_global_span_v1_0_2_stride_1;
    t_ndarray global_span_v1_1_1 = {.shape = NULL};
    void* bound_global_span_v1_1_1;
    int64_t bound_global_span_v1_1_1_shape_1;
    int64_t bound_global_span_v1_1_1_stride_1;
    t_ndarray global_span_v1_1_2 = {.shape = NULL};
    void* bound_global_span_v1_1_2;
    int64_t bound_global_span_v1_1_2_shape_1;
    int64_t bound_global_span_v1_1_2_stride_1;
    t_ndarray global_x1 = {.shape = NULL};
    void* bound_global_x1;
    int64_t bound_global_x1_shape_1;
    int64_t bound_global_x1_shape_2;
    int64_t bound_global_x1_stride_1;
    int64_t bound_global_x1_stride_2;
    t_ndarray global_x2 = {.shape = NULL};
    void* bound_global_x2;
    int64_t bound_global_x2_shape_1;
    int64_t bound_global_x2_shape_2;
    int64_t bound_global_x2_stride_1;
    int64_t bound_global_x2_stride_2;
    int64_t test_v1_0_p1;
    int64_t test_v1_0_p2;
    int64_t test_v1_1_p1;
    int64_t test_v1_1_p2;
    int64_t trial_u1_0_p1;
    int64_t trial_u1_0_p2;
    int64_t trial_u1_1_p1;
    int64_t trial_u1_1_p2;
    int64_t n_element_1;
    int64_t n_element_2;
    int64_t k1;
    int64_t k2;
    int64_t pad1;
    int64_t pad2;
    t_ndarray g_mat_u1_0_v1_0_e0el77nq = {.shape = NULL};
    void* bound_g_mat_u1_0_v1_0_e0el77nq;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_shape_1;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_shape_2;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_shape_3;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_shape_4;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_stride_1;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_stride_2;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_stride_3;
    int64_t bound_g_mat_u1_0_v1_0_e0el77nq_stride_4;
    t_ndarray g_mat_u1_1_v1_0_e0el77nq = {.shape = NULL};
    void* bound_g_mat_u1_1_v1_0_e0el77nq;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_shape_1;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_shape_2;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_shape_3;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_shape_4;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_stride_1;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_stride_2;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_stride_3;
    int64_t bound_g_mat_u1_1_v1_0_e0el77nq_stride_4;
    t_ndarray g_mat_u1_0_v1_1_e0el77nq = {.shape = NULL};
    void* bound_g_mat_u1_0_v1_1_e0el77nq;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_shape_1;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_shape_2;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_shape_3;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_shape_4;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_stride_1;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_stride_2;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_stride_3;
    int64_t bound_g_mat_u1_0_v1_1_e0el77nq_stride_4;
    t_ndarray g_mat_u1_1_v1_1_e0el77nq = {.shape = NULL};
    void* bound_g_mat_u1_1_v1_1_e0el77nq;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_shape_1;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_shape_2;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_shape_3;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_shape_4;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_stride_1;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_stride_2;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_stride_3;
    int64_t bound_g_mat_u1_1_v1_1_e0el77nq_stride_4;
    static char *kwlist[] = {
        "global_test_basis_v1_0_1",
        "global_test_basis_v1_0_2",
        "global_test_basis_v1_1_1",
        "global_test_basis_v1_1_2",
        "global_trial_basis_u1_0_1",
        "global_trial_basis_u1_0_2",
        "global_trial_basis_u1_1_1",
        "global_trial_basis_u1_1_2",
        "global_span_v1_0_1",
        "global_span_v1_0_2",
        "global_span_v1_1_1",
        "global_span_v1_1_2",
        "global_x1",
        "global_x2",
        "test_v1_0_p1",
        "test_v1_0_p2",
        "test_v1_1_p1",
        "test_v1_1_p2",
        "trial_u1_0_p1",
        "trial_u1_0_p2",
        "trial_u1_1_p1",
        "trial_u1_1_p2",
        "n_element_1",
        "n_element_2",
        "k1",
        "k2",
        "pad1",
        "pad2",
        "g_mat_u1_0_v1_0_e0el77nq",
        "g_mat_u1_1_v1_0_e0el77nq",
        "g_mat_u1_0_v1_1_e0el77nq",
        "g_mat_u1_1_v1_1_e0el77nq",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO", kwlist, &global_test_basis_v1_0_1_obj, &global_test_basis_v1_0_2_obj, &global_test_basis_v1_1_1_obj, &global_test_basis_v1_1_2_obj, &global_trial_basis_u1_0_1_obj, &global_trial_basis_u1_0_2_obj, &global_trial_basis_u1_1_1_obj, &global_trial_basis_u1_1_2_obj, &global_span_v1_0_1_obj, &global_span_v1_0_2_obj, &global_span_v1_1_1_obj, &global_span_v1_1_2_obj, &global_x1_obj, &global_x2_obj, &test_v1_0_p1_obj, &test_v1_0_p2_obj, &test_v1_1_p1_obj, &test_v1_1_p2_obj, &trial_u1_0_p1_obj, &trial_u1_0_p2_obj, &trial_u1_1_p1_obj, &trial_u1_1_p2_obj, &n_element_1_obj, &n_element_2_obj, &k1_obj, &k2_obj, &pad1_obj, &pad2_obj, &g_mat_u1_0_v1_0_e0el77nq_obj, &g_mat_u1_1_v1_0_e0el77nq_obj, &g_mat_u1_0_v1_1_e0el77nq_obj, &g_mat_u1_1_v1_1_e0el77nq_obj))
    {
        return NULL;
    }
    if (pyarray_check(global_test_basis_v1_0_1_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_test_basis_v1_0_1 = pyarray_to_ndarray(global_test_basis_v1_0_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_test_basis_v1_0_1");
        return NULL;
    }
    bound_global_test_basis_v1_0_1 = nd_data(&global_test_basis_v1_0_1);
    bound_global_test_basis_v1_0_1_shape_1 = nd_ndim(&global_test_basis_v1_0_1, INT64_C(0));
    bound_global_test_basis_v1_0_1_shape_2 = nd_ndim(&global_test_basis_v1_0_1, INT64_C(1));
    bound_global_test_basis_v1_0_1_shape_3 = nd_ndim(&global_test_basis_v1_0_1, INT64_C(2));
    bound_global_test_basis_v1_0_1_shape_4 = nd_ndim(&global_test_basis_v1_0_1, INT64_C(3));
    bound_global_test_basis_v1_0_1_stride_1 = nd_nstep_C(&global_test_basis_v1_0_1, INT64_C(0));
    bound_global_test_basis_v1_0_1_stride_2 = nd_nstep_C(&global_test_basis_v1_0_1, INT64_C(1));
    bound_global_test_basis_v1_0_1_stride_3 = nd_nstep_C(&global_test_basis_v1_0_1, INT64_C(2));
    bound_global_test_basis_v1_0_1_stride_4 = nd_nstep_C(&global_test_basis_v1_0_1, INT64_C(3));
    if (pyarray_check(global_test_basis_v1_0_2_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_test_basis_v1_0_2 = pyarray_to_ndarray(global_test_basis_v1_0_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_test_basis_v1_0_2");
        return NULL;
    }
    bound_global_test_basis_v1_0_2 = nd_data(&global_test_basis_v1_0_2);
    bound_global_test_basis_v1_0_2_shape_1 = nd_ndim(&global_test_basis_v1_0_2, INT64_C(0));
    bound_global_test_basis_v1_0_2_shape_2 = nd_ndim(&global_test_basis_v1_0_2, INT64_C(1));
    bound_global_test_basis_v1_0_2_shape_3 = nd_ndim(&global_test_basis_v1_0_2, INT64_C(2));
    bound_global_test_basis_v1_0_2_shape_4 = nd_ndim(&global_test_basis_v1_0_2, INT64_C(3));
    bound_global_test_basis_v1_0_2_stride_1 = nd_nstep_C(&global_test_basis_v1_0_2, INT64_C(0));
    bound_global_test_basis_v1_0_2_stride_2 = nd_nstep_C(&global_test_basis_v1_0_2, INT64_C(1));
    bound_global_test_basis_v1_0_2_stride_3 = nd_nstep_C(&global_test_basis_v1_0_2, INT64_C(2));
    bound_global_test_basis_v1_0_2_stride_4 = nd_nstep_C(&global_test_basis_v1_0_2, INT64_C(3));
    if (pyarray_check(global_test_basis_v1_1_1_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_test_basis_v1_1_1 = pyarray_to_ndarray(global_test_basis_v1_1_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_test_basis_v1_1_1");
        return NULL;
    }
    bound_global_test_basis_v1_1_1 = nd_data(&global_test_basis_v1_1_1);
    bound_global_test_basis_v1_1_1_shape_1 = nd_ndim(&global_test_basis_v1_1_1, INT64_C(0));
    bound_global_test_basis_v1_1_1_shape_2 = nd_ndim(&global_test_basis_v1_1_1, INT64_C(1));
    bound_global_test_basis_v1_1_1_shape_3 = nd_ndim(&global_test_basis_v1_1_1, INT64_C(2));
    bound_global_test_basis_v1_1_1_shape_4 = nd_ndim(&global_test_basis_v1_1_1, INT64_C(3));
    bound_global_test_basis_v1_1_1_stride_1 = nd_nstep_C(&global_test_basis_v1_1_1, INT64_C(0));
    bound_global_test_basis_v1_1_1_stride_2 = nd_nstep_C(&global_test_basis_v1_1_1, INT64_C(1));
    bound_global_test_basis_v1_1_1_stride_3 = nd_nstep_C(&global_test_basis_v1_1_1, INT64_C(2));
    bound_global_test_basis_v1_1_1_stride_4 = nd_nstep_C(&global_test_basis_v1_1_1, INT64_C(3));
    if (pyarray_check(global_test_basis_v1_1_2_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_test_basis_v1_1_2 = pyarray_to_ndarray(global_test_basis_v1_1_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_test_basis_v1_1_2");
        return NULL;
    }
    bound_global_test_basis_v1_1_2 = nd_data(&global_test_basis_v1_1_2);
    bound_global_test_basis_v1_1_2_shape_1 = nd_ndim(&global_test_basis_v1_1_2, INT64_C(0));
    bound_global_test_basis_v1_1_2_shape_2 = nd_ndim(&global_test_basis_v1_1_2, INT64_C(1));
    bound_global_test_basis_v1_1_2_shape_3 = nd_ndim(&global_test_basis_v1_1_2, INT64_C(2));
    bound_global_test_basis_v1_1_2_shape_4 = nd_ndim(&global_test_basis_v1_1_2, INT64_C(3));
    bound_global_test_basis_v1_1_2_stride_1 = nd_nstep_C(&global_test_basis_v1_1_2, INT64_C(0));
    bound_global_test_basis_v1_1_2_stride_2 = nd_nstep_C(&global_test_basis_v1_1_2, INT64_C(1));
    bound_global_test_basis_v1_1_2_stride_3 = nd_nstep_C(&global_test_basis_v1_1_2, INT64_C(2));
    bound_global_test_basis_v1_1_2_stride_4 = nd_nstep_C(&global_test_basis_v1_1_2, INT64_C(3));
    if (pyarray_check(global_trial_basis_u1_0_1_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_trial_basis_u1_0_1 = pyarray_to_ndarray(global_trial_basis_u1_0_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_trial_basis_u1_0_1");
        return NULL;
    }
    bound_global_trial_basis_u1_0_1 = nd_data(&global_trial_basis_u1_0_1);
    bound_global_trial_basis_u1_0_1_shape_1 = nd_ndim(&global_trial_basis_u1_0_1, INT64_C(0));
    bound_global_trial_basis_u1_0_1_shape_2 = nd_ndim(&global_trial_basis_u1_0_1, INT64_C(1));
    bound_global_trial_basis_u1_0_1_shape_3 = nd_ndim(&global_trial_basis_u1_0_1, INT64_C(2));
    bound_global_trial_basis_u1_0_1_shape_4 = nd_ndim(&global_trial_basis_u1_0_1, INT64_C(3));
    bound_global_trial_basis_u1_0_1_stride_1 = nd_nstep_C(&global_trial_basis_u1_0_1, INT64_C(0));
    bound_global_trial_basis_u1_0_1_stride_2 = nd_nstep_C(&global_trial_basis_u1_0_1, INT64_C(1));
    bound_global_trial_basis_u1_0_1_stride_3 = nd_nstep_C(&global_trial_basis_u1_0_1, INT64_C(2));
    bound_global_trial_basis_u1_0_1_stride_4 = nd_nstep_C(&global_trial_basis_u1_0_1, INT64_C(3));
    if (pyarray_check(global_trial_basis_u1_0_2_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_trial_basis_u1_0_2 = pyarray_to_ndarray(global_trial_basis_u1_0_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_trial_basis_u1_0_2");
        return NULL;
    }
    bound_global_trial_basis_u1_0_2 = nd_data(&global_trial_basis_u1_0_2);
    bound_global_trial_basis_u1_0_2_shape_1 = nd_ndim(&global_trial_basis_u1_0_2, INT64_C(0));
    bound_global_trial_basis_u1_0_2_shape_2 = nd_ndim(&global_trial_basis_u1_0_2, INT64_C(1));
    bound_global_trial_basis_u1_0_2_shape_3 = nd_ndim(&global_trial_basis_u1_0_2, INT64_C(2));
    bound_global_trial_basis_u1_0_2_shape_4 = nd_ndim(&global_trial_basis_u1_0_2, INT64_C(3));
    bound_global_trial_basis_u1_0_2_stride_1 = nd_nstep_C(&global_trial_basis_u1_0_2, INT64_C(0));
    bound_global_trial_basis_u1_0_2_stride_2 = nd_nstep_C(&global_trial_basis_u1_0_2, INT64_C(1));
    bound_global_trial_basis_u1_0_2_stride_3 = nd_nstep_C(&global_trial_basis_u1_0_2, INT64_C(2));
    bound_global_trial_basis_u1_0_2_stride_4 = nd_nstep_C(&global_trial_basis_u1_0_2, INT64_C(3));
    if (pyarray_check(global_trial_basis_u1_1_1_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_trial_basis_u1_1_1 = pyarray_to_ndarray(global_trial_basis_u1_1_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_trial_basis_u1_1_1");
        return NULL;
    }
    bound_global_trial_basis_u1_1_1 = nd_data(&global_trial_basis_u1_1_1);
    bound_global_trial_basis_u1_1_1_shape_1 = nd_ndim(&global_trial_basis_u1_1_1, INT64_C(0));
    bound_global_trial_basis_u1_1_1_shape_2 = nd_ndim(&global_trial_basis_u1_1_1, INT64_C(1));
    bound_global_trial_basis_u1_1_1_shape_3 = nd_ndim(&global_trial_basis_u1_1_1, INT64_C(2));
    bound_global_trial_basis_u1_1_1_shape_4 = nd_ndim(&global_trial_basis_u1_1_1, INT64_C(3));
    bound_global_trial_basis_u1_1_1_stride_1 = nd_nstep_C(&global_trial_basis_u1_1_1, INT64_C(0));
    bound_global_trial_basis_u1_1_1_stride_2 = nd_nstep_C(&global_trial_basis_u1_1_1, INT64_C(1));
    bound_global_trial_basis_u1_1_1_stride_3 = nd_nstep_C(&global_trial_basis_u1_1_1, INT64_C(2));
    bound_global_trial_basis_u1_1_1_stride_4 = nd_nstep_C(&global_trial_basis_u1_1_1, INT64_C(3));
    if (pyarray_check(global_trial_basis_u1_1_2_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_trial_basis_u1_1_2 = pyarray_to_ndarray(global_trial_basis_u1_1_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_trial_basis_u1_1_2");
        return NULL;
    }
    bound_global_trial_basis_u1_1_2 = nd_data(&global_trial_basis_u1_1_2);
    bound_global_trial_basis_u1_1_2_shape_1 = nd_ndim(&global_trial_basis_u1_1_2, INT64_C(0));
    bound_global_trial_basis_u1_1_2_shape_2 = nd_ndim(&global_trial_basis_u1_1_2, INT64_C(1));
    bound_global_trial_basis_u1_1_2_shape_3 = nd_ndim(&global_trial_basis_u1_1_2, INT64_C(2));
    bound_global_trial_basis_u1_1_2_shape_4 = nd_ndim(&global_trial_basis_u1_1_2, INT64_C(3));
    bound_global_trial_basis_u1_1_2_stride_1 = nd_nstep_C(&global_trial_basis_u1_1_2, INT64_C(0));
    bound_global_trial_basis_u1_1_2_stride_2 = nd_nstep_C(&global_trial_basis_u1_1_2, INT64_C(1));
    bound_global_trial_basis_u1_1_2_stride_3 = nd_nstep_C(&global_trial_basis_u1_1_2, INT64_C(2));
    bound_global_trial_basis_u1_1_2_stride_4 = nd_nstep_C(&global_trial_basis_u1_1_2, INT64_C(3));
    if (pyarray_check(global_span_v1_0_1_obj, NPY_LONG, INT64_C(1), NO_ORDER_CHECK))
    {
        global_span_v1_0_1 = pyarray_to_ndarray(global_span_v1_0_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument global_span_v1_0_1");
        return NULL;
    }
    bound_global_span_v1_0_1 = nd_data(&global_span_v1_0_1);
    bound_global_span_v1_0_1_shape_1 = nd_ndim(&global_span_v1_0_1, INT64_C(0));
    bound_global_span_v1_0_1_stride_1 = nd_nstep_F(&global_span_v1_0_1, INT64_C(0));
    if (pyarray_check(global_span_v1_0_2_obj, NPY_LONG, INT64_C(1), NO_ORDER_CHECK))
    {
        global_span_v1_0_2 = pyarray_to_ndarray(global_span_v1_0_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument global_span_v1_0_2");
        return NULL;
    }
    bound_global_span_v1_0_2 = nd_data(&global_span_v1_0_2);
    bound_global_span_v1_0_2_shape_1 = nd_ndim(&global_span_v1_0_2, INT64_C(0));
    bound_global_span_v1_0_2_stride_1 = nd_nstep_F(&global_span_v1_0_2, INT64_C(0));
    if (pyarray_check(global_span_v1_1_1_obj, NPY_LONG, INT64_C(1), NO_ORDER_CHECK))
    {
        global_span_v1_1_1 = pyarray_to_ndarray(global_span_v1_1_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument global_span_v1_1_1");
        return NULL;
    }
    bound_global_span_v1_1_1 = nd_data(&global_span_v1_1_1);
    bound_global_span_v1_1_1_shape_1 = nd_ndim(&global_span_v1_1_1, INT64_C(0));
    bound_global_span_v1_1_1_stride_1 = nd_nstep_F(&global_span_v1_1_1, INT64_C(0));
    if (pyarray_check(global_span_v1_1_2_obj, NPY_LONG, INT64_C(1), NO_ORDER_CHECK))
    {
        global_span_v1_1_2 = pyarray_to_ndarray(global_span_v1_1_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument global_span_v1_1_2");
        return NULL;
    }
    bound_global_span_v1_1_2 = nd_data(&global_span_v1_1_2);
    bound_global_span_v1_1_2_shape_1 = nd_ndim(&global_span_v1_1_2, INT64_C(0));
    bound_global_span_v1_1_2_stride_1 = nd_nstep_F(&global_span_v1_1_2, INT64_C(0));
    if (pyarray_check(global_x1_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_x1 = pyarray_to_ndarray(global_x1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_x1");
        return NULL;
    }
    bound_global_x1 = nd_data(&global_x1);
    bound_global_x1_shape_1 = nd_ndim(&global_x1, INT64_C(0));
    bound_global_x1_shape_2 = nd_ndim(&global_x1, INT64_C(1));
    bound_global_x1_stride_1 = nd_nstep_C(&global_x1, INT64_C(0));
    bound_global_x1_stride_2 = nd_nstep_C(&global_x1, INT64_C(1));
    if (pyarray_check(global_x2_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        global_x2 = pyarray_to_ndarray(global_x2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument global_x2");
        return NULL;
    }
    bound_global_x2 = nd_data(&global_x2);
    bound_global_x2_shape_1 = nd_ndim(&global_x2, INT64_C(0));
    bound_global_x2_shape_2 = nd_ndim(&global_x2, INT64_C(1));
    bound_global_x2_stride_1 = nd_nstep_C(&global_x2, INT64_C(0));
    bound_global_x2_stride_2 = nd_nstep_C(&global_x2, INT64_C(1));
    if (PyIs_Int64(test_v1_0_p1_obj))
    {
        test_v1_0_p1 = PyInt64_to_Int64(test_v1_0_p1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument test_v1_0_p1");
        return NULL;
    }
    if (PyIs_Int64(test_v1_0_p2_obj))
    {
        test_v1_0_p2 = PyInt64_to_Int64(test_v1_0_p2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument test_v1_0_p2");
        return NULL;
    }
    if (PyIs_Int64(test_v1_1_p1_obj))
    {
        test_v1_1_p1 = PyInt64_to_Int64(test_v1_1_p1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument test_v1_1_p1");
        return NULL;
    }
    if (PyIs_Int64(test_v1_1_p2_obj))
    {
        test_v1_1_p2 = PyInt64_to_Int64(test_v1_1_p2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument test_v1_1_p2");
        return NULL;
    }
    if (PyIs_Int64(trial_u1_0_p1_obj))
    {
        trial_u1_0_p1 = PyInt64_to_Int64(trial_u1_0_p1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument trial_u1_0_p1");
        return NULL;
    }
    if (PyIs_Int64(trial_u1_0_p2_obj))
    {
        trial_u1_0_p2 = PyInt64_to_Int64(trial_u1_0_p2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument trial_u1_0_p2");
        return NULL;
    }
    if (PyIs_Int64(trial_u1_1_p1_obj))
    {
        trial_u1_1_p1 = PyInt64_to_Int64(trial_u1_1_p1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument trial_u1_1_p1");
        return NULL;
    }
    if (PyIs_Int64(trial_u1_1_p2_obj))
    {
        trial_u1_1_p2 = PyInt64_to_Int64(trial_u1_1_p2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument trial_u1_1_p2");
        return NULL;
    }
    if (PyIs_Int64(n_element_1_obj))
    {
        n_element_1 = PyInt64_to_Int64(n_element_1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n_element_1");
        return NULL;
    }
    if (PyIs_Int64(n_element_2_obj))
    {
        n_element_2 = PyInt64_to_Int64(n_element_2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument n_element_2");
        return NULL;
    }
    if (PyIs_Int64(k1_obj))
    {
        k1 = PyInt64_to_Int64(k1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument k1");
        return NULL;
    }
    if (PyIs_Int64(k2_obj))
    {
        k2 = PyInt64_to_Int64(k2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument k2");
        return NULL;
    }
    if (PyIs_Int64(pad1_obj))
    {
        pad1 = PyInt64_to_Int64(pad1_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument pad1");
        return NULL;
    }
    if (PyIs_Int64(pad2_obj))
    {
        pad2 = PyInt64_to_Int64(pad2_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument pad2");
        return NULL;
    }
    if (pyarray_check(g_mat_u1_0_v1_0_e0el77nq_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        g_mat_u1_0_v1_0_e0el77nq = pyarray_to_ndarray(g_mat_u1_0_v1_0_e0el77nq_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument g_mat_u1_0_v1_0_e0el77nq");
        return NULL;
    }
    bound_g_mat_u1_0_v1_0_e0el77nq = nd_data(&g_mat_u1_0_v1_0_e0el77nq);
    bound_g_mat_u1_0_v1_0_e0el77nq_shape_1 = nd_ndim(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(0));
    bound_g_mat_u1_0_v1_0_e0el77nq_shape_2 = nd_ndim(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(1));
    bound_g_mat_u1_0_v1_0_e0el77nq_shape_3 = nd_ndim(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(2));
    bound_g_mat_u1_0_v1_0_e0el77nq_shape_4 = nd_ndim(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(3));
    bound_g_mat_u1_0_v1_0_e0el77nq_stride_1 = nd_nstep_C(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(0));
    bound_g_mat_u1_0_v1_0_e0el77nq_stride_2 = nd_nstep_C(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(1));
    bound_g_mat_u1_0_v1_0_e0el77nq_stride_3 = nd_nstep_C(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(2));
    bound_g_mat_u1_0_v1_0_e0el77nq_stride_4 = nd_nstep_C(&g_mat_u1_0_v1_0_e0el77nq, INT64_C(3));
    if (pyarray_check(g_mat_u1_1_v1_0_e0el77nq_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        g_mat_u1_1_v1_0_e0el77nq = pyarray_to_ndarray(g_mat_u1_1_v1_0_e0el77nq_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument g_mat_u1_1_v1_0_e0el77nq");
        return NULL;
    }
    bound_g_mat_u1_1_v1_0_e0el77nq = nd_data(&g_mat_u1_1_v1_0_e0el77nq);
    bound_g_mat_u1_1_v1_0_e0el77nq_shape_1 = nd_ndim(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(0));
    bound_g_mat_u1_1_v1_0_e0el77nq_shape_2 = nd_ndim(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(1));
    bound_g_mat_u1_1_v1_0_e0el77nq_shape_3 = nd_ndim(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(2));
    bound_g_mat_u1_1_v1_0_e0el77nq_shape_4 = nd_ndim(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(3));
    bound_g_mat_u1_1_v1_0_e0el77nq_stride_1 = nd_nstep_C(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(0));
    bound_g_mat_u1_1_v1_0_e0el77nq_stride_2 = nd_nstep_C(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(1));
    bound_g_mat_u1_1_v1_0_e0el77nq_stride_3 = nd_nstep_C(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(2));
    bound_g_mat_u1_1_v1_0_e0el77nq_stride_4 = nd_nstep_C(&g_mat_u1_1_v1_0_e0el77nq, INT64_C(3));
    if (pyarray_check(g_mat_u1_0_v1_1_e0el77nq_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        g_mat_u1_0_v1_1_e0el77nq = pyarray_to_ndarray(g_mat_u1_0_v1_1_e0el77nq_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument g_mat_u1_0_v1_1_e0el77nq");
        return NULL;
    }
    bound_g_mat_u1_0_v1_1_e0el77nq = nd_data(&g_mat_u1_0_v1_1_e0el77nq);
    bound_g_mat_u1_0_v1_1_e0el77nq_shape_1 = nd_ndim(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(0));
    bound_g_mat_u1_0_v1_1_e0el77nq_shape_2 = nd_ndim(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(1));
    bound_g_mat_u1_0_v1_1_e0el77nq_shape_3 = nd_ndim(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(2));
    bound_g_mat_u1_0_v1_1_e0el77nq_shape_4 = nd_ndim(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(3));
    bound_g_mat_u1_0_v1_1_e0el77nq_stride_1 = nd_nstep_C(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(0));
    bound_g_mat_u1_0_v1_1_e0el77nq_stride_2 = nd_nstep_C(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(1));
    bound_g_mat_u1_0_v1_1_e0el77nq_stride_3 = nd_nstep_C(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(2));
    bound_g_mat_u1_0_v1_1_e0el77nq_stride_4 = nd_nstep_C(&g_mat_u1_0_v1_1_e0el77nq, INT64_C(3));
    if (pyarray_check(g_mat_u1_1_v1_1_e0el77nq_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        g_mat_u1_1_v1_1_e0el77nq = pyarray_to_ndarray(g_mat_u1_1_v1_1_e0el77nq_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument g_mat_u1_1_v1_1_e0el77nq");
        return NULL;
    }
    bound_g_mat_u1_1_v1_1_e0el77nq = nd_data(&g_mat_u1_1_v1_1_e0el77nq);
    bound_g_mat_u1_1_v1_1_e0el77nq_shape_1 = nd_ndim(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(0));
    bound_g_mat_u1_1_v1_1_e0el77nq_shape_2 = nd_ndim(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(1));
    bound_g_mat_u1_1_v1_1_e0el77nq_shape_3 = nd_ndim(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(2));
    bound_g_mat_u1_1_v1_1_e0el77nq_shape_4 = nd_ndim(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(3));
    bound_g_mat_u1_1_v1_1_e0el77nq_stride_1 = nd_nstep_C(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(0));
    bound_g_mat_u1_1_v1_1_e0el77nq_stride_2 = nd_nstep_C(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(1));
    bound_g_mat_u1_1_v1_1_e0el77nq_stride_3 = nd_nstep_C(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(2));
    bound_g_mat_u1_1_v1_1_e0el77nq_stride_4 = nd_nstep_C(&g_mat_u1_1_v1_1_e0el77nq, INT64_C(3));
    bind_c_assemble_matrix_e0el77nq(bound_global_test_basis_v1_0_1, bound_global_test_basis_v1_0_1_shape_1, bound_global_test_basis_v1_0_1_shape_2, bound_global_test_basis_v1_0_1_shape_3, bound_global_test_basis_v1_0_1_shape_4, bound_global_test_basis_v1_0_1_stride_1, bound_global_test_basis_v1_0_1_stride_2, bound_global_test_basis_v1_0_1_stride_3, bound_global_test_basis_v1_0_1_stride_4, bound_global_test_basis_v1_0_2, bound_global_test_basis_v1_0_2_shape_1, bound_global_test_basis_v1_0_2_shape_2, bound_global_test_basis_v1_0_2_shape_3, bound_global_test_basis_v1_0_2_shape_4, bound_global_test_basis_v1_0_2_stride_1, bound_global_test_basis_v1_0_2_stride_2, bound_global_test_basis_v1_0_2_stride_3, bound_global_test_basis_v1_0_2_stride_4, bound_global_test_basis_v1_1_1, bound_global_test_basis_v1_1_1_shape_1, bound_global_test_basis_v1_1_1_shape_2, bound_global_test_basis_v1_1_1_shape_3, bound_global_test_basis_v1_1_1_shape_4, bound_global_test_basis_v1_1_1_stride_1, bound_global_test_basis_v1_1_1_stride_2, bound_global_test_basis_v1_1_1_stride_3, bound_global_test_basis_v1_1_1_stride_4, bound_global_test_basis_v1_1_2, bound_global_test_basis_v1_1_2_shape_1, bound_global_test_basis_v1_1_2_shape_2, bound_global_test_basis_v1_1_2_shape_3, bound_global_test_basis_v1_1_2_shape_4, bound_global_test_basis_v1_1_2_stride_1, bound_global_test_basis_v1_1_2_stride_2, bound_global_test_basis_v1_1_2_stride_3, bound_global_test_basis_v1_1_2_stride_4, bound_global_trial_basis_u1_0_1, bound_global_trial_basis_u1_0_1_shape_1, bound_global_trial_basis_u1_0_1_shape_2, bound_global_trial_basis_u1_0_1_shape_3, bound_global_trial_basis_u1_0_1_shape_4, bound_global_trial_basis_u1_0_1_stride_1, bound_global_trial_basis_u1_0_1_stride_2, bound_global_trial_basis_u1_0_1_stride_3, bound_global_trial_basis_u1_0_1_stride_4, bound_global_trial_basis_u1_0_2, bound_global_trial_basis_u1_0_2_shape_1, bound_global_trial_basis_u1_0_2_shape_2, bound_global_trial_basis_u1_0_2_shape_3, bound_global_trial_basis_u1_0_2_shape_4, bound_global_trial_basis_u1_0_2_stride_1, bound_global_trial_basis_u1_0_2_stride_2, bound_global_trial_basis_u1_0_2_stride_3, bound_global_trial_basis_u1_0_2_stride_4, bound_global_trial_basis_u1_1_1, bound_global_trial_basis_u1_1_1_shape_1, bound_global_trial_basis_u1_1_1_shape_2, bound_global_trial_basis_u1_1_1_shape_3, bound_global_trial_basis_u1_1_1_shape_4, bound_global_trial_basis_u1_1_1_stride_1, bound_global_trial_basis_u1_1_1_stride_2, bound_global_trial_basis_u1_1_1_stride_3, bound_global_trial_basis_u1_1_1_stride_4, bound_global_trial_basis_u1_1_2, bound_global_trial_basis_u1_1_2_shape_1, bound_global_trial_basis_u1_1_2_shape_2, bound_global_trial_basis_u1_1_2_shape_3, bound_global_trial_basis_u1_1_2_shape_4, bound_global_trial_basis_u1_1_2_stride_1, bound_global_trial_basis_u1_1_2_stride_2, bound_global_trial_basis_u1_1_2_stride_3, bound_global_trial_basis_u1_1_2_stride_4, bound_global_span_v1_0_1, bound_global_span_v1_0_1_shape_1, bound_global_span_v1_0_1_stride_1, bound_global_span_v1_0_2, bound_global_span_v1_0_2_shape_1, bound_global_span_v1_0_2_stride_1, bound_global_span_v1_1_1, bound_global_span_v1_1_1_shape_1, bound_global_span_v1_1_1_stride_1, bound_global_span_v1_1_2, bound_global_span_v1_1_2_shape_1, bound_global_span_v1_1_2_stride_1, bound_global_x1, bound_global_x1_shape_1, bound_global_x1_shape_2, bound_global_x1_stride_1, bound_global_x1_stride_2, bound_global_x2, bound_global_x2_shape_1, bound_global_x2_shape_2, bound_global_x2_stride_1, bound_global_x2_stride_2, test_v1_0_p1, test_v1_0_p2, test_v1_1_p1, test_v1_1_p2, trial_u1_0_p1, trial_u1_0_p2, trial_u1_1_p1, trial_u1_1_p2, n_element_1, n_element_2, k1, k2, pad1, pad2, bound_g_mat_u1_0_v1_0_e0el77nq, bound_g_mat_u1_0_v1_0_e0el77nq_shape_1, bound_g_mat_u1_0_v1_0_e0el77nq_shape_2, bound_g_mat_u1_0_v1_0_e0el77nq_shape_3, bound_g_mat_u1_0_v1_0_e0el77nq_shape_4, bound_g_mat_u1_0_v1_0_e0el77nq_stride_1, bound_g_mat_u1_0_v1_0_e0el77nq_stride_2, bound_g_mat_u1_0_v1_0_e0el77nq_stride_3, bound_g_mat_u1_0_v1_0_e0el77nq_stride_4, bound_g_mat_u1_1_v1_0_e0el77nq, bound_g_mat_u1_1_v1_0_e0el77nq_shape_1, bound_g_mat_u1_1_v1_0_e0el77nq_shape_2, bound_g_mat_u1_1_v1_0_e0el77nq_shape_3, bound_g_mat_u1_1_v1_0_e0el77nq_shape_4, bound_g_mat_u1_1_v1_0_e0el77nq_stride_1, bound_g_mat_u1_1_v1_0_e0el77nq_stride_2, bound_g_mat_u1_1_v1_0_e0el77nq_stride_3, bound_g_mat_u1_1_v1_0_e0el77nq_stride_4, bound_g_mat_u1_0_v1_1_e0el77nq, bound_g_mat_u1_0_v1_1_e0el77nq_shape_1, bound_g_mat_u1_0_v1_1_e0el77nq_shape_2, bound_g_mat_u1_0_v1_1_e0el77nq_shape_3, bound_g_mat_u1_0_v1_1_e0el77nq_shape_4, bound_g_mat_u1_0_v1_1_e0el77nq_stride_1, bound_g_mat_u1_0_v1_1_e0el77nq_stride_2, bound_g_mat_u1_0_v1_1_e0el77nq_stride_3, bound_g_mat_u1_0_v1_1_e0el77nq_stride_4, bound_g_mat_u1_1_v1_1_e0el77nq, bound_g_mat_u1_1_v1_1_e0el77nq_shape_1, bound_g_mat_u1_1_v1_1_e0el77nq_shape_2, bound_g_mat_u1_1_v1_1_e0el77nq_shape_3, bound_g_mat_u1_1_v1_1_e0el77nq_shape_4, bound_g_mat_u1_1_v1_1_e0el77nq_stride_1, bound_g_mat_u1_1_v1_1_e0el77nq_stride_2, bound_g_mat_u1_1_v1_1_e0el77nq_stride_3, bound_g_mat_u1_1_v1_1_e0el77nq_stride_4);
    free_pointer(&global_test_basis_v1_0_1);
    free_pointer(&global_test_basis_v1_0_2);
    free_pointer(&global_test_basis_v1_1_1);
    free_pointer(&global_test_basis_v1_1_2);
    free_pointer(&global_trial_basis_u1_0_1);
    free_pointer(&global_trial_basis_u1_0_2);
    free_pointer(&global_trial_basis_u1_1_1);
    free_pointer(&global_trial_basis_u1_1_2);
    free_pointer(&global_span_v1_0_1);
    free_pointer(&global_span_v1_0_2);
    free_pointer(&global_span_v1_1_1);
    free_pointer(&global_span_v1_1_2);
    free_pointer(&global_x1);
    free_pointer(&global_x2);
    free_pointer(&g_mat_u1_0_v1_0_e0el77nq);
    free_pointer(&g_mat_u1_1_v1_0_e0el77nq);
    free_pointer(&g_mat_u1_0_v1_1_e0el77nq);
    free_pointer(&g_mat_u1_1_v1_1_e0el77nq);
    Py_INCREF(Py_None);
    return Py_None;
}
/*........................................*/

/*........................................*/

static PyMethodDef dependencies_e0el77nq_mx2lvwh9fed2_methods[] = {
    {
        "assemble_matrix_e0el77nq",
        (PyCFunction)bind_c_assemble_matrix_e0el77nq_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef dependencies_e0el77nq_mx2lvwh9fed2_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "dependencies_e0el77nq_mx2lvwh9fed2",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    dependencies_e0el77nq_mx2lvwh9fed2_methods,
};

/*........................................*/

PyMODINIT_FUNC PyInit_dependencies_e0el77nq_mx2lvwh9fed2(void)
{
    PyObject* mod;
    static void* Pydependencies_e0el77nq_mx2lvwh9fed2_API[0];
    PyObject* c_api_object_0001;
    mod = PyModule_Create(&dependencies_e0el77nq_mx2lvwh9fed2_module);
    if (mod == NULL)
    {
        return NULL;
    }
    c_api_object_0001 = PyCapsule_New((void *)Pydependencies_e0el77nq_mx2lvwh9fed2_API, "dependencies_e0el77nq_mx2lvwh9fed2._C_API", NULL);
    if (PyModule_AddObject(mod, "_C_API", c_api_object_0001) < INT64_C(0))
    {
        Py_DECREF(mod);
        return NULL;
    }
    Py_INCREF(c_api_object_0001);
    import_array();
    return mod;
}
