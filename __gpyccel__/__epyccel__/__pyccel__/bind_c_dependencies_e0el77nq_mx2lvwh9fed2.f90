module bind_c_dependencies_e0el77nq_mx2lvwh9fed2

  use dependencies_e0el77nq_mx2lvwh9fed2

  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T , C_F_Pointer , c_ptr
  implicit none

  contains

  !........................................
  subroutine bind_c_assemble_matrix_e0el77nq( &
        bound_global_test_basis_v1_0_1, &
        bound_global_test_basis_v1_0_1_shape_1, &
        bound_global_test_basis_v1_0_1_shape_2, &
        bound_global_test_basis_v1_0_1_shape_3, &
        bound_global_test_basis_v1_0_1_shape_4, &
        bound_global_test_basis_v1_0_1_stride_1, &
        bound_global_test_basis_v1_0_1_stride_2, &
        bound_global_test_basis_v1_0_1_stride_3, &
        bound_global_test_basis_v1_0_1_stride_4, &
        bound_global_test_basis_v1_0_2, &
        bound_global_test_basis_v1_0_2_shape_1, &
        bound_global_test_basis_v1_0_2_shape_2, &
        bound_global_test_basis_v1_0_2_shape_3, &
        bound_global_test_basis_v1_0_2_shape_4, &
        bound_global_test_basis_v1_0_2_stride_1, &
        bound_global_test_basis_v1_0_2_stride_2, &
        bound_global_test_basis_v1_0_2_stride_3, &
        bound_global_test_basis_v1_0_2_stride_4, &
        bound_global_test_basis_v1_1_1, &
        bound_global_test_basis_v1_1_1_shape_1, &
        bound_global_test_basis_v1_1_1_shape_2, &
        bound_global_test_basis_v1_1_1_shape_3, &
        bound_global_test_basis_v1_1_1_shape_4, &
        bound_global_test_basis_v1_1_1_stride_1, &
        bound_global_test_basis_v1_1_1_stride_2, &
        bound_global_test_basis_v1_1_1_stride_3, &
        bound_global_test_basis_v1_1_1_stride_4, &
        bound_global_test_basis_v1_1_2, &
        bound_global_test_basis_v1_1_2_shape_1, &
        bound_global_test_basis_v1_1_2_shape_2, &
        bound_global_test_basis_v1_1_2_shape_3, &
        bound_global_test_basis_v1_1_2_shape_4, &
        bound_global_test_basis_v1_1_2_stride_1, &
        bound_global_test_basis_v1_1_2_stride_2, &
        bound_global_test_basis_v1_1_2_stride_3, &
        bound_global_test_basis_v1_1_2_stride_4, &
        bound_global_trial_basis_u1_0_1, &
        bound_global_trial_basis_u1_0_1_shape_1, &
        bound_global_trial_basis_u1_0_1_shape_2, &
        bound_global_trial_basis_u1_0_1_shape_3, &
        bound_global_trial_basis_u1_0_1_shape_4, &
        bound_global_trial_basis_u1_0_1_stride_1, &
        bound_global_trial_basis_u1_0_1_stride_2, &
        bound_global_trial_basis_u1_0_1_stride_3, &
        bound_global_trial_basis_u1_0_1_stride_4, &
        bound_global_trial_basis_u1_0_2, &
        bound_global_trial_basis_u1_0_2_shape_1, &
        bound_global_trial_basis_u1_0_2_shape_2, &
        bound_global_trial_basis_u1_0_2_shape_3, &
        bound_global_trial_basis_u1_0_2_shape_4, &
        bound_global_trial_basis_u1_0_2_stride_1, &
        bound_global_trial_basis_u1_0_2_stride_2, &
        bound_global_trial_basis_u1_0_2_stride_3, &
        bound_global_trial_basis_u1_0_2_stride_4, &
        bound_global_trial_basis_u1_1_1, &
        bound_global_trial_basis_u1_1_1_shape_1, &
        bound_global_trial_basis_u1_1_1_shape_2, &
        bound_global_trial_basis_u1_1_1_shape_3, &
        bound_global_trial_basis_u1_1_1_shape_4, &
        bound_global_trial_basis_u1_1_1_stride_1, &
        bound_global_trial_basis_u1_1_1_stride_2, &
        bound_global_trial_basis_u1_1_1_stride_3, &
        bound_global_trial_basis_u1_1_1_stride_4, &
        bound_global_trial_basis_u1_1_2, &
        bound_global_trial_basis_u1_1_2_shape_1, &
        bound_global_trial_basis_u1_1_2_shape_2, &
        bound_global_trial_basis_u1_1_2_shape_3, &
        bound_global_trial_basis_u1_1_2_shape_4, &
        bound_global_trial_basis_u1_1_2_stride_1, &
        bound_global_trial_basis_u1_1_2_stride_2, &
        bound_global_trial_basis_u1_1_2_stride_3, &
        bound_global_trial_basis_u1_1_2_stride_4, &
        bound_global_span_v1_0_1, bound_global_span_v1_0_1_shape_1, &
        bound_global_span_v1_0_1_stride_1, bound_global_span_v1_0_2, &
        bound_global_span_v1_0_2_shape_1, &
        bound_global_span_v1_0_2_stride_1, bound_global_span_v1_1_1, &
        bound_global_span_v1_1_1_shape_1, &
        bound_global_span_v1_1_1_stride_1, bound_global_span_v1_1_2, &
        bound_global_span_v1_1_2_shape_1, &
        bound_global_span_v1_1_2_stride_1, bound_global_x1, &
        bound_global_x1_shape_1, bound_global_x1_shape_2, &
        bound_global_x1_stride_1, bound_global_x1_stride_2, &
        bound_global_x2, bound_global_x2_shape_1, &
        bound_global_x2_shape_2, bound_global_x2_stride_1, &
        bound_global_x2_stride_2, test_v1_0_p1, test_v1_0_p2, &
        test_v1_1_p1, test_v1_1_p2, trial_u1_0_p1, trial_u1_0_p2, &
        trial_u1_1_p1, trial_u1_1_p2, n_element_1, n_element_2, k1, k2, &
        pad1, pad2, bound_g_mat_u1_0_v1_0_e0el77nq, &
        bound_g_mat_u1_0_v1_0_e0el77nq_shape_1, &
        bound_g_mat_u1_0_v1_0_e0el77nq_shape_2, &
        bound_g_mat_u1_0_v1_0_e0el77nq_shape_3, &
        bound_g_mat_u1_0_v1_0_e0el77nq_shape_4, &
        bound_g_mat_u1_0_v1_0_e0el77nq_stride_1, &
        bound_g_mat_u1_0_v1_0_e0el77nq_stride_2, &
        bound_g_mat_u1_0_v1_0_e0el77nq_stride_3, &
        bound_g_mat_u1_0_v1_0_e0el77nq_stride_4, &
        bound_g_mat_u1_1_v1_0_e0el77nq, &
        bound_g_mat_u1_1_v1_0_e0el77nq_shape_1, &
        bound_g_mat_u1_1_v1_0_e0el77nq_shape_2, &
        bound_g_mat_u1_1_v1_0_e0el77nq_shape_3, &
        bound_g_mat_u1_1_v1_0_e0el77nq_shape_4, &
        bound_g_mat_u1_1_v1_0_e0el77nq_stride_1, &
        bound_g_mat_u1_1_v1_0_e0el77nq_stride_2, &
        bound_g_mat_u1_1_v1_0_e0el77nq_stride_3, &
        bound_g_mat_u1_1_v1_0_e0el77nq_stride_4, &
        bound_g_mat_u1_0_v1_1_e0el77nq, &
        bound_g_mat_u1_0_v1_1_e0el77nq_shape_1, &
        bound_g_mat_u1_0_v1_1_e0el77nq_shape_2, &
        bound_g_mat_u1_0_v1_1_e0el77nq_shape_3, &
        bound_g_mat_u1_0_v1_1_e0el77nq_shape_4, &
        bound_g_mat_u1_0_v1_1_e0el77nq_stride_1, &
        bound_g_mat_u1_0_v1_1_e0el77nq_stride_2, &
        bound_g_mat_u1_0_v1_1_e0el77nq_stride_3, &
        bound_g_mat_u1_0_v1_1_e0el77nq_stride_4, &
        bound_g_mat_u1_1_v1_1_e0el77nq, &
        bound_g_mat_u1_1_v1_1_e0el77nq_shape_1, &
        bound_g_mat_u1_1_v1_1_e0el77nq_shape_2, &
        bound_g_mat_u1_1_v1_1_e0el77nq_shape_3, &
        bound_g_mat_u1_1_v1_1_e0el77nq_shape_4, &
        bound_g_mat_u1_1_v1_1_e0el77nq_stride_1, &
        bound_g_mat_u1_1_v1_1_e0el77nq_stride_2, &
        bound_g_mat_u1_1_v1_1_e0el77nq_stride_3, &
        bound_g_mat_u1_1_v1_1_e0el77nq_stride_4) bind(c)

    implicit none

    type(c_ptr), value :: bound_global_test_basis_v1_0_1
    integer(i64), value :: bound_global_test_basis_v1_0_1_shape_1
    integer(i64), value :: bound_global_test_basis_v1_0_1_shape_2
    integer(i64), value :: bound_global_test_basis_v1_0_1_shape_3
    integer(i64), value :: bound_global_test_basis_v1_0_1_shape_4
    integer(i64), value :: bound_global_test_basis_v1_0_1_stride_1
    integer(i64), value :: bound_global_test_basis_v1_0_1_stride_2
    integer(i64), value :: bound_global_test_basis_v1_0_1_stride_3
    integer(i64), value :: bound_global_test_basis_v1_0_1_stride_4
    type(c_ptr), value :: bound_global_test_basis_v1_0_2
    integer(i64), value :: bound_global_test_basis_v1_0_2_shape_1
    integer(i64), value :: bound_global_test_basis_v1_0_2_shape_2
    integer(i64), value :: bound_global_test_basis_v1_0_2_shape_3
    integer(i64), value :: bound_global_test_basis_v1_0_2_shape_4
    integer(i64), value :: bound_global_test_basis_v1_0_2_stride_1
    integer(i64), value :: bound_global_test_basis_v1_0_2_stride_2
    integer(i64), value :: bound_global_test_basis_v1_0_2_stride_3
    integer(i64), value :: bound_global_test_basis_v1_0_2_stride_4
    type(c_ptr), value :: bound_global_test_basis_v1_1_1
    integer(i64), value :: bound_global_test_basis_v1_1_1_shape_1
    integer(i64), value :: bound_global_test_basis_v1_1_1_shape_2
    integer(i64), value :: bound_global_test_basis_v1_1_1_shape_3
    integer(i64), value :: bound_global_test_basis_v1_1_1_shape_4
    integer(i64), value :: bound_global_test_basis_v1_1_1_stride_1
    integer(i64), value :: bound_global_test_basis_v1_1_1_stride_2
    integer(i64), value :: bound_global_test_basis_v1_1_1_stride_3
    integer(i64), value :: bound_global_test_basis_v1_1_1_stride_4
    type(c_ptr), value :: bound_global_test_basis_v1_1_2
    integer(i64), value :: bound_global_test_basis_v1_1_2_shape_1
    integer(i64), value :: bound_global_test_basis_v1_1_2_shape_2
    integer(i64), value :: bound_global_test_basis_v1_1_2_shape_3
    integer(i64), value :: bound_global_test_basis_v1_1_2_shape_4
    integer(i64), value :: bound_global_test_basis_v1_1_2_stride_1
    integer(i64), value :: bound_global_test_basis_v1_1_2_stride_2
    integer(i64), value :: bound_global_test_basis_v1_1_2_stride_3
    integer(i64), value :: bound_global_test_basis_v1_1_2_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u1_0_1
    integer(i64), value :: bound_global_trial_basis_u1_0_1_shape_1
    integer(i64), value :: bound_global_trial_basis_u1_0_1_shape_2
    integer(i64), value :: bound_global_trial_basis_u1_0_1_shape_3
    integer(i64), value :: bound_global_trial_basis_u1_0_1_shape_4
    integer(i64), value :: bound_global_trial_basis_u1_0_1_stride_1
    integer(i64), value :: bound_global_trial_basis_u1_0_1_stride_2
    integer(i64), value :: bound_global_trial_basis_u1_0_1_stride_3
    integer(i64), value :: bound_global_trial_basis_u1_0_1_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u1_0_2
    integer(i64), value :: bound_global_trial_basis_u1_0_2_shape_1
    integer(i64), value :: bound_global_trial_basis_u1_0_2_shape_2
    integer(i64), value :: bound_global_trial_basis_u1_0_2_shape_3
    integer(i64), value :: bound_global_trial_basis_u1_0_2_shape_4
    integer(i64), value :: bound_global_trial_basis_u1_0_2_stride_1
    integer(i64), value :: bound_global_trial_basis_u1_0_2_stride_2
    integer(i64), value :: bound_global_trial_basis_u1_0_2_stride_3
    integer(i64), value :: bound_global_trial_basis_u1_0_2_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u1_1_1
    integer(i64), value :: bound_global_trial_basis_u1_1_1_shape_1
    integer(i64), value :: bound_global_trial_basis_u1_1_1_shape_2
    integer(i64), value :: bound_global_trial_basis_u1_1_1_shape_3
    integer(i64), value :: bound_global_trial_basis_u1_1_1_shape_4
    integer(i64), value :: bound_global_trial_basis_u1_1_1_stride_1
    integer(i64), value :: bound_global_trial_basis_u1_1_1_stride_2
    integer(i64), value :: bound_global_trial_basis_u1_1_1_stride_3
    integer(i64), value :: bound_global_trial_basis_u1_1_1_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u1_1_2
    integer(i64), value :: bound_global_trial_basis_u1_1_2_shape_1
    integer(i64), value :: bound_global_trial_basis_u1_1_2_shape_2
    integer(i64), value :: bound_global_trial_basis_u1_1_2_shape_3
    integer(i64), value :: bound_global_trial_basis_u1_1_2_shape_4
    integer(i64), value :: bound_global_trial_basis_u1_1_2_stride_1
    integer(i64), value :: bound_global_trial_basis_u1_1_2_stride_2
    integer(i64), value :: bound_global_trial_basis_u1_1_2_stride_3
    integer(i64), value :: bound_global_trial_basis_u1_1_2_stride_4
    type(c_ptr), value :: bound_global_span_v1_0_1
    integer(i64), value :: bound_global_span_v1_0_1_shape_1
    integer(i64), value :: bound_global_span_v1_0_1_stride_1
    type(c_ptr), value :: bound_global_span_v1_0_2
    integer(i64), value :: bound_global_span_v1_0_2_shape_1
    integer(i64), value :: bound_global_span_v1_0_2_stride_1
    type(c_ptr), value :: bound_global_span_v1_1_1
    integer(i64), value :: bound_global_span_v1_1_1_shape_1
    integer(i64), value :: bound_global_span_v1_1_1_stride_1
    type(c_ptr), value :: bound_global_span_v1_1_2
    integer(i64), value :: bound_global_span_v1_1_2_shape_1
    integer(i64), value :: bound_global_span_v1_1_2_stride_1
    type(c_ptr), value :: bound_global_x1
    integer(i64), value :: bound_global_x1_shape_1
    integer(i64), value :: bound_global_x1_shape_2
    integer(i64), value :: bound_global_x1_stride_1
    integer(i64), value :: bound_global_x1_stride_2
    type(c_ptr), value :: bound_global_x2
    integer(i64), value :: bound_global_x2_shape_1
    integer(i64), value :: bound_global_x2_shape_2
    integer(i64), value :: bound_global_x2_stride_1
    integer(i64), value :: bound_global_x2_stride_2
    integer(i64), value :: test_v1_0_p1
    integer(i64), value :: test_v1_0_p2
    integer(i64), value :: test_v1_1_p1
    integer(i64), value :: test_v1_1_p2
    integer(i64), value :: trial_u1_0_p1
    integer(i64), value :: trial_u1_0_p2
    integer(i64), value :: trial_u1_1_p1
    integer(i64), value :: trial_u1_1_p2
    integer(i64), value :: n_element_1
    integer(i64), value :: n_element_2
    integer(i64), value :: k1
    integer(i64), value :: k2
    integer(i64), value :: pad1
    integer(i64), value :: pad2
    type(c_ptr), value :: bound_g_mat_u1_0_v1_0_e0el77nq
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_shape_1
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_shape_2
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_shape_3
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_shape_4
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_stride_1
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_stride_2
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_stride_3
    integer(i64), value :: bound_g_mat_u1_0_v1_0_e0el77nq_stride_4
    type(c_ptr), value :: bound_g_mat_u1_1_v1_0_e0el77nq
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_shape_1
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_shape_2
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_shape_3
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_shape_4
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_stride_1
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_stride_2
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_stride_3
    integer(i64), value :: bound_g_mat_u1_1_v1_0_e0el77nq_stride_4
    type(c_ptr), value :: bound_g_mat_u1_0_v1_1_e0el77nq
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_shape_1
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_shape_2
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_shape_3
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_shape_4
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_stride_1
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_stride_2
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_stride_3
    integer(i64), value :: bound_g_mat_u1_0_v1_1_e0el77nq_stride_4
    type(c_ptr), value :: bound_g_mat_u1_1_v1_1_e0el77nq
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_shape_1
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_shape_2
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_shape_3
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_shape_4
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_stride_1
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_stride_2
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_stride_3
    integer(i64), value :: bound_g_mat_u1_1_v1_1_e0el77nq_stride_4
    real(f64), pointer :: global_test_basis_v1_0_1(:,:,:,:)
    real(f64), pointer :: global_test_basis_v1_0_2(:,:,:,:)
    real(f64), pointer :: global_test_basis_v1_1_1(:,:,:,:)
    real(f64), pointer :: global_test_basis_v1_1_2(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u1_0_1(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u1_0_2(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u1_1_1(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u1_1_2(:,:,:,:)
    integer(i64), pointer :: global_span_v1_0_1(:)
    integer(i64), pointer :: global_span_v1_0_2(:)
    integer(i64), pointer :: global_span_v1_1_1(:)
    integer(i64), pointer :: global_span_v1_1_2(:)
    real(f64), pointer :: global_x1(:,:)
    real(f64), pointer :: global_x2(:,:)
    real(f64), pointer :: g_mat_u1_0_v1_0_e0el77nq(:,:,:,:)
    real(f64), pointer :: g_mat_u1_1_v1_0_e0el77nq(:,:,:,:)
    real(f64), pointer :: g_mat_u1_0_v1_1_e0el77nq(:,:,:,:)
    real(f64), pointer :: g_mat_u1_1_v1_1_e0el77nq(:,:,:,:)

    call C_F_Pointer(bound_global_test_basis_v1_0_1, &
          global_test_basis_v1_0_1, [ &
          bound_global_test_basis_v1_0_1_shape_4 * &
          bound_global_test_basis_v1_0_1_stride_4, &
          bound_global_test_basis_v1_0_1_shape_3 * &
          bound_global_test_basis_v1_0_1_stride_3, &
          bound_global_test_basis_v1_0_1_shape_2 * &
          bound_global_test_basis_v1_0_1_stride_2, &
          bound_global_test_basis_v1_0_1_shape_1 * &
          bound_global_test_basis_v1_0_1_stride_1])
    call C_F_Pointer(bound_global_test_basis_v1_0_2, &
          global_test_basis_v1_0_2, [ &
          bound_global_test_basis_v1_0_2_shape_4 * &
          bound_global_test_basis_v1_0_2_stride_4, &
          bound_global_test_basis_v1_0_2_shape_3 * &
          bound_global_test_basis_v1_0_2_stride_3, &
          bound_global_test_basis_v1_0_2_shape_2 * &
          bound_global_test_basis_v1_0_2_stride_2, &
          bound_global_test_basis_v1_0_2_shape_1 * &
          bound_global_test_basis_v1_0_2_stride_1])
    call C_F_Pointer(bound_global_test_basis_v1_1_1, &
          global_test_basis_v1_1_1, [ &
          bound_global_test_basis_v1_1_1_shape_4 * &
          bound_global_test_basis_v1_1_1_stride_4, &
          bound_global_test_basis_v1_1_1_shape_3 * &
          bound_global_test_basis_v1_1_1_stride_3, &
          bound_global_test_basis_v1_1_1_shape_2 * &
          bound_global_test_basis_v1_1_1_stride_2, &
          bound_global_test_basis_v1_1_1_shape_1 * &
          bound_global_test_basis_v1_1_1_stride_1])
    call C_F_Pointer(bound_global_test_basis_v1_1_2, &
          global_test_basis_v1_1_2, [ &
          bound_global_test_basis_v1_1_2_shape_4 * &
          bound_global_test_basis_v1_1_2_stride_4, &
          bound_global_test_basis_v1_1_2_shape_3 * &
          bound_global_test_basis_v1_1_2_stride_3, &
          bound_global_test_basis_v1_1_2_shape_2 * &
          bound_global_test_basis_v1_1_2_stride_2, &
          bound_global_test_basis_v1_1_2_shape_1 * &
          bound_global_test_basis_v1_1_2_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u1_0_1, &
          global_trial_basis_u1_0_1, [ &
          bound_global_trial_basis_u1_0_1_shape_4 * &
          bound_global_trial_basis_u1_0_1_stride_4, &
          bound_global_trial_basis_u1_0_1_shape_3 * &
          bound_global_trial_basis_u1_0_1_stride_3, &
          bound_global_trial_basis_u1_0_1_shape_2 * &
          bound_global_trial_basis_u1_0_1_stride_2, &
          bound_global_trial_basis_u1_0_1_shape_1 * &
          bound_global_trial_basis_u1_0_1_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u1_0_2, &
          global_trial_basis_u1_0_2, [ &
          bound_global_trial_basis_u1_0_2_shape_4 * &
          bound_global_trial_basis_u1_0_2_stride_4, &
          bound_global_trial_basis_u1_0_2_shape_3 * &
          bound_global_trial_basis_u1_0_2_stride_3, &
          bound_global_trial_basis_u1_0_2_shape_2 * &
          bound_global_trial_basis_u1_0_2_stride_2, &
          bound_global_trial_basis_u1_0_2_shape_1 * &
          bound_global_trial_basis_u1_0_2_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u1_1_1, &
          global_trial_basis_u1_1_1, [ &
          bound_global_trial_basis_u1_1_1_shape_4 * &
          bound_global_trial_basis_u1_1_1_stride_4, &
          bound_global_trial_basis_u1_1_1_shape_3 * &
          bound_global_trial_basis_u1_1_1_stride_3, &
          bound_global_trial_basis_u1_1_1_shape_2 * &
          bound_global_trial_basis_u1_1_1_stride_2, &
          bound_global_trial_basis_u1_1_1_shape_1 * &
          bound_global_trial_basis_u1_1_1_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u1_1_2, &
          global_trial_basis_u1_1_2, [ &
          bound_global_trial_basis_u1_1_2_shape_4 * &
          bound_global_trial_basis_u1_1_2_stride_4, &
          bound_global_trial_basis_u1_1_2_shape_3 * &
          bound_global_trial_basis_u1_1_2_stride_3, &
          bound_global_trial_basis_u1_1_2_shape_2 * &
          bound_global_trial_basis_u1_1_2_stride_2, &
          bound_global_trial_basis_u1_1_2_shape_1 * &
          bound_global_trial_basis_u1_1_2_stride_1])
    call C_F_Pointer(bound_global_span_v1_0_1, global_span_v1_0_1, [ &
          bound_global_span_v1_0_1_shape_1 * &
          bound_global_span_v1_0_1_stride_1])
    call C_F_Pointer(bound_global_span_v1_0_2, global_span_v1_0_2, [ &
          bound_global_span_v1_0_2_shape_1 * &
          bound_global_span_v1_0_2_stride_1])
    call C_F_Pointer(bound_global_span_v1_1_1, global_span_v1_1_1, [ &
          bound_global_span_v1_1_1_shape_1 * &
          bound_global_span_v1_1_1_stride_1])
    call C_F_Pointer(bound_global_span_v1_1_2, global_span_v1_1_2, [ &
          bound_global_span_v1_1_2_shape_1 * &
          bound_global_span_v1_1_2_stride_1])
    call C_F_Pointer(bound_global_x1, global_x1, [ &
          bound_global_x1_shape_2 * bound_global_x1_stride_2, &
          bound_global_x1_shape_1 * bound_global_x1_stride_1])
    call C_F_Pointer(bound_global_x2, global_x2, [ &
          bound_global_x2_shape_2 * bound_global_x2_stride_2, &
          bound_global_x2_shape_1 * bound_global_x2_stride_1])
    call C_F_Pointer(bound_g_mat_u1_0_v1_0_e0el77nq, &
          g_mat_u1_0_v1_0_e0el77nq, [ &
          bound_g_mat_u1_0_v1_0_e0el77nq_shape_4 * &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_4, &
          bound_g_mat_u1_0_v1_0_e0el77nq_shape_3 * &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_3, &
          bound_g_mat_u1_0_v1_0_e0el77nq_shape_2 * &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_2, &
          bound_g_mat_u1_0_v1_0_e0el77nq_shape_1 * &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_1])
    call C_F_Pointer(bound_g_mat_u1_1_v1_0_e0el77nq, &
          g_mat_u1_1_v1_0_e0el77nq, [ &
          bound_g_mat_u1_1_v1_0_e0el77nq_shape_4 * &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_4, &
          bound_g_mat_u1_1_v1_0_e0el77nq_shape_3 * &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_3, &
          bound_g_mat_u1_1_v1_0_e0el77nq_shape_2 * &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_2, &
          bound_g_mat_u1_1_v1_0_e0el77nq_shape_1 * &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_1])
    call C_F_Pointer(bound_g_mat_u1_0_v1_1_e0el77nq, &
          g_mat_u1_0_v1_1_e0el77nq, [ &
          bound_g_mat_u1_0_v1_1_e0el77nq_shape_4 * &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_4, &
          bound_g_mat_u1_0_v1_1_e0el77nq_shape_3 * &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_3, &
          bound_g_mat_u1_0_v1_1_e0el77nq_shape_2 * &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_2, &
          bound_g_mat_u1_0_v1_1_e0el77nq_shape_1 * &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_1])
    call C_F_Pointer(bound_g_mat_u1_1_v1_1_e0el77nq, &
          g_mat_u1_1_v1_1_e0el77nq, [ &
          bound_g_mat_u1_1_v1_1_e0el77nq_shape_4 * &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_4, &
          bound_g_mat_u1_1_v1_1_e0el77nq_shape_3 * &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_3, &
          bound_g_mat_u1_1_v1_1_e0el77nq_shape_2 * &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_2, &
          bound_g_mat_u1_1_v1_1_e0el77nq_shape_1 * &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_1])
    call assemble_matrix_e0el77nq(global_test_basis_v1_0_1 = &
          global_test_basis_v1_0_1(1_i64:: &
          bound_global_test_basis_v1_0_1_stride_4, 1_i64:: &
          bound_global_test_basis_v1_0_1_stride_3, 1_i64:: &
          bound_global_test_basis_v1_0_1_stride_2, 1_i64:: &
          bound_global_test_basis_v1_0_1_stride_1), &
          global_test_basis_v1_0_2 = global_test_basis_v1_0_2(1_i64:: &
          bound_global_test_basis_v1_0_2_stride_4, 1_i64:: &
          bound_global_test_basis_v1_0_2_stride_3, 1_i64:: &
          bound_global_test_basis_v1_0_2_stride_2, 1_i64:: &
          bound_global_test_basis_v1_0_2_stride_1), &
          global_test_basis_v1_1_1 = global_test_basis_v1_1_1(1_i64:: &
          bound_global_test_basis_v1_1_1_stride_4, 1_i64:: &
          bound_global_test_basis_v1_1_1_stride_3, 1_i64:: &
          bound_global_test_basis_v1_1_1_stride_2, 1_i64:: &
          bound_global_test_basis_v1_1_1_stride_1), &
          global_test_basis_v1_1_2 = global_test_basis_v1_1_2(1_i64:: &
          bound_global_test_basis_v1_1_2_stride_4, 1_i64:: &
          bound_global_test_basis_v1_1_2_stride_3, 1_i64:: &
          bound_global_test_basis_v1_1_2_stride_2, 1_i64:: &
          bound_global_test_basis_v1_1_2_stride_1), &
          global_trial_basis_u1_0_1 = global_trial_basis_u1_0_1(1_i64:: &
          bound_global_trial_basis_u1_0_1_stride_4, 1_i64:: &
          bound_global_trial_basis_u1_0_1_stride_3, 1_i64:: &
          bound_global_trial_basis_u1_0_1_stride_2, 1_i64:: &
          bound_global_trial_basis_u1_0_1_stride_1), &
          global_trial_basis_u1_0_2 = global_trial_basis_u1_0_2(1_i64:: &
          bound_global_trial_basis_u1_0_2_stride_4, 1_i64:: &
          bound_global_trial_basis_u1_0_2_stride_3, 1_i64:: &
          bound_global_trial_basis_u1_0_2_stride_2, 1_i64:: &
          bound_global_trial_basis_u1_0_2_stride_1), &
          global_trial_basis_u1_1_1 = global_trial_basis_u1_1_1(1_i64:: &
          bound_global_trial_basis_u1_1_1_stride_4, 1_i64:: &
          bound_global_trial_basis_u1_1_1_stride_3, 1_i64:: &
          bound_global_trial_basis_u1_1_1_stride_2, 1_i64:: &
          bound_global_trial_basis_u1_1_1_stride_1), &
          global_trial_basis_u1_1_2 = global_trial_basis_u1_1_2(1_i64:: &
          bound_global_trial_basis_u1_1_2_stride_4, 1_i64:: &
          bound_global_trial_basis_u1_1_2_stride_3, 1_i64:: &
          bound_global_trial_basis_u1_1_2_stride_2, 1_i64:: &
          bound_global_trial_basis_u1_1_2_stride_1), global_span_v1_0_1 &
          = global_span_v1_0_1(1_i64::bound_global_span_v1_0_1_stride_1 &
          ), global_span_v1_0_2 = global_span_v1_0_2(1_i64:: &
          bound_global_span_v1_0_2_stride_1), global_span_v1_1_1 = &
          global_span_v1_1_1(1_i64::bound_global_span_v1_1_1_stride_1), &
          global_span_v1_1_2 = global_span_v1_1_2(1_i64:: &
          bound_global_span_v1_1_2_stride_1), global_x1 = global_x1( &
          1_i64::bound_global_x1_stride_2, 1_i64:: &
          bound_global_x1_stride_1), global_x2 = global_x2(1_i64:: &
          bound_global_x2_stride_2, 1_i64::bound_global_x2_stride_1), &
          test_v1_0_p1 = test_v1_0_p1, test_v1_0_p2 = test_v1_0_p2, &
          test_v1_1_p1 = test_v1_1_p1, test_v1_1_p2 = test_v1_1_p2, &
          trial_u1_0_p1 = trial_u1_0_p1, trial_u1_0_p2 = trial_u1_0_p2, &
          trial_u1_1_p1 = trial_u1_1_p1, trial_u1_1_p2 = trial_u1_1_p2, &
          n_element_1 = n_element_1, n_element_2 = n_element_2, k1 = k1 &
          , k2 = k2, pad1 = pad1, pad2 = pad2, g_mat_u1_0_v1_0_e0el77nq &
          = g_mat_u1_0_v1_0_e0el77nq(1_i64:: &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_4, 1_i64:: &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_3, 1_i64:: &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_2, 1_i64:: &
          bound_g_mat_u1_0_v1_0_e0el77nq_stride_1), &
          g_mat_u1_1_v1_0_e0el77nq = g_mat_u1_1_v1_0_e0el77nq(1_i64:: &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_4, 1_i64:: &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_3, 1_i64:: &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_2, 1_i64:: &
          bound_g_mat_u1_1_v1_0_e0el77nq_stride_1), &
          g_mat_u1_0_v1_1_e0el77nq = g_mat_u1_0_v1_1_e0el77nq(1_i64:: &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_4, 1_i64:: &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_3, 1_i64:: &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_2, 1_i64:: &
          bound_g_mat_u1_0_v1_1_e0el77nq_stride_1), &
          g_mat_u1_1_v1_1_e0el77nq = g_mat_u1_1_v1_1_e0el77nq(1_i64:: &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_4, 1_i64:: &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_3, 1_i64:: &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_2, 1_i64:: &
          bound_g_mat_u1_1_v1_1_e0el77nq_stride_1))

  end subroutine bind_c_assemble_matrix_e0el77nq
  !........................................

end module bind_c_dependencies_e0el77nq_mx2lvwh9fed2
