module bind_c_dependencies_wti60kr7_6feykyplisif

  use dependencies_wti60kr7_6feykyplisif

  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T , C_F_Pointer , c_ptr
  implicit none

  contains

  !........................................
  subroutine bind_c_assemble_matrix_wti60kr7( &
        bound_global_test_basis_v2_1, &
        bound_global_test_basis_v2_1_shape_1, &
        bound_global_test_basis_v2_1_shape_2, &
        bound_global_test_basis_v2_1_shape_3, &
        bound_global_test_basis_v2_1_shape_4, &
        bound_global_test_basis_v2_1_stride_1, &
        bound_global_test_basis_v2_1_stride_2, &
        bound_global_test_basis_v2_1_stride_3, &
        bound_global_test_basis_v2_1_stride_4, &
        bound_global_test_basis_v2_2, &
        bound_global_test_basis_v2_2_shape_1, &
        bound_global_test_basis_v2_2_shape_2, &
        bound_global_test_basis_v2_2_shape_3, &
        bound_global_test_basis_v2_2_shape_4, &
        bound_global_test_basis_v2_2_stride_1, &
        bound_global_test_basis_v2_2_stride_2, &
        bound_global_test_basis_v2_2_stride_3, &
        bound_global_test_basis_v2_2_stride_4, &
        bound_global_trial_basis_u2_1, &
        bound_global_trial_basis_u2_1_shape_1, &
        bound_global_trial_basis_u2_1_shape_2, &
        bound_global_trial_basis_u2_1_shape_3, &
        bound_global_trial_basis_u2_1_shape_4, &
        bound_global_trial_basis_u2_1_stride_1, &
        bound_global_trial_basis_u2_1_stride_2, &
        bound_global_trial_basis_u2_1_stride_3, &
        bound_global_trial_basis_u2_1_stride_4, &
        bound_global_trial_basis_u2_2, &
        bound_global_trial_basis_u2_2_shape_1, &
        bound_global_trial_basis_u2_2_shape_2, &
        bound_global_trial_basis_u2_2_shape_3, &
        bound_global_trial_basis_u2_2_shape_4, &
        bound_global_trial_basis_u2_2_stride_1, &
        bound_global_trial_basis_u2_2_stride_2, &
        bound_global_trial_basis_u2_2_stride_3, &
        bound_global_trial_basis_u2_2_stride_4, bound_global_span_v2_1, &
        bound_global_span_v2_1_shape_1, bound_global_span_v2_1_stride_1 &
        , bound_global_span_v2_2, bound_global_span_v2_2_shape_1, &
        bound_global_span_v2_2_stride_1, bound_global_x1, &
        bound_global_x1_shape_1, bound_global_x1_shape_2, &
        bound_global_x1_stride_1, bound_global_x1_stride_2, &
        bound_global_x2, bound_global_x2_shape_1, &
        bound_global_x2_shape_2, bound_global_x2_stride_1, &
        bound_global_x2_stride_2, test_v2_p1, test_v2_p2, trial_u2_p1, &
        trial_u2_p2, n_element_1, n_element_2, k1, k2, pad1, pad2, &
        bound_g_mat_u2_v2_wti60kr7, bound_g_mat_u2_v2_wti60kr7_shape_1, &
        bound_g_mat_u2_v2_wti60kr7_shape_2, &
        bound_g_mat_u2_v2_wti60kr7_shape_3, &
        bound_g_mat_u2_v2_wti60kr7_shape_4, &
        bound_g_mat_u2_v2_wti60kr7_stride_1, &
        bound_g_mat_u2_v2_wti60kr7_stride_2, &
        bound_g_mat_u2_v2_wti60kr7_stride_3, &
        bound_g_mat_u2_v2_wti60kr7_stride_4) bind(c)

    implicit none

    type(c_ptr), value :: bound_global_test_basis_v2_1
    integer(i64), value :: bound_global_test_basis_v2_1_shape_1
    integer(i64), value :: bound_global_test_basis_v2_1_shape_2
    integer(i64), value :: bound_global_test_basis_v2_1_shape_3
    integer(i64), value :: bound_global_test_basis_v2_1_shape_4
    integer(i64), value :: bound_global_test_basis_v2_1_stride_1
    integer(i64), value :: bound_global_test_basis_v2_1_stride_2
    integer(i64), value :: bound_global_test_basis_v2_1_stride_3
    integer(i64), value :: bound_global_test_basis_v2_1_stride_4
    type(c_ptr), value :: bound_global_test_basis_v2_2
    integer(i64), value :: bound_global_test_basis_v2_2_shape_1
    integer(i64), value :: bound_global_test_basis_v2_2_shape_2
    integer(i64), value :: bound_global_test_basis_v2_2_shape_3
    integer(i64), value :: bound_global_test_basis_v2_2_shape_4
    integer(i64), value :: bound_global_test_basis_v2_2_stride_1
    integer(i64), value :: bound_global_test_basis_v2_2_stride_2
    integer(i64), value :: bound_global_test_basis_v2_2_stride_3
    integer(i64), value :: bound_global_test_basis_v2_2_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u2_1
    integer(i64), value :: bound_global_trial_basis_u2_1_shape_1
    integer(i64), value :: bound_global_trial_basis_u2_1_shape_2
    integer(i64), value :: bound_global_trial_basis_u2_1_shape_3
    integer(i64), value :: bound_global_trial_basis_u2_1_shape_4
    integer(i64), value :: bound_global_trial_basis_u2_1_stride_1
    integer(i64), value :: bound_global_trial_basis_u2_1_stride_2
    integer(i64), value :: bound_global_trial_basis_u2_1_stride_3
    integer(i64), value :: bound_global_trial_basis_u2_1_stride_4
    type(c_ptr), value :: bound_global_trial_basis_u2_2
    integer(i64), value :: bound_global_trial_basis_u2_2_shape_1
    integer(i64), value :: bound_global_trial_basis_u2_2_shape_2
    integer(i64), value :: bound_global_trial_basis_u2_2_shape_3
    integer(i64), value :: bound_global_trial_basis_u2_2_shape_4
    integer(i64), value :: bound_global_trial_basis_u2_2_stride_1
    integer(i64), value :: bound_global_trial_basis_u2_2_stride_2
    integer(i64), value :: bound_global_trial_basis_u2_2_stride_3
    integer(i64), value :: bound_global_trial_basis_u2_2_stride_4
    type(c_ptr), value :: bound_global_span_v2_1
    integer(i64), value :: bound_global_span_v2_1_shape_1
    integer(i64), value :: bound_global_span_v2_1_stride_1
    type(c_ptr), value :: bound_global_span_v2_2
    integer(i64), value :: bound_global_span_v2_2_shape_1
    integer(i64), value :: bound_global_span_v2_2_stride_1
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
    integer(i64), value :: test_v2_p1
    integer(i64), value :: test_v2_p2
    integer(i64), value :: trial_u2_p1
    integer(i64), value :: trial_u2_p2
    integer(i64), value :: n_element_1
    integer(i64), value :: n_element_2
    integer(i64), value :: k1
    integer(i64), value :: k2
    integer(i64), value :: pad1
    integer(i64), value :: pad2
    type(c_ptr), value :: bound_g_mat_u2_v2_wti60kr7
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_shape_1
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_shape_2
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_shape_3
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_shape_4
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_stride_1
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_stride_2
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_stride_3
    integer(i64), value :: bound_g_mat_u2_v2_wti60kr7_stride_4
    real(f64), pointer :: global_test_basis_v2_1(:,:,:,:)
    real(f64), pointer :: global_test_basis_v2_2(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u2_1(:,:,:,:)
    real(f64), pointer :: global_trial_basis_u2_2(:,:,:,:)
    integer(i64), pointer :: global_span_v2_1(:)
    integer(i64), pointer :: global_span_v2_2(:)
    real(f64), pointer :: global_x1(:,:)
    real(f64), pointer :: global_x2(:,:)
    real(f64), pointer :: g_mat_u2_v2_wti60kr7(:,:,:,:)

    call C_F_Pointer(bound_global_test_basis_v2_1, &
          global_test_basis_v2_1, [bound_global_test_basis_v2_1_shape_4 &
          * bound_global_test_basis_v2_1_stride_4, &
          bound_global_test_basis_v2_1_shape_3 * &
          bound_global_test_basis_v2_1_stride_3, &
          bound_global_test_basis_v2_1_shape_2 * &
          bound_global_test_basis_v2_1_stride_2, &
          bound_global_test_basis_v2_1_shape_1 * &
          bound_global_test_basis_v2_1_stride_1])
    call C_F_Pointer(bound_global_test_basis_v2_2, &
          global_test_basis_v2_2, [bound_global_test_basis_v2_2_shape_4 &
          * bound_global_test_basis_v2_2_stride_4, &
          bound_global_test_basis_v2_2_shape_3 * &
          bound_global_test_basis_v2_2_stride_3, &
          bound_global_test_basis_v2_2_shape_2 * &
          bound_global_test_basis_v2_2_stride_2, &
          bound_global_test_basis_v2_2_shape_1 * &
          bound_global_test_basis_v2_2_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u2_1, &
          global_trial_basis_u2_1, [ &
          bound_global_trial_basis_u2_1_shape_4 * &
          bound_global_trial_basis_u2_1_stride_4, &
          bound_global_trial_basis_u2_1_shape_3 * &
          bound_global_trial_basis_u2_1_stride_3, &
          bound_global_trial_basis_u2_1_shape_2 * &
          bound_global_trial_basis_u2_1_stride_2, &
          bound_global_trial_basis_u2_1_shape_1 * &
          bound_global_trial_basis_u2_1_stride_1])
    call C_F_Pointer(bound_global_trial_basis_u2_2, &
          global_trial_basis_u2_2, [ &
          bound_global_trial_basis_u2_2_shape_4 * &
          bound_global_trial_basis_u2_2_stride_4, &
          bound_global_trial_basis_u2_2_shape_3 * &
          bound_global_trial_basis_u2_2_stride_3, &
          bound_global_trial_basis_u2_2_shape_2 * &
          bound_global_trial_basis_u2_2_stride_2, &
          bound_global_trial_basis_u2_2_shape_1 * &
          bound_global_trial_basis_u2_2_stride_1])
    call C_F_Pointer(bound_global_span_v2_1, global_span_v2_1, [ &
          bound_global_span_v2_1_shape_1 * &
          bound_global_span_v2_1_stride_1])
    call C_F_Pointer(bound_global_span_v2_2, global_span_v2_2, [ &
          bound_global_span_v2_2_shape_1 * &
          bound_global_span_v2_2_stride_1])
    call C_F_Pointer(bound_global_x1, global_x1, [ &
          bound_global_x1_shape_2 * bound_global_x1_stride_2, &
          bound_global_x1_shape_1 * bound_global_x1_stride_1])
    call C_F_Pointer(bound_global_x2, global_x2, [ &
          bound_global_x2_shape_2 * bound_global_x2_stride_2, &
          bound_global_x2_shape_1 * bound_global_x2_stride_1])
    call C_F_Pointer(bound_g_mat_u2_v2_wti60kr7, g_mat_u2_v2_wti60kr7, [ &
          bound_g_mat_u2_v2_wti60kr7_shape_4 * &
          bound_g_mat_u2_v2_wti60kr7_stride_4, &
          bound_g_mat_u2_v2_wti60kr7_shape_3 * &
          bound_g_mat_u2_v2_wti60kr7_stride_3, &
          bound_g_mat_u2_v2_wti60kr7_shape_2 * &
          bound_g_mat_u2_v2_wti60kr7_stride_2, &
          bound_g_mat_u2_v2_wti60kr7_shape_1 * &
          bound_g_mat_u2_v2_wti60kr7_stride_1])
    call assemble_matrix_wti60kr7(global_test_basis_v2_1 = &
          global_test_basis_v2_1(1_i64:: &
          bound_global_test_basis_v2_1_stride_4, 1_i64:: &
          bound_global_test_basis_v2_1_stride_3, 1_i64:: &
          bound_global_test_basis_v2_1_stride_2, 1_i64:: &
          bound_global_test_basis_v2_1_stride_1), &
          global_test_basis_v2_2 = global_test_basis_v2_2(1_i64:: &
          bound_global_test_basis_v2_2_stride_4, 1_i64:: &
          bound_global_test_basis_v2_2_stride_3, 1_i64:: &
          bound_global_test_basis_v2_2_stride_2, 1_i64:: &
          bound_global_test_basis_v2_2_stride_1), &
          global_trial_basis_u2_1 = global_trial_basis_u2_1(1_i64:: &
          bound_global_trial_basis_u2_1_stride_4, 1_i64:: &
          bound_global_trial_basis_u2_1_stride_3, 1_i64:: &
          bound_global_trial_basis_u2_1_stride_2, 1_i64:: &
          bound_global_trial_basis_u2_1_stride_1), &
          global_trial_basis_u2_2 = global_trial_basis_u2_2(1_i64:: &
          bound_global_trial_basis_u2_2_stride_4, 1_i64:: &
          bound_global_trial_basis_u2_2_stride_3, 1_i64:: &
          bound_global_trial_basis_u2_2_stride_2, 1_i64:: &
          bound_global_trial_basis_u2_2_stride_1), global_span_v2_1 = &
          global_span_v2_1(1_i64::bound_global_span_v2_1_stride_1), &
          global_span_v2_2 = global_span_v2_2(1_i64:: &
          bound_global_span_v2_2_stride_1), global_x1 = global_x1(1_i64 &
          ::bound_global_x1_stride_2, 1_i64::bound_global_x1_stride_1), &
          global_x2 = global_x2(1_i64::bound_global_x2_stride_2, 1_i64 &
          ::bound_global_x2_stride_1), test_v2_p1 = test_v2_p1, &
          test_v2_p2 = test_v2_p2, trial_u2_p1 = trial_u2_p1, &
          trial_u2_p2 = trial_u2_p2, n_element_1 = n_element_1, &
          n_element_2 = n_element_2, k1 = k1, k2 = k2, pad1 = pad1, &
          pad2 = pad2, g_mat_u2_v2_wti60kr7 = g_mat_u2_v2_wti60kr7( &
          1_i64::bound_g_mat_u2_v2_wti60kr7_stride_4, 1_i64:: &
          bound_g_mat_u2_v2_wti60kr7_stride_3, 1_i64:: &
          bound_g_mat_u2_v2_wti60kr7_stride_2, 1_i64:: &
          bound_g_mat_u2_v2_wti60kr7_stride_1))

  end subroutine bind_c_assemble_matrix_wti60kr7
  !........................................

end module bind_c_dependencies_wti60kr7_6feykyplisif
