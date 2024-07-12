module dependencies_e0el77nq_mx2lvwh9fed2


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine assemble_matrix_e0el77nq(global_test_basis_v1_0_1, &
        global_test_basis_v1_0_2, global_test_basis_v1_1_1, &
        global_test_basis_v1_1_2, global_trial_basis_u1_0_1, &
        global_trial_basis_u1_0_2, global_trial_basis_u1_1_1, &
        global_trial_basis_u1_1_2, global_span_v1_0_1, &
        global_span_v1_0_2, global_span_v1_1_1, global_span_v1_1_2, &
        global_x1, global_x2, test_v1_0_p1, test_v1_0_p2, test_v1_1_p1, &
        test_v1_1_p2, trial_u1_0_p1, trial_u1_0_p2, trial_u1_1_p1, &
        trial_u1_1_p2, n_element_1, n_element_2, k1, k2, pad1, pad2, &
        g_mat_u1_0_v1_0_e0el77nq, g_mat_u1_1_v1_0_e0el77nq, &
        g_mat_u1_0_v1_1_e0el77nq, g_mat_u1_1_v1_1_e0el77nq)

    implicit none

    real(f64), intent(in) :: global_test_basis_v1_0_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_test_basis_v1_0_2(0:,0:,0:,0:)
    real(f64), intent(in) :: global_test_basis_v1_1_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_test_basis_v1_1_2(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u1_0_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u1_0_2(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u1_1_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u1_1_2(0:,0:,0:,0:)
    integer(i64), intent(in) :: global_span_v1_0_1(0:)
    integer(i64), intent(in) :: global_span_v1_0_2(0:)
    integer(i64), intent(in) :: global_span_v1_1_1(0:)
    integer(i64), intent(in) :: global_span_v1_1_2(0:)
    real(f64), intent(in) :: global_x1(0:,0:)
    real(f64), intent(in) :: global_x2(0:,0:)
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
    real(f64), intent(inout) :: g_mat_u1_0_v1_0_e0el77nq(0:,0:,0:,0:)
    real(f64), intent(inout) :: g_mat_u1_1_v1_0_e0el77nq(0:,0:,0:,0:)
    real(f64), intent(inout) :: g_mat_u1_0_v1_1_e0el77nq(0:,0:,0:,0:)
    real(f64), intent(inout) :: g_mat_u1_1_v1_1_e0el77nq(0:,0:,0:,0:)
    real(f64), allocatable :: local_x1(:)
    real(f64), allocatable :: local_x2(:)
    real(f64), allocatable :: l_mat_u1_0_v1_0_e0el77nq(:,:,:,:)
    real(f64), allocatable :: l_mat_u1_0_v1_1_e0el77nq(:,:,:,:)
    real(f64), allocatable :: l_mat_u1_1_v1_0_e0el77nq(:,:,:,:)
    real(f64), allocatable :: l_mat_u1_1_v1_1_e0el77nq(:,:,:,:)
    integer(i64) :: i_element_1
    integer(i64) :: span_v1_0_1
    integer(i64) :: span_v1_1_1
    integer(i64) :: i_element_2
    integer(i64) :: span_v1_0_2
    integer(i64) :: span_v1_1_2
    integer(i64) :: i_quad_1
    real(f64) :: x1
    integer(i64) :: i_quad_2
    real(f64) :: x2
    integer(i64) :: i_basis_1
    integer(i64) :: i_basis_2
    integer(i64) :: j_basis_1
    real(f64) :: v1_0_1
    real(f64) :: v1_0_1_x1
    real(f64) :: u1_0_1
    real(f64) :: u1_0_1_x1
    integer(i64) :: j_basis_2
    real(f64) :: v1_0_2
    real(f64) :: v1_0_2_x2
    real(f64) :: u1_0_2
    real(f64) :: u1_0_2_x2
    real(f64) :: v1_0
    real(f64) :: v1_0_x2
    real(f64) :: v1_0_x1
    real(f64) :: u1_0
    real(f64) :: u1_0_x2
    real(f64) :: u1_0_x1
    real(f64) :: temp_v1_0_u1_0_0
    real(f64) :: temp_v1_0_u1_0_1
    real(f64) :: temp_v1_0_u1_0_2
    real(f64) :: temp_v1_0_u1_0_3
    real(f64) :: temp_v1_0_u1_0_4
    real(f64) :: temp_v1_0_u1_0_5
    real(f64) :: temp_v1_0_u1_0_6
    real(f64) :: temp_v1_0_u1_0_7
    real(f64) :: temp_v1_0_u1_0_8
    real(f64) :: contribution_v1_0_u1_0_e0el77nq
    real(f64) :: u1_1_1
    real(f64) :: u1_1_1_x1
    real(f64) :: u1_1_2
    real(f64) :: u1_1_2_x2
    real(f64) :: u1_1
    real(f64) :: u1_1_x2
    real(f64) :: u1_1_x1
    real(f64) :: temp_v1_0_u1_1_0
    real(f64) :: temp_v1_0_u1_1_1
    real(f64) :: temp_v1_0_u1_1_2
    real(f64) :: temp_v1_0_u1_1_3
    real(f64) :: temp_v1_0_u1_1_4
    real(f64) :: temp_v1_0_u1_1_5
    real(f64) :: temp_v1_0_u1_1_6
    real(f64) :: temp_v1_0_u1_1_7
    real(f64) :: temp_v1_0_u1_1_8
    real(f64) :: temp_v1_0_u1_1_9
    real(f64) :: contribution_v1_0_u1_1_e0el77nq
    real(f64) :: v1_1_1
    real(f64) :: v1_1_1_x1
    real(f64) :: v1_1_2
    real(f64) :: v1_1_2_x2
    real(f64) :: v1_1
    real(f64) :: v1_1_x2
    real(f64) :: v1_1_x1
    real(f64) :: temp_v1_1_u1_0_0
    real(f64) :: temp_v1_1_u1_0_1
    real(f64) :: temp_v1_1_u1_0_2
    real(f64) :: temp_v1_1_u1_0_3
    real(f64) :: temp_v1_1_u1_0_4
    real(f64) :: temp_v1_1_u1_0_5
    real(f64) :: temp_v1_1_u1_0_6
    real(f64) :: temp_v1_1_u1_0_7
    real(f64) :: temp_v1_1_u1_0_8
    real(f64) :: temp_v1_1_u1_0_9
    real(f64) :: contribution_v1_1_u1_0_e0el77nq
    real(f64) :: temp_v1_1_u1_1_0
    real(f64) :: temp_v1_1_u1_1_1
    real(f64) :: temp_v1_1_u1_1_2
    real(f64) :: temp_v1_1_u1_1_3
    real(f64) :: temp_v1_1_u1_1_4
    real(f64) :: temp_v1_1_u1_1_5
    real(f64) :: temp_v1_1_u1_1_6
    real(f64) :: temp_v1_1_u1_1_7
    real(f64) :: temp_v1_1_u1_1_8
    real(f64) :: contribution_v1_1_u1_1_e0el77nq

    allocate(local_x1(0:size(global_x1, 1_i64, i64) - 1_i64))
    local_x1 = 0.0_f64
    allocate(local_x2(0:size(global_x2, 1_i64, i64) - 1_i64))
    local_x2 = 0.0_f64
    allocate(l_mat_u1_0_v1_0_e0el77nq(0:6_i64, 0:4_i64, 0:3_i64, 0:2_i64 &
          ))
    l_mat_u1_0_v1_0_e0el77nq = 0.0_f64
    allocate(l_mat_u1_0_v1_1_e0el77nq(0:6_i64, 0:6_i64, 0:2_i64, 0:3_i64 &
          ))
    l_mat_u1_0_v1_1_e0el77nq = 0.0_f64
    allocate(l_mat_u1_1_v1_0_e0el77nq(0:6_i64, 0:6_i64, 0:3_i64, 0:2_i64 &
          ))
    l_mat_u1_1_v1_0_e0el77nq = 0.0_f64
    allocate(l_mat_u1_1_v1_1_e0el77nq(0:4_i64, 0:6_i64, 0:2_i64, 0:3_i64 &
          ))
    l_mat_u1_1_v1_1_e0el77nq = 0.0_f64
    do i_element_1 = 0_i64, n_element_1 - 1_i64
      local_x1(:) = global_x1(:, i_element_1)
      span_v1_0_1 = global_span_v1_0_1(i_element_1)
      span_v1_1_1 = global_span_v1_1_1(i_element_1)
      do i_element_2 = 0_i64, n_element_2 - 1_i64
        local_x2(:) = global_x2(:, i_element_2)
        span_v1_0_2 = global_span_v1_0_2(i_element_2)
        span_v1_1_2 = global_span_v1_1_2(i_element_2)
        l_mat_u1_0_v1_0_e0el77nq(:, :, :, :) = 0.0_f64
        l_mat_u1_1_v1_0_e0el77nq(:, :, :, :) = 0.0_f64
        l_mat_u1_0_v1_1_e0el77nq(:, :, :, :) = 0.0_f64
        l_mat_u1_1_v1_1_e0el77nq(:, :, :, :) = 0.0_f64
        do i_quad_1 = 0_i64, 3_i64
          x1 = local_x1(i_quad_1)
          do i_quad_2 = 0_i64, 3_i64
            x2 = local_x2(i_quad_2)
            do i_basis_1 = 0_i64, 2_i64
              do i_basis_2 = 0_i64, 3_i64
                do j_basis_1 = 0_i64, 2_i64
                  v1_0_1 = global_test_basis_v1_0_1(i_quad_1, 0_i64, &
                        i_basis_1, i_element_1)
                  v1_0_1_x1 = global_test_basis_v1_0_1(i_quad_1, 1_i64, &
                        i_basis_1, i_element_1)
                  u1_0_1 = global_trial_basis_u1_0_1(i_quad_1, 0_i64, &
                        j_basis_1, i_element_1)
                  u1_0_1_x1 = global_trial_basis_u1_0_1(i_quad_1, 1_i64, &
                        j_basis_1, i_element_1)
                  do j_basis_2 = 0_i64, 3_i64
                    v1_0_2 = global_test_basis_v1_0_2(i_quad_2, 0_i64, &
                          i_basis_2, i_element_2)
                    v1_0_2_x2 = global_test_basis_v1_0_2(i_quad_2, 1_i64 &
                          , i_basis_2, i_element_2)
                    u1_0_2 = global_trial_basis_u1_0_2(i_quad_2, 0_i64, &
                          j_basis_2, i_element_2)
                    u1_0_2_x2 = global_trial_basis_u1_0_2(i_quad_2, &
                          1_i64, j_basis_2, i_element_2)
                    v1_0 = v1_0_1 * v1_0_2
                    v1_0_x2 = v1_0_1 * v1_0_2_x2
                    v1_0_x1 = v1_0_1_x1 * v1_0_2
                    u1_0 = u1_0_1 * u1_0_2
                    u1_0_x2 = u1_0_1 * u1_0_2_x2
                    u1_0_x1 = u1_0_1_x1 * u1_0_2
                    temp_v1_0_u1_0_0 = 2_i64 * 3.141592653589793_f64
                    temp_v1_0_u1_0_1 = temp_v1_0_u1_0_0 * x2
                    temp_v1_0_u1_0_2 = temp_v1_0_u1_0_0 * x1
                    temp_v1_0_u1_0_3 = 0.25_f64 * sin(temp_v1_0_u1_0_1) &
                          * cos(temp_v1_0_u1_0_2)
                    temp_v1_0_u1_0_4 = sin(temp_v1_0_u1_0_2)
                    temp_v1_0_u1_0_5 = cos(temp_v1_0_u1_0_1)
                    temp_v1_0_u1_0_6 = 0.25_f64 * temp_v1_0_u1_0_4 * &
                          temp_v1_0_u1_0_5
                    temp_v1_0_u1_0_7 = temp_v1_0_u1_0_6 + 1.0_f64
                    temp_v1_0_u1_0_8 = u1_0 * v1_0 / (temp_v1_0_u1_0_3 + &
                          temp_v1_0_u1_0_6 + 1_i64) ** 2_i64
                    contribution_v1_0_u1_0_e0el77nq = 4.0_f64 * ( &
                          0.015625_f64 * (temp_v1_0_u1_0_4 * &
                          temp_v1_0_u1_0_4) * (temp_v1_0_u1_0_5 * &
                          temp_v1_0_u1_0_5) * temp_v1_0_u1_0_8 + &
                          0.25_f64 * (temp_v1_0_u1_0_7 * &
                          temp_v1_0_u1_0_7) * temp_v1_0_u1_0_8) * abs(( &
                          temp_v1_0_u1_0_3 + temp_v1_0_u1_0_7))
                    l_mat_u1_0_v1_0_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 2_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) = &
                          l_mat_u1_0_v1_0_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 2_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) + &
                          contribution_v1_0_u1_0_e0el77nq
                  end do
                end do
              end do
            end do
            do i_basis_1 = 0_i64, 2_i64
              do i_basis_2 = 0_i64, 3_i64
                do j_basis_1 = 0_i64, 3_i64
                  v1_0_1 = global_test_basis_v1_0_1(i_quad_1, 0_i64, &
                        i_basis_1, i_element_1)
                  v1_0_1_x1 = global_test_basis_v1_0_1(i_quad_1, 1_i64, &
                        i_basis_1, i_element_1)
                  u1_1_1 = global_trial_basis_u1_1_1(i_quad_1, 0_i64, &
                        j_basis_1, i_element_1)
                  u1_1_1_x1 = global_trial_basis_u1_1_1(i_quad_1, 1_i64, &
                        j_basis_1, i_element_1)
                  do j_basis_2 = 0_i64, 2_i64
                    v1_0_2 = global_test_basis_v1_0_2(i_quad_2, 0_i64, &
                          i_basis_2, i_element_2)
                    v1_0_2_x2 = global_test_basis_v1_0_2(i_quad_2, 1_i64 &
                          , i_basis_2, i_element_2)
                    u1_1_2 = global_trial_basis_u1_1_2(i_quad_2, 0_i64, &
                          j_basis_2, i_element_2)
                    u1_1_2_x2 = global_trial_basis_u1_1_2(i_quad_2, &
                          1_i64, j_basis_2, i_element_2)
                    v1_0 = v1_0_1 * v1_0_2
                    v1_0_x2 = v1_0_1 * v1_0_2_x2
                    v1_0_x1 = v1_0_1_x1 * v1_0_2
                    u1_1 = u1_1_1 * u1_1_2
                    u1_1_x2 = u1_1_1 * u1_1_2_x2
                    u1_1_x1 = u1_1_1_x1 * u1_1_2
                    temp_v1_0_u1_1_0 = 2_i64 * 3.141592653589793_f64
                    temp_v1_0_u1_1_1 = temp_v1_0_u1_1_0 * x1
                    temp_v1_0_u1_1_2 = temp_v1_0_u1_1_0 * x2
                    temp_v1_0_u1_1_3 = sin(temp_v1_0_u1_1_2) * cos( &
                          temp_v1_0_u1_1_1)
                    temp_v1_0_u1_1_4 = sin(temp_v1_0_u1_1_1) * cos( &
                          temp_v1_0_u1_1_2)
                    temp_v1_0_u1_1_5 = 0.25_f64 * temp_v1_0_u1_1_4
                    temp_v1_0_u1_1_6 = temp_v1_0_u1_1_5 + 1.0_f64
                    temp_v1_0_u1_1_7 = 0.5_f64 * temp_v1_0_u1_1_3
                    temp_v1_0_u1_1_8 = temp_v1_0_u1_1_7 + 2.0_f64
                    temp_v1_0_u1_1_9 = u1_1 * v1_0 / ((0.5_f64 * &
                          temp_v1_0_u1_1_4 + temp_v1_0_u1_1_8) * ( &
                          1.0_f64 * temp_v1_0_u1_1_3 + 1.0_f64 * &
                          temp_v1_0_u1_1_4 + 4.0_f64))
                    contribution_v1_0_u1_1_e0el77nq = 4.0_f64 * (( &
                          -temp_v1_0_u1_1_5) * temp_v1_0_u1_1_8 * &
                          temp_v1_0_u1_1_9 - temp_v1_0_u1_1_6 * &
                          temp_v1_0_u1_1_7 * temp_v1_0_u1_1_9) * abs(( &
                          0.25_f64 * temp_v1_0_u1_1_3 + &
                          temp_v1_0_u1_1_6))
                    l_mat_u1_1_v1_0_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) = &
                          l_mat_u1_1_v1_0_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) + &
                          contribution_v1_0_u1_1_e0el77nq
                  end do
                end do
              end do
            end do
            do i_basis_1 = 0_i64, 3_i64
              do i_basis_2 = 0_i64, 2_i64
                do j_basis_1 = 0_i64, 2_i64
                  v1_1_1 = global_test_basis_v1_1_1(i_quad_1, 0_i64, &
                        i_basis_1, i_element_1)
                  v1_1_1_x1 = global_test_basis_v1_1_1(i_quad_1, 1_i64, &
                        i_basis_1, i_element_1)
                  u1_0_1 = global_trial_basis_u1_0_1(i_quad_1, 0_i64, &
                        j_basis_1, i_element_1)
                  u1_0_1_x1 = global_trial_basis_u1_0_1(i_quad_1, 1_i64, &
                        j_basis_1, i_element_1)
                  do j_basis_2 = 0_i64, 3_i64
                    v1_1_2 = global_test_basis_v1_1_2(i_quad_2, 0_i64, &
                          i_basis_2, i_element_2)
                    v1_1_2_x2 = global_test_basis_v1_1_2(i_quad_2, 1_i64 &
                          , i_basis_2, i_element_2)
                    u1_0_2 = global_trial_basis_u1_0_2(i_quad_2, 0_i64, &
                          j_basis_2, i_element_2)
                    u1_0_2_x2 = global_trial_basis_u1_0_2(i_quad_2, &
                          1_i64, j_basis_2, i_element_2)
                    v1_1 = v1_1_1 * v1_1_2
                    v1_1_x2 = v1_1_1 * v1_1_2_x2
                    v1_1_x1 = v1_1_1_x1 * v1_1_2
                    u1_0 = u1_0_1 * u1_0_2
                    u1_0_x2 = u1_0_1 * u1_0_2_x2
                    u1_0_x1 = u1_0_1_x1 * u1_0_2
                    temp_v1_1_u1_0_0 = 2_i64 * 3.141592653589793_f64
                    temp_v1_1_u1_0_1 = temp_v1_1_u1_0_0 * x1
                    temp_v1_1_u1_0_2 = temp_v1_1_u1_0_0 * x2
                    temp_v1_1_u1_0_3 = sin(temp_v1_1_u1_0_2) * cos( &
                          temp_v1_1_u1_0_1)
                    temp_v1_1_u1_0_4 = sin(temp_v1_1_u1_0_1) * cos( &
                          temp_v1_1_u1_0_2)
                    temp_v1_1_u1_0_5 = 0.25_f64 * temp_v1_1_u1_0_4
                    temp_v1_1_u1_0_6 = temp_v1_1_u1_0_5 + 1.0_f64
                    temp_v1_1_u1_0_7 = 0.5_f64 * temp_v1_1_u1_0_3
                    temp_v1_1_u1_0_8 = temp_v1_1_u1_0_7 + 2.0_f64
                    temp_v1_1_u1_0_9 = u1_0 * v1_1 / ((0.5_f64 * &
                          temp_v1_1_u1_0_4 + temp_v1_1_u1_0_8) * ( &
                          1.0_f64 * temp_v1_1_u1_0_3 + 1.0_f64 * &
                          temp_v1_1_u1_0_4 + 4.0_f64))
                    contribution_v1_1_u1_0_e0el77nq = 4.0_f64 * (( &
                          -temp_v1_1_u1_0_5) * temp_v1_1_u1_0_8 * &
                          temp_v1_1_u1_0_9 - temp_v1_1_u1_0_6 * &
                          temp_v1_1_u1_0_7 * temp_v1_1_u1_0_9) * abs(( &
                          0.25_f64 * temp_v1_1_u1_0_3 + &
                          temp_v1_1_u1_0_6))
                    l_mat_u1_0_v1_1_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) = &
                          l_mat_u1_0_v1_1_e0el77nq(3_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) + &
                          contribution_v1_1_u1_0_e0el77nq
                  end do
                end do
              end do
            end do
            do i_basis_1 = 0_i64, 3_i64
              do i_basis_2 = 0_i64, 2_i64
                do j_basis_1 = 0_i64, 3_i64
                  v1_1_1 = global_test_basis_v1_1_1(i_quad_1, 0_i64, &
                        i_basis_1, i_element_1)
                  v1_1_1_x1 = global_test_basis_v1_1_1(i_quad_1, 1_i64, &
                        i_basis_1, i_element_1)
                  u1_1_1 = global_trial_basis_u1_1_1(i_quad_1, 0_i64, &
                        j_basis_1, i_element_1)
                  u1_1_1_x1 = global_trial_basis_u1_1_1(i_quad_1, 1_i64, &
                        j_basis_1, i_element_1)
                  do j_basis_2 = 0_i64, 2_i64
                    v1_1_2 = global_test_basis_v1_1_2(i_quad_2, 0_i64, &
                          i_basis_2, i_element_2)
                    v1_1_2_x2 = global_test_basis_v1_1_2(i_quad_2, 1_i64 &
                          , i_basis_2, i_element_2)
                    u1_1_2 = global_trial_basis_u1_1_2(i_quad_2, 0_i64, &
                          j_basis_2, i_element_2)
                    u1_1_2_x2 = global_trial_basis_u1_1_2(i_quad_2, &
                          1_i64, j_basis_2, i_element_2)
                    v1_1 = v1_1_1 * v1_1_2
                    v1_1_x2 = v1_1_1 * v1_1_2_x2
                    v1_1_x1 = v1_1_1_x1 * v1_1_2
                    u1_1 = u1_1_1 * u1_1_2
                    u1_1_x2 = u1_1_1 * u1_1_2_x2
                    u1_1_x1 = u1_1_1_x1 * u1_1_2
                    temp_v1_1_u1_1_0 = 2_i64 * 3.141592653589793_f64
                    temp_v1_1_u1_1_1 = temp_v1_1_u1_1_0 * x1
                    temp_v1_1_u1_1_2 = temp_v1_1_u1_1_0 * x2
                    temp_v1_1_u1_1_3 = 0.25_f64 * sin(temp_v1_1_u1_1_1) &
                          * cos(temp_v1_1_u1_1_2)
                    temp_v1_1_u1_1_4 = sin(temp_v1_1_u1_1_2)
                    temp_v1_1_u1_1_5 = cos(temp_v1_1_u1_1_1)
                    temp_v1_1_u1_1_6 = 0.25_f64 * temp_v1_1_u1_1_4 * &
                          temp_v1_1_u1_1_5
                    temp_v1_1_u1_1_7 = temp_v1_1_u1_1_6 + 1_i64
                    temp_v1_1_u1_1_8 = u1_1 * v1_1 / (temp_v1_1_u1_1_3 + &
                          temp_v1_1_u1_1_7) ** 2_i64
                    contribution_v1_1_u1_1_e0el77nq = 4.0_f64 * ( &
                          0.015625_f64 * (temp_v1_1_u1_1_4 * &
                          temp_v1_1_u1_1_4) * (temp_v1_1_u1_1_5 * &
                          temp_v1_1_u1_1_5) * temp_v1_1_u1_1_8 + &
                          0.25_f64 * (temp_v1_1_u1_1_7 * &
                          temp_v1_1_u1_1_7) * temp_v1_1_u1_1_8) * abs(( &
                          temp_v1_1_u1_1_3 + temp_v1_1_u1_1_6 + 1.0_f64 &
                          ))
                    l_mat_u1_1_v1_1_e0el77nq(2_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) = &
                          l_mat_u1_1_v1_1_e0el77nq(2_i64 - i_basis_2 + &
                          j_basis_2, 3_i64 - i_basis_1 + j_basis_1, &
                          i_basis_2, i_basis_1) + &
                          contribution_v1_1_u1_1_e0el77nq
                  end do
                end do
              end do
            end do
          end do
        end do
        g_mat_u1_0_v1_0_e0el77nq(:, :, pad2 + span_v1_0_2 - test_v1_0_p2 &
              :1_i64 + pad2 + span_v1_0_2 - 1_i64, pad1 + span_v1_0_1 - &
              test_v1_0_p1:1_i64 + pad1 + span_v1_0_1 - 1_i64) = &
              g_mat_u1_0_v1_0_e0el77nq(:, :, pad2 + span_v1_0_2 - &
              test_v1_0_p2:1_i64 + pad2 + span_v1_0_2 - 1_i64, pad1 + &
              span_v1_0_1 - test_v1_0_p1:1_i64 + pad1 + span_v1_0_1 - &
              1_i64) + l_mat_u1_0_v1_0_e0el77nq(:, :, :, :)
        g_mat_u1_1_v1_0_e0el77nq(:, :, pad2 + span_v1_0_2 - test_v1_0_p2 &
              :1_i64 + pad2 + span_v1_0_2 - 1_i64, pad1 + span_v1_0_1 - &
              test_v1_0_p1:1_i64 + pad1 + span_v1_0_1 - 1_i64) = &
              g_mat_u1_1_v1_0_e0el77nq(:, :, pad2 + span_v1_0_2 - &
              test_v1_0_p2:1_i64 + pad2 + span_v1_0_2 - 1_i64, pad1 + &
              span_v1_0_1 - test_v1_0_p1:1_i64 + pad1 + span_v1_0_1 - &
              1_i64) + l_mat_u1_1_v1_0_e0el77nq(:, :, :, :)
        g_mat_u1_0_v1_1_e0el77nq(:, :, pad2 + span_v1_1_2 - test_v1_1_p2 &
              :1_i64 + pad2 + span_v1_1_2 - 1_i64, pad1 + span_v1_1_1 - &
              test_v1_1_p1:1_i64 + pad1 + span_v1_1_1 - 1_i64) = &
              g_mat_u1_0_v1_1_e0el77nq(:, :, pad2 + span_v1_1_2 - &
              test_v1_1_p2:1_i64 + pad2 + span_v1_1_2 - 1_i64, pad1 + &
              span_v1_1_1 - test_v1_1_p1:1_i64 + pad1 + span_v1_1_1 - &
              1_i64) + l_mat_u1_0_v1_1_e0el77nq(:, :, :, :)
        g_mat_u1_1_v1_1_e0el77nq(:, :, pad2 + span_v1_1_2 - test_v1_1_p2 &
              :1_i64 + pad2 + span_v1_1_2 - 1_i64, pad1 + span_v1_1_1 - &
              test_v1_1_p1:1_i64 + pad1 + span_v1_1_1 - 1_i64) = &
              g_mat_u1_1_v1_1_e0el77nq(:, :, pad2 + span_v1_1_2 - &
              test_v1_1_p2:1_i64 + pad2 + span_v1_1_2 - 1_i64, pad1 + &
              span_v1_1_1 - test_v1_1_p1:1_i64 + pad1 + span_v1_1_1 - &
              1_i64) + l_mat_u1_1_v1_1_e0el77nq(:, :, :, :)
      end do
    end do
    if (allocated(l_mat_u1_0_v1_0_e0el77nq)) then
      deallocate(l_mat_u1_0_v1_0_e0el77nq)
    end if
    if (allocated(local_x1)) then
      deallocate(local_x1)
    end if
    if (allocated(l_mat_u1_0_v1_1_e0el77nq)) then
      deallocate(l_mat_u1_0_v1_1_e0el77nq)
    end if
    if (allocated(l_mat_u1_1_v1_0_e0el77nq)) then
      deallocate(l_mat_u1_1_v1_0_e0el77nq)
    end if
    if (allocated(l_mat_u1_1_v1_1_e0el77nq)) then
      deallocate(l_mat_u1_1_v1_1_e0el77nq)
    end if
    if (allocated(local_x2)) then
      deallocate(local_x2)
    end if
    return

  end subroutine assemble_matrix_e0el77nq
  !........................................

end module dependencies_e0el77nq_mx2lvwh9fed2
