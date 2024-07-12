module dependencies_wti60kr7_6feykyplisif


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine assemble_matrix_wti60kr7(global_test_basis_v2_1, &
        global_test_basis_v2_2, global_trial_basis_u2_1, &
        global_trial_basis_u2_2, global_span_v2_1, global_span_v2_2, &
        global_x1, global_x2, test_v2_p1, test_v2_p2, trial_u2_p1, &
        trial_u2_p2, n_element_1, n_element_2, k1, k2, pad1, pad2, &
        g_mat_u2_v2_wti60kr7)

    implicit none

    real(f64), intent(in) :: global_test_basis_v2_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_test_basis_v2_2(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u2_1(0:,0:,0:,0:)
    real(f64), intent(in) :: global_trial_basis_u2_2(0:,0:,0:,0:)
    integer(i64), intent(in) :: global_span_v2_1(0:)
    integer(i64), intent(in) :: global_span_v2_2(0:)
    real(f64), intent(in) :: global_x1(0:,0:)
    real(f64), intent(in) :: global_x2(0:,0:)
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
    real(f64), intent(inout) :: g_mat_u2_v2_wti60kr7(0:,0:,0:,0:)
    real(f64), allocatable :: local_x1(:)
    real(f64), allocatable :: local_x2(:)
    real(f64), allocatable :: l_mat_u2_v2_wti60kr7(:,:,:,:)
    integer(i64) :: i_element_1
    integer(i64) :: span_v2_1
    integer(i64) :: i_element_2
    integer(i64) :: span_v2_2
    integer(i64) :: i_quad_1
    real(f64) :: x1
    integer(i64) :: i_quad_2
    real(f64) :: x2
    integer(i64) :: i_basis_1
    integer(i64) :: i_basis_2
    integer(i64) :: j_basis_1
    real(f64) :: v2_1
    real(f64) :: v2_1_x1
    real(f64) :: u2_1
    real(f64) :: u2_1_x1
    integer(i64) :: j_basis_2
    real(f64) :: v2_2
    real(f64) :: v2_2_x2
    real(f64) :: u2_2
    real(f64) :: u2_2_x2
    real(f64) :: v2
    real(f64) :: v2_x2
    real(f64) :: v2_x1
    real(f64) :: u2
    real(f64) :: u2_x2
    real(f64) :: u2_x1
    real(f64) :: temp_v2_u2_0
    real(f64) :: temp_v2_u2_1
    real(f64) :: temp_v2_u2_2
    real(f64) :: temp_v2_u2_3
    real(f64) :: contribution_v2_u2_wti60kr7

    allocate(local_x1(0:size(global_x1, 1_i64, i64) - 1_i64))
    local_x1 = 0.0_f64
    allocate(local_x2(0:size(global_x2, 1_i64, i64) - 1_i64))
    local_x2 = 0.0_f64
    allocate(l_mat_u2_v2_wti60kr7(0:4_i64, 0:4_i64, 0:2_i64, 0:2_i64))
    l_mat_u2_v2_wti60kr7 = 0.0_f64
    do i_element_1 = 0_i64, n_element_1 - 1_i64
      local_x1(:) = global_x1(:, i_element_1)
      span_v2_1 = global_span_v2_1(i_element_1)
      do i_element_2 = 0_i64, n_element_2 - 1_i64
        local_x2(:) = global_x2(:, i_element_2)
        span_v2_2 = global_span_v2_2(i_element_2)
        l_mat_u2_v2_wti60kr7(:, :, :, :) = 0.0_f64
        do i_quad_1 = 0_i64, 3_i64
          x1 = local_x1(i_quad_1)
          do i_quad_2 = 0_i64, 3_i64
            x2 = local_x2(i_quad_2)
            do i_basis_1 = 0_i64, 2_i64
              do i_basis_2 = 0_i64, 2_i64
                do j_basis_1 = 0_i64, 2_i64
                  v2_1 = global_test_basis_v2_1(i_quad_1, 0_i64, &
                        i_basis_1, i_element_1)
                  v2_1_x1 = global_test_basis_v2_1(i_quad_1, 1_i64, &
                        i_basis_1, i_element_1)
                  u2_1 = global_trial_basis_u2_1(i_quad_1, 0_i64, &
                        j_basis_1, i_element_1)
                  u2_1_x1 = global_trial_basis_u2_1(i_quad_1, 1_i64, &
                        j_basis_1, i_element_1)
                  do j_basis_2 = 0_i64, 2_i64
                    v2_2 = global_test_basis_v2_2(i_quad_2, 0_i64, &
                          i_basis_2, i_element_2)
                    v2_2_x2 = global_test_basis_v2_2(i_quad_2, 1_i64, &
                          i_basis_2, i_element_2)
                    u2_2 = global_trial_basis_u2_2(i_quad_2, 0_i64, &
                          j_basis_2, i_element_2)
                    u2_2_x2 = global_trial_basis_u2_2(i_quad_2, 1_i64, &
                          j_basis_2, i_element_2)
                    v2 = v2_1 * v2_2
                    v2_x2 = v2_1 * v2_2_x2
                    v2_x1 = v2_1_x1 * v2_2
                    u2 = u2_1 * u2_2
                    u2_x2 = u2_1 * u2_2_x2
                    u2_x1 = u2_1_x1 * u2_2
                    temp_v2_u2_0 = 2_i64 * 3.141592653589793_f64
                    temp_v2_u2_1 = temp_v2_u2_0 * x1
                    temp_v2_u2_2 = temp_v2_u2_0 * x2
                    temp_v2_u2_3 = (0.25_f64 * sin(temp_v2_u2_1) * cos( &
                          temp_v2_u2_2) + 0.25_f64 * sin(temp_v2_u2_2) &
                          * cos(temp_v2_u2_1) + 1.0_f64) ** 2_i64
                    contribution_v2_u2_wti60kr7 = 0.25_f64 * u2 * v2 / &
                          sqrt(temp_v2_u2_3)
                    l_mat_u2_v2_wti60kr7(2_i64 - i_basis_2 + j_basis_2, &
                          2_i64 - i_basis_1 + j_basis_1, i_basis_2, &
                          i_basis_1) = l_mat_u2_v2_wti60kr7(2_i64 - &
                          i_basis_2 + j_basis_2, 2_i64 - i_basis_1 + &
                          j_basis_1, i_basis_2, i_basis_1) + &
                          contribution_v2_u2_wti60kr7
                  end do
                end do
              end do
            end do
          end do
        end do
        g_mat_u2_v2_wti60kr7(:, :, pad2 + span_v2_2 - test_v2_p2:1_i64 + &
              pad2 + span_v2_2 - 1_i64, pad1 + span_v2_1 - test_v2_p1: &
              1_i64 + pad1 + span_v2_1 - 1_i64) = g_mat_u2_v2_wti60kr7( &
              :, :, pad2 + span_v2_2 - test_v2_p2:1_i64 + pad2 + &
              span_v2_2 - 1_i64, pad1 + span_v2_1 - test_v2_p1:1_i64 + &
              pad1 + span_v2_1 - 1_i64) + l_mat_u2_v2_wti60kr7(:, :, :, &
              :)
      end do
    end do
    if (allocated(l_mat_u2_v2_wti60kr7)) then
      deallocate(l_mat_u2_v2_wti60kr7)
    end if
    if (allocated(local_x2)) then
      deallocate(local_x2)
    end if
    if (allocated(local_x1)) then
      deallocate(local_x1)
    end if
    return

  end subroutine assemble_matrix_wti60kr7
  !........................................

end module dependencies_wti60kr7_6feykyplisif
