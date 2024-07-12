module dependencies_vga8r193_6xftl47r79t3


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine lo_dot_vga8r193(mat00, mat01, mat10, mat11, x0, x1, out0, &
        out1, s00_1, s00_2, s01_1, s01_2, s10_1, s10_2, s11_1, s11_2, &
        n00_1, n00_2, n01_1, n01_2, n10_1, n10_2, n11_1, n11_2, ne00_1, &
        ne00_2, ne01_1, ne01_2, ne10_1, ne10_2, ne11_1, ne11_2)

    implicit none

    real(f64), intent(in) :: mat00(0:,0:,0:,0:)
    real(f64), intent(in) :: mat01(0:,0:,0:,0:)
    real(f64), intent(in) :: mat10(0:,0:,0:,0:)
    real(f64), intent(in) :: mat11(0:,0:,0:,0:)
    real(f64), intent(in) :: x0(0:,0:)
    real(f64), intent(in) :: x1(0:,0:)
    real(f64), intent(inout) :: out0(0:,0:)
    real(f64), intent(inout) :: out1(0:,0:)
    integer(i64), value :: s00_1
    integer(i64), value :: s00_2
    integer(i64), value :: s01_1
    integer(i64), value :: s01_2
    integer(i64), value :: s10_1
    integer(i64), value :: s10_2
    integer(i64), value :: s11_1
    integer(i64), value :: s11_2
    integer(i64), value :: n00_1
    integer(i64), value :: n00_2
    integer(i64), value :: n01_1
    integer(i64), value :: n01_2
    integer(i64), value :: n10_1
    integer(i64), value :: n10_2
    integer(i64), value :: n11_1
    integer(i64), value :: n11_2
    integer(i64), value :: ne00_1
    integer(i64), value :: ne00_2
    integer(i64), value :: ne01_1
    integer(i64), value :: ne01_2
    integer(i64), value :: ne10_1
    integer(i64), value :: ne10_2
    integer(i64), value :: ne11_1
    integer(i64), value :: ne11_2
    integer(i64) :: i1
    integer(i64) :: i2
    real(f64) :: v00
    integer(i64) :: k1
    integer(i64) :: k2
    real(f64) :: v11
    real(f64) :: v01
    real(f64) :: v10

    do i1 = 0_i64, n00_1 - 1_i64
      do i2 = 0_i64, n00_2 - 1_i64
        v00 = 0.0_f64
        do k1 = 0_i64, 4_i64
          do k2 = 0_i64, 6_i64
            v00 = v00 + mat00(k2, k1, 3_i64 + i2, 3_i64 + i1) * x0(i2 + &
                  k2, 1_i64 + i1 + k1)
          end do
        end do
        out0(3_i64 + i2, 3_i64 + i1) = v00
      end do
    end do
    do i1 = 0_i64, ne00_1 - 1_i64
      do i2 = 0_i64, n00_2 - 1_i64
        v00 = 0.0_f64
        do k1 = 0_i64, 4_i64 - i1 - 1_i64
          do k2 = 0_i64, 6_i64
            v00 = v00 + x0(i2 + k2, 1_i64 + i1 + k1 + n00_1) * mat00(k2, &
                  k1, 3_i64 + i2, 3_i64 + i1 + n00_1)
          end do
        end do
        out0(3_i64 + i2, 3_i64 + i1 + n00_1) = v00
      end do
    end do
    do i1 = 0_i64, n00_1 + ne00_1 - 1_i64
      do i2 = 0_i64, ne00_2 - 1_i64
        v00 = 0.0_f64
        do k1 = 0_i64, 5_i64 - maxval([0_i64, i1 - n00_1 + 1_i64]) - &
              1_i64
          do k2 = 0_i64, 6_i64 - i2 - 1_i64
            v00 = v00 + x0(i2 + k2 + n00_2, 1_i64 + i1 + k1) * mat00(k2, &
                  k1, 3_i64 + i2 + n00_2, 3_i64 + i1)
          end do
        end do
        out0(3_i64 + i2 + n00_2, 3_i64 + i1) = v00
      end do
    end do
    do i1 = 0_i64, n11_1 - 1_i64
      do i2 = 0_i64, n11_2 - 1_i64
        v11 = 0.0_f64
        do k1 = 0_i64, 6_i64
          do k2 = 0_i64, 4_i64
            v11 = v11 + mat11(k2, k1, 3_i64 + i2, 3_i64 + i1) * x1(1_i64 &
                  + i2 + k2, i1 + k1)
          end do
        end do
        out1(3_i64 + i2, 3_i64 + i1) = v11
      end do
    end do
    do i1 = 0_i64, ne11_1 - 1_i64
      do i2 = 0_i64, n11_2 - 1_i64
        v11 = 0.0_f64
        do k1 = 0_i64, 6_i64 - i1 - 1_i64
          do k2 = 0_i64, 4_i64
            v11 = v11 + x1(1_i64 + i2 + k2, i1 + k1 + n11_1) * mat11(k2, &
                  k1, 3_i64 + i2, 3_i64 + i1 + n11_1)
          end do
        end do
        out1(3_i64 + i2, 3_i64 + i1 + n11_1) = v11
      end do
    end do
    do i1 = 0_i64, n11_1 + ne11_1 - 1_i64
      do i2 = 0_i64, ne11_2 - 1_i64
        v11 = 0.0_f64
        do k1 = 0_i64, 7_i64 - maxval([0_i64, i1 - n11_1 + 1_i64]) - &
              1_i64
          do k2 = 0_i64, 4_i64 - i2 - 1_i64
            v11 = v11 + x1(1_i64 + i2 + k2 + n11_2, i1 + k1) * mat11(k2, &
                  k1, 3_i64 + i2 + n11_2, 3_i64 + i1)
          end do
        end do
        out1(3_i64 + i2 + n11_2, 3_i64 + i1) = v11
      end do
    end do
    do i1 = 0_i64, n01_1 - 1_i64
      do i2 = 0_i64, n01_2 - 1_i64
        v01 = 0.0_f64
        do k1 = 0_i64, 6_i64
          do k2 = 0_i64, 6_i64
            v01 = v01 + mat01(k2, k1, 3_i64 + i2, 3_i64 + i1) * x1(i2 + &
                  k2, i1 + k1)
          end do
        end do
        out0(3_i64 + i2, 3_i64 + i1) = out0(3_i64 + i2, 3_i64 + i1) + &
              v01
      end do
    end do
    do i1 = 0_i64, ne01_1 - 1_i64
      do i2 = 0_i64, n01_2 - 1_i64
        v01 = 0.0_f64
        do k1 = 0_i64, 6_i64 - i1 - 1_i64
          do k2 = 0_i64, 6_i64
            v01 = v01 + x1(i2 + k2, i1 + k1 + n01_1) * mat01(k2, k1, &
                  3_i64 + i2, 3_i64 + i1 + n01_1)
          end do
        end do
        out0(3_i64 + i2, 3_i64 + i1 + n01_1) = out0(3_i64 + i2, 3_i64 + &
              i1 + n01_1) + v01
      end do
    end do
    do i1 = 0_i64, n01_1 + ne01_1 - 1_i64
      do i2 = 0_i64, ne01_2 - 1_i64
        v01 = 0.0_f64
        do k1 = 0_i64, 7_i64 - maxval([0_i64, i1 - n01_1 + 1_i64]) - &
              1_i64
          do k2 = 0_i64, 6_i64 - i2 - 1_i64
            v01 = v01 + x1(i2 + k2 + n01_2, i1 + k1) * mat01(k2, k1, &
                  3_i64 + i2 + n01_2, 3_i64 + i1)
          end do
        end do
        out0(3_i64 + i2 + n01_2, 3_i64 + i1) = out0(3_i64 + i2 + n01_2, &
              3_i64 + i1) + v01
      end do
    end do
    do i1 = 0_i64, n10_1 - 1_i64
      do i2 = 0_i64, n10_2 - 1_i64
        v10 = 0.0_f64
        do k1 = 0_i64, 6_i64
          do k2 = 0_i64, 6_i64
            v10 = v10 + mat10(k2, k1, 3_i64 + i2, 3_i64 + i1) * x0(i2 + &
                  k2, i1 + k1)
          end do
        end do
        out1(3_i64 + i2, 3_i64 + i1) = out1(3_i64 + i2, 3_i64 + i1) + &
              v10
      end do
    end do
    do i1 = 0_i64, ne10_1 - 1_i64
      do i2 = 0_i64, n10_2 - 1_i64
        v10 = 0.0_f64
        do k1 = 0_i64, 6_i64 - i1 - 1_i64
          do k2 = 0_i64, 6_i64
            v10 = v10 + x0(i2 + k2, i1 + k1 + n10_1) * mat10(k2, k1, &
                  3_i64 + i2, 3_i64 + i1 + n10_1)
          end do
        end do
        out1(3_i64 + i2, 3_i64 + i1 + n10_1) = out1(3_i64 + i2, 3_i64 + &
              i1 + n10_1) + v10
      end do
    end do
    do i1 = 0_i64, n10_1 + ne10_1 - 1_i64
      do i2 = 0_i64, ne10_2 - 1_i64
        v10 = 0.0_f64
        do k1 = 0_i64, 7_i64 - maxval([0_i64, i1 - n10_1 + 1_i64]) - &
              1_i64
          do k2 = 0_i64, 6_i64 - i2 - 1_i64
            v10 = v10 + x0(i2 + k2 + n10_2, i1 + k1) * mat10(k2, k1, &
                  3_i64 + i2 + n10_2, 3_i64 + i1)
          end do
        end do
        out1(3_i64 + i2 + n10_2, 3_i64 + i1) = out1(3_i64 + i2 + n10_2, &
              3_i64 + i1) + v10
      end do
    end do
    return

  end subroutine lo_dot_vga8r193
  !........................................

end module dependencies_vga8r193_6xftl47r79t3
