module bind_c_dependencies_vga8r193_6xftl47r79t3

  use dependencies_vga8r193_6xftl47r79t3

  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T , C_F_Pointer , c_ptr
  implicit none

  contains

  !........................................
  subroutine bind_c_lo_dot_vga8r193(bound_mat00, bound_mat00_shape_1, &
        bound_mat00_shape_2, bound_mat00_shape_3, bound_mat00_shape_4, &
        bound_mat00_stride_1, bound_mat00_stride_2, &
        bound_mat00_stride_3, bound_mat00_stride_4, bound_mat01, &
        bound_mat01_shape_1, bound_mat01_shape_2, bound_mat01_shape_3, &
        bound_mat01_shape_4, bound_mat01_stride_1, bound_mat01_stride_2 &
        , bound_mat01_stride_3, bound_mat01_stride_4, bound_mat10, &
        bound_mat10_shape_1, bound_mat10_shape_2, bound_mat10_shape_3, &
        bound_mat10_shape_4, bound_mat10_stride_1, bound_mat10_stride_2 &
        , bound_mat10_stride_3, bound_mat10_stride_4, bound_mat11, &
        bound_mat11_shape_1, bound_mat11_shape_2, bound_mat11_shape_3, &
        bound_mat11_shape_4, bound_mat11_stride_1, bound_mat11_stride_2 &
        , bound_mat11_stride_3, bound_mat11_stride_4, bound_x0, &
        bound_x0_shape_1, bound_x0_shape_2, bound_x0_stride_1, &
        bound_x0_stride_2, bound_x1, bound_x1_shape_1, bound_x1_shape_2 &
        , bound_x1_stride_1, bound_x1_stride_2, bound_out0, &
        bound_out0_shape_1, bound_out0_shape_2, bound_out0_stride_1, &
        bound_out0_stride_2, bound_out1, bound_out1_shape_1, &
        bound_out1_shape_2, bound_out1_stride_1, bound_out1_stride_2, &
        s00_1, s00_2, s01_1, s01_2, s10_1, s10_2, s11_1, s11_2, n00_1, &
        n00_2, n01_1, n01_2, n10_1, n10_2, n11_1, n11_2, ne00_1, ne00_2 &
        , ne01_1, ne01_2, ne10_1, ne10_2, ne11_1, ne11_2) bind(c)

    implicit none

    type(c_ptr), value :: bound_mat00
    integer(i64), value :: bound_mat00_shape_1
    integer(i64), value :: bound_mat00_shape_2
    integer(i64), value :: bound_mat00_shape_3
    integer(i64), value :: bound_mat00_shape_4
    integer(i64), value :: bound_mat00_stride_1
    integer(i64), value :: bound_mat00_stride_2
    integer(i64), value :: bound_mat00_stride_3
    integer(i64), value :: bound_mat00_stride_4
    type(c_ptr), value :: bound_mat01
    integer(i64), value :: bound_mat01_shape_1
    integer(i64), value :: bound_mat01_shape_2
    integer(i64), value :: bound_mat01_shape_3
    integer(i64), value :: bound_mat01_shape_4
    integer(i64), value :: bound_mat01_stride_1
    integer(i64), value :: bound_mat01_stride_2
    integer(i64), value :: bound_mat01_stride_3
    integer(i64), value :: bound_mat01_stride_4
    type(c_ptr), value :: bound_mat10
    integer(i64), value :: bound_mat10_shape_1
    integer(i64), value :: bound_mat10_shape_2
    integer(i64), value :: bound_mat10_shape_3
    integer(i64), value :: bound_mat10_shape_4
    integer(i64), value :: bound_mat10_stride_1
    integer(i64), value :: bound_mat10_stride_2
    integer(i64), value :: bound_mat10_stride_3
    integer(i64), value :: bound_mat10_stride_4
    type(c_ptr), value :: bound_mat11
    integer(i64), value :: bound_mat11_shape_1
    integer(i64), value :: bound_mat11_shape_2
    integer(i64), value :: bound_mat11_shape_3
    integer(i64), value :: bound_mat11_shape_4
    integer(i64), value :: bound_mat11_stride_1
    integer(i64), value :: bound_mat11_stride_2
    integer(i64), value :: bound_mat11_stride_3
    integer(i64), value :: bound_mat11_stride_4
    type(c_ptr), value :: bound_x0
    integer(i64), value :: bound_x0_shape_1
    integer(i64), value :: bound_x0_shape_2
    integer(i64), value :: bound_x0_stride_1
    integer(i64), value :: bound_x0_stride_2
    type(c_ptr), value :: bound_x1
    integer(i64), value :: bound_x1_shape_1
    integer(i64), value :: bound_x1_shape_2
    integer(i64), value :: bound_x1_stride_1
    integer(i64), value :: bound_x1_stride_2
    type(c_ptr), value :: bound_out0
    integer(i64), value :: bound_out0_shape_1
    integer(i64), value :: bound_out0_shape_2
    integer(i64), value :: bound_out0_stride_1
    integer(i64), value :: bound_out0_stride_2
    type(c_ptr), value :: bound_out1
    integer(i64), value :: bound_out1_shape_1
    integer(i64), value :: bound_out1_shape_2
    integer(i64), value :: bound_out1_stride_1
    integer(i64), value :: bound_out1_stride_2
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
    real(f64), pointer :: mat00(:,:,:,:)
    real(f64), pointer :: mat01(:,:,:,:)
    real(f64), pointer :: mat10(:,:,:,:)
    real(f64), pointer :: mat11(:,:,:,:)
    real(f64), pointer :: x0(:,:)
    real(f64), pointer :: x1(:,:)
    real(f64), pointer :: out0(:,:)
    real(f64), pointer :: out1(:,:)

    call C_F_Pointer(bound_mat00, mat00, [bound_mat00_shape_4 * &
          bound_mat00_stride_4,bound_mat00_shape_3 * &
          bound_mat00_stride_3,bound_mat00_shape_2 * &
          bound_mat00_stride_2,bound_mat00_shape_1 * &
          bound_mat00_stride_1])
    call C_F_Pointer(bound_mat01, mat01, [bound_mat01_shape_4 * &
          bound_mat01_stride_4,bound_mat01_shape_3 * &
          bound_mat01_stride_3,bound_mat01_shape_2 * &
          bound_mat01_stride_2,bound_mat01_shape_1 * &
          bound_mat01_stride_1])
    call C_F_Pointer(bound_mat10, mat10, [bound_mat10_shape_4 * &
          bound_mat10_stride_4,bound_mat10_shape_3 * &
          bound_mat10_stride_3,bound_mat10_shape_2 * &
          bound_mat10_stride_2,bound_mat10_shape_1 * &
          bound_mat10_stride_1])
    call C_F_Pointer(bound_mat11, mat11, [bound_mat11_shape_4 * &
          bound_mat11_stride_4,bound_mat11_shape_3 * &
          bound_mat11_stride_3,bound_mat11_shape_2 * &
          bound_mat11_stride_2,bound_mat11_shape_1 * &
          bound_mat11_stride_1])
    call C_F_Pointer(bound_x0, x0, [bound_x0_shape_2 * bound_x0_stride_2 &
          ,bound_x0_shape_1 * bound_x0_stride_1])
    call C_F_Pointer(bound_x1, x1, [bound_x1_shape_2 * bound_x1_stride_2 &
          ,bound_x1_shape_1 * bound_x1_stride_1])
    call C_F_Pointer(bound_out0, out0, [bound_out0_shape_2 * &
          bound_out0_stride_2,bound_out0_shape_1 * bound_out0_stride_1] &
          )
    call C_F_Pointer(bound_out1, out1, [bound_out1_shape_2 * &
          bound_out1_stride_2,bound_out1_shape_1 * bound_out1_stride_1] &
          )
    call lo_dot_vga8r193(mat00 = mat00(1_i64::bound_mat00_stride_4, &
          1_i64::bound_mat00_stride_3, 1_i64::bound_mat00_stride_2, &
          1_i64::bound_mat00_stride_1), mat01 = mat01(1_i64:: &
          bound_mat01_stride_4, 1_i64::bound_mat01_stride_3, 1_i64:: &
          bound_mat01_stride_2, 1_i64::bound_mat01_stride_1), mat10 = &
          mat10(1_i64::bound_mat10_stride_4, 1_i64:: &
          bound_mat10_stride_3, 1_i64::bound_mat10_stride_2, 1_i64:: &
          bound_mat10_stride_1), mat11 = mat11(1_i64:: &
          bound_mat11_stride_4, 1_i64::bound_mat11_stride_3, 1_i64:: &
          bound_mat11_stride_2, 1_i64::bound_mat11_stride_1), x0 = x0( &
          1_i64::bound_x0_stride_2, 1_i64::bound_x0_stride_1), x1 = x1( &
          1_i64::bound_x1_stride_2, 1_i64::bound_x1_stride_1), out0 = &
          out0(1_i64::bound_out0_stride_2, 1_i64::bound_out0_stride_1), &
          out1 = out1(1_i64::bound_out1_stride_2, 1_i64:: &
          bound_out1_stride_1), s00_1 = s00_1, s00_2 = s00_2, s01_1 = &
          s01_1, s01_2 = s01_2, s10_1 = s10_1, s10_2 = s10_2, s11_1 = &
          s11_1, s11_2 = s11_2, n00_1 = n00_1, n00_2 = n00_2, n01_1 = &
          n01_1, n01_2 = n01_2, n10_1 = n10_1, n10_2 = n10_2, n11_1 = &
          n11_1, n11_2 = n11_2, ne00_1 = ne00_1, ne00_2 = ne00_2, &
          ne01_1 = ne01_1, ne01_2 = ne01_2, ne10_1 = ne10_1, ne10_2 = &
          ne10_2, ne11_1 = ne11_1, ne11_2 = ne11_2)

  end subroutine bind_c_lo_dot_vga8r193
  !........................................

end module bind_c_dependencies_vga8r193_6xftl47r79t3
