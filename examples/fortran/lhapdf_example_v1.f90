program lhapdf_example
    use iso_c_binding
    use pineappl

    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type(pineappl_grid) :: grid
    type(pineappl_lumi) :: channels
    type(pineappl_kinematics) :: kinematics(3)
    type(pineappl_interp_tuples) :: interpolations(3)

    type(pineappl_xfx) :: xfx(2)
    type(pineappl_alphas) :: alphas

    integer(kind(pineappl_reweight_meth)) :: q2_reweight
    integer(kind(pineappl_reweight_meth)) :: x_reweight
    integer(kind(pineappl_map)) :: q2_mapping
    integer(kind(pineappl_map)) :: x_mapping
    integer(kind(pineappl_interp_meth)) :: interpolation_meth

    integer, target :: flags(2)

    channels = pineappl_channels_new()
    call pineappl_channels_add(channels, 3, 2, [0, 0, 1, -1, 2, -2], [1.0_dp, 1.0_dp, 1.0_dp])

    kinematics = [&
        pineappl_kinematics(pineappl_scale, 0), &
        pineappl_kinematics(pineappl_x, 0), &
        pineappl_kinematics(pineappl_x, 1) &
    ]

    q2_reweight = pineappl_no_reweight
    x_reweight = pineappl_applgrid_x
    q2_mapping = pineappl_applgrid_h0
    x_mapping = pineappl_applgrid_f2
    interpolation_meth = pineappl_lagrange
    interpolations = [ &
        pineappl_interp_tuples(1e2, 1e8, 40, 3, q2_reweight, q2_mapping, interpolation_meth), &
        pineappl_interp_tuples(2e-7, 1.0, 50, 3, x_reweight, x_mapping, interpolation_meth), &
        pineappl_interp_tuples(2e-7, 1.0, 50, 3, x_reweight, x_mapping, interpolation_meth) &
    ]

    grid = pineappl_grid_new2(pineappl_pdg, channels, 1, [2_1, 0_1, 0_1, 0_1, 0_1], 2, &
        [0.0_dp, 1.0_dp, 2.0_dp], 2, [pineappl_unpol_pdf, pineappl_unpol_pdf], [2212, 2212], kinematics, interpolations, [1, 1, 0])

    call pineappl_grid_fill_all2(grid, 0, 0.5_dp, [100.0_dp, 0.5_dp, 0.5_dp], [0.5_dp, 0.5_dp, 0.5_dp])
    call pineappl_grid_fill_all2(grid, 0, 1.5_dp, [100.0_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 1.5_dp])

    call lhapdf_initpdfset_byname(0, "nCTEQ15_1_1")
    ! call lhapdf_initpdfset_byname(0, "nCTEQ15FullNuc_208_82")
    call lhapdf_initpdfset_byname(1, "nCTEQ15FullNuc_208_82")

    ! write(*, *) "xfx_test1: ", xfx_test1(0, 0.5_dp, 100.0_dp, c_null_ptr)

    ! calling pineappl_grid_convolve without any flags
    xfx = pineappl_xfx(xfx_test1)
    alphas = pineappl_alphas(alphas_test1)
    write(*, *) "first pineappl_grid_convolve: "
    write(*, *) pineappl_grid_convolve(grid, [xfx, xfx], alphas, &
        [.true.], [.true.], [0, 1], 1, [1.0_dp, 1.0_dp, 1.0_dp])

    ! calling pineappl_grid_convolve with two integer flags that are used in xfx_test2 and alphas_test2 to determine the set and member indices
    xfx = pineappl_xfx(xfx_test2)
    alphas = pineappl_alphas(alphas_test2)
    flags = [1, 0]
    write(*, *) "second pineappl_grid_convolve: "
    write(*, *) pineappl_grid_convolve(grid, [xfx, xfx], alphas, &
        [.true.], [.true.], [0, 1], 1, [1.0_dp, 1.0_dp, 1.0_dp], c_loc(flags(1)))
contains

    ! Passing a Fortran procedure to C needs the iso_c_binding
    function xfx_test1(pdg_id, x, q2, state) bind(c)
        use iso_c_binding

        implicit none

        integer(c_int32_t), value, intent(in) :: pdg_id
        real(c_double), value, intent(in)     :: x, q2
        type(c_ptr), value, intent(in)        :: state
        real(c_double)                        :: xfx_test1

        call lhapdf_xfxq2(0, 0, pdg_id, x, q2, xfx_test1)
    end function

    function xfx_test2(pdg_id, x, q2, state) bind(c)
        use iso_c_binding

        implicit none

        integer(c_int32_t), value, intent(in) :: pdg_id
        real(c_double), value, intent(in)     :: x, q2
        type(c_ptr), value, intent(in)        :: state
        real(c_double)                        :: xfx_test2

        integer, pointer :: flags(:)

        call c_f_pointer(state, flags, [2])

        call lhapdf_xfxq2(flags(1), flags(2), pdg_id, x, q2, xfx_test2)
    end function

    function alphas_test1(q2, state) bind(c)
        use iso_c_binding

        implicit none

        real(c_double), value, intent(in) :: q2
        type(c_ptr), value, intent(in)    :: state
        real(c_double)                    :: alphas_test1

        call lhapdf_alphasq2(0, 0, q2, alphas_test1)
    end function

    function alphas_test2(q2, state) bind(c)
        use iso_c_binding

        implicit none

        real(c_double), value, intent(in) :: q2
        type(c_ptr), value, intent(in)    :: state
        real(c_double)                    :: alphas_test2

        integer, pointer :: flags(:)

        call c_f_pointer(state, flags, [2])

        call lhapdf_alphasq2(flags(1), flags(2), q2, alphas_test2)
    end function

end program lhapdf_example
