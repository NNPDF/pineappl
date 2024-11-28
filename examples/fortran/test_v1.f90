program test_pineappl
    use pineappl
    use iso_c_binding
    
    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type(pineappl_lumi) :: channels, channels2
    type(pineappl_grid) :: grid, grid2
    type(pineappl_kinematics) :: kinematics(3)
    type(pineappl_interp_tuples) :: interpolations(3)

    real(dp), allocatable :: result(:), bin_limits_left(:), bin_limits_right(:), bin_normalizations(:)

    integer(kind(pineappl_reweight_meth)) :: q2_reweight
    integer(kind(pineappl_reweight_meth)) :: x_reweight
    integer(kind(pineappl_map)) :: q2_mapping
    integer(kind(pineappl_map)) :: x_mapping
    integer(kind(pineappl_interp_meth)) :: interpolation_meth

    type (pineappl_xfx)    :: xfx1, xfx2
    type (pineappl_alphas) :: alphas

    channels = pineappl_channels_new()
    call pineappl_channels_add(channels, 3, 2, [0, 0, 1, -1, 2, -2], [1.0_dp, 1.0_dp, 1.0_dp])
    
    if (pineappl_lumi_count(channels) /= 1) then
        write(*, *) "pineappl_lumi_count(): ", pineappl_lumi_count(channels)
        error stop "error: pineappl_lumi_count"
    end if

    if (pineappl_lumi_combinations(channels, 0) /= 3) then
        write(*, *) "pineappl_lumi_combinations(): ", pineappl_lumi_combinations(channels, 0)
        error stop "error: pineappl_lumi_combinations"
    end if

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
    grid = pineappl_grid_new2(pineappl_pdg, channels, 1, [2_1, 0_1, 0_1, 0_1], 2, [0.0_dp, 1.0_dp, 2.0_dp], &
        2, [pineappl_unpol_pdf, pineappl_unpol_pdf], [2212, 2212], kinematics, interpolations, [1, 1, 0])

    if (pineappl_grid_order_count(grid) /= 1) then
        write(*, *) "pineappl_grid_order_count(): ", pineappl_grid_order_count(grid)
        error stop "error: pineappl_grid_order_count"
    end if

    if (any(pineappl_grid_order_params(grid) /= [2, 0, 0, 0])) then
        write(*, *) "pineappl_grid_order_params(): ", pineappl_grid_order_params(grid)
        error stop "error: pineappl_grid_order_params"
    end if

    call pineappl_grid_fill2(grid, 0, 0.5_dp, 0, [100.0_dp, 0.5_dp, 0.5_dp], 14.0_dp)
    call pineappl_grid_fill_all2(grid, 0, 0.5_dp, [100.0_dp, 0.5_dp, 0.5_dp], [15.0_dp, 16.0_dp])
    call pineappl_grid_fill_array2(grid, [0, 0], [1.5_dp, 1.5_dp], &
        [100.0_dp, 0.4_dp, 0.6_dp, 110.0_dp, 0.6_dp, 0.4_dp], [0, 0], [20.0_dp, 21.0_dp])

    if (pineappl_grid_bin_count(grid) /= 2) then
        write(*, *) "pineappl_grid_bin_count(): ", pineappl_grid_bin_count(grid)
        error stop "error: pineappl_grid_bin_count"
    end if

    if (pineappl_grid_bin_dimensions(grid) /= 1) then
        write(*, *) "pineappl_grid_bin_dimensions(): ", pineappl_grid_bin_dimensions(grid)
        error stop "error: pineappl_grid_bin_dimensions"
    end if

    bin_limits_left = pineappl_grid_bin_limits_left(grid, 0)
    if (any(abs(bin_limits_left - [0.0_dp, 1.0_dp]) > 1e-10)) then
        write(*, *) "pineappl_grid_bin_limits_left(): ", abs(bin_limits_left - [0.0_dp, 1.0_dp]) < 1e-6
        error stop "error: pineappl_grid_bin_limits_left"
    end if

    bin_limits_right = pineappl_grid_bin_limits_right(grid, 0)
    if (any(abs(bin_limits_right - [1.0_dp, 2.0_dp]) > 1e-10)) then
        write(*, *) "pineappl_grid_bin_limits_right(): ", bin_limits_right
        error stop "error: pineappl_grid_bin_limits_right"
    end if

    bin_normalizations = pineappl_grid_bin_normalizations(grid)
    if (any(abs(bin_normalizations - [1.0_dp, 1.0_dp]) > 1e-10)) then
        write(*, *) "pineappl_grid_bin_normalizations(): ", bin_normalizations
        error stop "error: pineappl_grid_bin_normalizations"
    end if

    grid2 = pineappl_grid_clone(grid)

    call pineappl_grid_delete(grid2)

    channels2 = pineappl_grid_lumi(grid)

    if (pineappl_lumi_count(channels2) /= 1) then
        write(*, *) "pineappl_lumi_count(): ", pineappl_lumi_count(channels2)
        error stop "error: pineappl_lumi_count"
    end if

    if (pineappl_lumi_combinations(channels2, 0) /= 3) then
        write(*, *) "pineappl_lumi_combinations(): ", pineappl_lumi_combinations(channels2, 0)
        error stop "error: pineappl_lumi_combinations"
    end if

    grid2 = pineappl_grid_new2(pineappl_pdg, channels, 1, [2_1, 0_1, 0_1, 0_1], 1, [2.0_dp, 3.0_dp], &
    2, [pineappl_unpol_pdf, pineappl_unpol_pdf], [2212, 2212], kinematics, interpolations, [1, 1, 0])

    call pineappl_grid_merge_and_delete(grid, grid2)

    if (pineappl_grid_order_count(grid) /= 1) then
        write(*, *) "pineappl_grid_order_count(): ", pineappl_grid_order_count(grid)
        error stop "error: pineappl_grid_order_count"
    end if

    call pineappl_grid_merge_bins(grid, 2, 3)

    if (pineappl_grid_order_count(grid) /= 1) then
        write(*, *) "pineappl_grid_order_count(): ", pineappl_grid_order_count(grid)
        error stop "error: pineappl_grid_order_count"
    end if

    call pineappl_grid_optimize_using(grid, int(b'11111'))

    if (pineappl_grid_order_count(grid) /= 1) then
        write(*, *) "pineappl_grid_order_count(): ", pineappl_grid_order_count(grid)
        error stop "error: pineappl_grid_order_count"
    end if

    if (any(pineappl_grid_order_params(grid) /= [2, 0, 0, 0])) then
        write(*, *) "pineappl_grid_order_params(): ", pineappl_grid_order_params(grid)
        error stop "error: pineappl_grid_order_params"
    end if

    call pineappl_grid_scale(grid, 0.5_dp)

    call pineappl_grid_scale_by_bin(grid, [2.0_dp, 2.0_dp])

    call pineappl_grid_scale_by_order(grid, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)

    call pineappl_grid_set_key_value(grid, "set_key_value", "set_key_value: success")

    ! at this point we have the bins [0, 1, 3]
    call pineappl_grid_set_remapper(grid, 2, [1.0_dp, 1.0_dp], [0.0_dp, 1.0_dp, 10.0_dp, 11.0_dp, 1.0_dp, 3.0_dp, 11.0_dp, 13.0_dp])

    call pineappl_grid_split_lumi(grid)

    xfx1 = pineappl_xfx(xfx1_test)
    xfx2 = pineappl_xfx(xfx2_test)
    alphas = pineappl_alphas(alphas_test)

    result = pineappl_grid_convolve_with_one(grid, 2212, xfx1, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp)
    if (any(result > 0 .neqv. [.true., .true., .false.])) then
        write(*, *) "pineappl_grid_convolve_with_one(): ", result
        error stop "error: pineappl_grid_convolve_with_one"
    end if

    result = pineappl_grid_convolve_with_two(grid, 2212, xfx1, 2212, xfx2, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp)
    if (any(result < 0 .neqv. [.true., .true., .false.])) then
        write(*, *) "pineappl_grid_convolve_with_two(): ", result
        error stop "error: pineappl_grid_convolve_with_two"
    end if

    result = pineappl_grid_convolve(grid, [xfx1, xfx2], alphas, [.true., .true.], [.true., .true.], &
        [0, 1, 2], 1, [1.0_dp, 1.0_dp])
    if (any(result < 0 .neqv. [.true., .true., .false.])) then
        write(*, *) "pineappl_grid_convolve_with_two(): ", result
        error stop "error: pineappl_grid_convolve_with_two"
    end if

    call pineappl_lumi_delete(channels)

    call pineappl_grid_delete(grid)

contains

    function xfx1_test(pdg_id, x, q2, state) bind(c)
        use iso_c_binding

        implicit none
        
        integer(c_int32_t), value, intent(in) :: pdg_id
        real(c_double), value, intent(in)     :: x, q2
        type(c_ptr), value, intent(in)        :: state
        real(c_double)                        :: xfx1_test

        xfx1_test = x
    end function

    function xfx2_test(pdg_id, x, q2, state) bind(c)
        use iso_c_binding

        implicit none
        
        integer(c_int32_t), value, intent(in) :: pdg_id
        real(c_double), value, intent(in)     :: x, q2
        type(c_ptr), value, intent(in)        :: state
        real(c_double)                        :: xfx2_test

        xfx2_test = -x
    end function

    function alphas_test(q2, state) bind(c)
        use iso_c_binding

        implicit none
        
        real(c_double), value, intent(in) :: q2
        type(c_ptr), value, intent(in)    :: state
        real(c_double)                    :: alphas_test

        alphas_test = q2
    end function

end program test_pineappl
