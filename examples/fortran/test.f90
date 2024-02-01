program test_pineappl
    use pineappl
    use iso_c_binding
    
    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type(pineappl_lumi) :: lumi, lumi2
    type(pineappl_grid) :: grid, grid2
    type(pineappl_keyval) :: key_vals, key_vals2

    real(dp), allocatable :: result(:), bin_limits_left(:), bin_limits_right(:), bin_normalizations(:)

    character(len=:), allocatable :: string

    procedure (pineappl_xfx), pointer    :: xfx1, xfx2
    procedure (pineappl_alphas), pointer :: alphas

    lumi = pineappl_lumi_new()
    call pineappl_lumi_add(lumi, 2, [0, 0, 1, -1], [1.0_dp, 1.0_dp])

    if (pineappl_lumi_count(lumi) /= 1) then
        write(*, *) "pineappl_lumi_count(): ", pineappl_lumi_count(lumi)
        error stop "error: pineappl_lumi_count"
    end if

    if (pineappl_lumi_combinations(lumi, 0) /= 2) then
        write(*, *) "pineappl_lumi_combinations(): ", pineappl_lumi_combinations(lumi, 0)
        error stop "error: pineappl_lumi_combinations"
    end if

    key_vals = pineappl_keyval_new()

    grid = pineappl_grid_new(lumi, 1, [2, 0, 0, 0], 2, [0.0_dp, 1.0_dp, 2.0_dp], key_vals)

    call pineappl_grid_fill_all(grid, 0.5_dp, 0.5_dp, 100.0_dp, 0, 0.5_dp, [15.0_dp, 16.0_dp])
    call pineappl_grid_fill_array(grid, [0.4_dp, 0.6_dp], [0.6_dp, 0.4_dp], [100.0_dp, 110.0_dp], &
       [0, 0], [0.5_dp, 1.5_dp], [0, 0], [20.0_dp, 21.0_dp])

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

    lumi2 = pineappl_grid_lumi(grid)

    grid2 = pineappl_grid_new(lumi, 1, [2, 0, 0, 0], 1, [2.0_dp, 3.0_dp], key_vals)

    call pineappl_grid_merge_and_delete(grid, grid2)

    call pineappl_grid_merge_bins(grid, 2, 3)

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

    call pineappl_keyval_set_bool(key_vals, "set_bool", .true.)

    call pineappl_keyval_set_double(key_vals, "set_double", 1.0_dp)

    call pineappl_keyval_set_int(key_vals, "set_int", 1)

    call pineappl_keyval_set_string(key_vals, "set_string", "set_string: success")

    if (pineappl_keyval_bool(key_vals, "set_bool") .neqv. .true.) then
        write(*, *) "pineappl_keyval_bool(): ", pineappl_keyval_bool(key_vals, "set_bool")
        error stop "error: pineappl_keyval_bool"
    end if

    if (abs(pineappl_keyval_double(key_vals, "set_double") - 1.0_dp) > 1e-10) then
        write(*, *) "pineappl_keyval_double(): ", pineappl_keyval_double(key_vals, "set_double")
        error stop "error: pineappl_keyval_double"
    end if

    if (pineappl_keyval_int(key_vals, "set_int") /= 1) then
        write(*, *) "pineappl_keyval_int(): ", pineappl_keyval_int(key_vals, "set_int")
        error stop "error: pineappl_keyval_int"
    end if

    string = pineappl_keyval_string(key_vals, "set_string")
    if (string /= "set_string: success") then
        write(*, *) "pineappl_keyval_string(): ", pineappl_keyval_string(key_vals, "set_string")
        error stop "error: pineappl_keyval_string"
    end if

    xfx1 => xfx1_test
    xfx2 => xfx2_test
    alphas => alphas_test

    result = pineappl_grid_convolute_with_one(grid, 2212, xfx1, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp)
    if (any(result > 0 .neqv. [.true., .true., .false.])) then
        write(*, *) "pineappl_grid_convolute_with_one(): ", result
        error stop "error: pineappl_grid_convolute_with_one"
    end if

    result = pineappl_grid_convolute_with_two(grid, 2212, xfx1, 2212, xfx2, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp)
    if (any(result < 0 .neqv. [.true., .true., .false.])) then
        write(*, *) "pineappl_grid_convolute_with_two(): ", result
        error stop "error: pineappl_grid_convolute_with_two"
    end if

    call pineappl_keyval_delete(key_vals)

    call pineappl_keyval_delete(key_vals2)

    call pineappl_lumi_delete(lumi)

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