program lhapdf_example
    use iso_c_binding
    use pineappl

    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type(pineappl_grid) :: grid
    type(pineappl_lumi) :: lumi
    type(pineappl_keyval) :: key_vals

    procedure (pineappl_xfx), pointer :: xfx
    procedure (pineappl_alphas), pointer :: alphas

    integer, target :: flags(2)

    lumi = pineappl_lumi_new()
    call pineappl_lumi_add(lumi, 2, [0, 0, 1, -1, 2, -2], [1.0_dp, 1.0_dp, 1.0_dp])

    key_vals = pineappl_keyval_new()
    grid = pineappl_grid_new(lumi, 1, [2, 0, 0, 0], 2, [0.0_dp, 1.0_dp, 2.0_dp], key_vals)

    call pineappl_grid_fill_all(grid, 0.5_dp, 0.5_dp, 100.0_dp, 0, 0.5_dp, [0.5_dp, 0.5_dp, 0.5_dp])
    call pineappl_grid_fill_all(grid, 0.5_dp, 0.5_dp, 100.0_dp, 0, 1.5_dp, [1.5_dp, 1.5_dp, 1.5_dp])

    call lhapdf_initpdfset_byname(0, "nCTEQ15_1_1")
    call lhapdf_initpdfset_byname(1, "nCTEQ15FullNuc_208_82")
    
    ! calling pineappl_grid_convolve without any flags
    xfx => xfx_test1
    alphas => alphas_test1
    write(*, *) "first pineappl_grid_convolve_with_one: "
    write(*, *) pineappl_grid_convolve_with_one(grid, 2212, xfx, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp)

    ! calling pineappl_grid_convolve with two integer flags that are used in xfx_test2 and alphas_test2 to determine the set and member indices
    xfx => xfx_test2
    alphas => alphas_test2
    flags = [1, 0]
    write(*, *) "second pineappl_grid_convolve_with_one: "
    write(*, *) pineappl_grid_convolve_with_one(grid, 2212, xfx, alphas, &
        [.true., .true.], [.true., .true.], 1.0_dp, 1.0_dp, c_loc(flags(1)))
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

        call lhapdf_alphasq2(0, 0, q2, alphas_test2)
    end function

end program lhapdf_example
