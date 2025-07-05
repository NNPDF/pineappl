module pineappl
    use iso_c_binding
    use iso_fortran_env

    implicit none

    integer, parameter, private :: dp = kind(0.0d0)

    type pineappl_grid
        type (c_ptr) :: ptr = c_null_ptr
    end type

    type pineappl_keyval
        type (c_ptr) :: ptr = c_null_ptr
    end type

    type pineappl_lumi
        type (c_ptr) :: ptr = c_null_ptr
    end type

    type pineappl_channels
        type (c_ptr) :: ptr = c_null_ptr
    end type


    ! As a workaround for typing Fortran enums, we define the name of the enum as the last enum value. This way, variables can be declared as, e.g. for pineappl_conv_type, integer(kind(pineappl_conv_type)). The compiler doesn't check that a value is from the right enum, but it clarifies the code for the user.

    enum, bind(c) ! :: pineappl_conv_type
        enumerator :: pineappl_conv_type_unpol_pdf
        enumerator :: pineappl_conv_type_pol_pdf
        enumerator :: pineappl_conv_type_unpol_ff
        enumerator :: pineappl_conv_type_pol_ff

        enumerator :: pineappl_conv_type
    end enum

    enum, bind(c) ! :: pineappl_interp_meth
        enumerator :: pineappl_interp_meth_lagrange

        enumerator :: pineappl_interp_meth
    end enum

    enum, bind(c) ! :: pineappl_map
        enumerator :: pineappl_map_applgrid_f2
        enumerator :: pineappl_map_applgrid_h0

        enumerator :: pineappl_map
    end enum

    enum, bind(c) ! :: pineappl_pid_basis
        enumerator :: pineappl_pid_basis_pdg
        enumerator :: pineappl_pid_basis_evol

        enumerator :: pineappl_pid_basis
    end enum

    enum, bind(c) ! :: pineappl_reweight_meth
        enumerator :: pineappl_reweight_meth_applgrid_x
        enumerator :: pineappl_reweight_meth_no_reweight

        enumerator :: pineappl_reweight_meth
    end enum

    enum, bind(c) ! :: pineappl_kinematics_tag
        enumerator :: pineappl_kinematics_tag_scale
        enumerator :: pineappl_kinematics_tag_x

        enumerator :: pineappl_kinematics_tag
    end enum

    enum, bind(c) ! :: pineappl_scale_func_form_tag
        enumerator :: pineappl_scale_func_form_tag_no_scale
        enumerator :: pineappl_scale_func_form_tag_scale
        enumerator :: pineappl_scale_func_form_tag_quadratic_sum
        enumerator :: pineappl_scale_func_form_tag_quadratic_mean
        enumerator :: pineappl_scale_func_form_tag_quadratic_sum_over4
        enumerator :: pineappl_scale_func_form_tag_linear_mean
        enumerator :: pineappl_scale_func_form_tag_linear_sum
        enumerator :: pineappl_scale_func_form_tag_scale_max
        enumerator :: pineappl_scale_func_form_tag_scale_min
        enumerator :: pineappl_scale_func_form_tag_prod
        enumerator :: pineappl_scale_func_form_tag_s2plus_s1half
        enumerator :: pineappl_scale_func_form_tag_pow4_sum
        enumerator :: pineappl_scale_func_form_tag_wgt_avg
        enumerator :: pineappl_scale_func_form_tag_s2plus_s1fourth
        enumerator :: pineappl_scale_func_form_tag_exp_prod2

        enumerator :: pineappl_scale_func_form_tag
    end enum

    ! The Kinematics struct is a tuple-like struct in the PineAPPL Rust code, which is realized as a C union. Fortran does not support unions, but fortunately the union is only for storing ints, so we just use an integer variable for `index`
    type, bind(c) :: pineappl_kinematics
        integer(c_int) :: tag
        integer(c_size_t) :: index
    end type

    ! Implement the ScaleFuncForm struct which is also a tuple-like struct ine PineAPPL Rust code. The `pineappl_scale_func_form_body` objects have to defined with two fields - if not required, the value(s) will be ignored.
    type, bind(c) :: pineappl_scale_func_form_body
        integer(c_size_t) :: index_0 ! index_0 maps to C union field _0
        integer(c_size_t) :: index_1 ! index_1 maps to C union field _1
    end type

    type, bind(c) :: pineappl_scale_func_form
        integer(c_int) :: tag
        type(pineappl_scale_func_form_body) :: body
    end type

    type, bind(c) :: pineappl_conv
        integer(c_int) :: conv_type
        integer(c_int32_t) :: pid
    end type

    type, bind(c) :: pineappl_interp
        real(c_double) :: min
        real(c_double) :: max
        integer(c_size_t) :: nodes
        integer(c_size_t) :: order
        integer(c_int) :: reweighting_method
        integer(c_int) :: mapping
        integer(c_int) :: interpolation_method
    end type

    type :: pineappl_xfx
        procedure (pineappl_xfx_proc), pointer, nopass :: proc
    end type

    type :: pineappl_alphas
        procedure (pineappl_alphas_proc), pointer, nopass :: proc
    end type

    abstract interface
        function pineappl_xfx_proc(pdg_id, x, q2, state) bind(c)
            use iso_c_binding

            implicit none

            integer(c_int32_t), value, intent(in) :: pdg_id
            real(c_double), value, intent(in)     :: x, q2
            type (c_ptr), value, intent(in)       :: state
            real(c_double)                        :: pineappl_xfx_proc
        end function

        function pineappl_alphas_proc(q2, state) bind(c)
            use iso_c_binding

            implicit none

            real(c_double), value, intent(in) :: q2
            type (c_ptr), value, intent(in)   :: state
            real(c_double)                    :: pineappl_alphas_proc
        end function
    end interface

    interface
        function strlen(s) bind(c, name="strlen")
            use iso_c_binding

            implicit none

            type (c_ptr), value :: s
            integer (c_size_t)  :: strlen
        end function strlen

        subroutine channels_add(channels, combinations, pdg_id_combinations, factors) &
            bind(c, name = 'pineappl_channels_add')

            use iso_c_binding
            type (c_ptr), value       :: channels
            integer (c_size_t), value :: combinations
            integer (c_int32_t)       :: pdg_id_combinations(*)
            real (c_double)           :: factors(*)
        end subroutine

        type (c_ptr) function channels_new(convolutions) bind(c, name = 'pineappl_channels_new')
            use iso_c_binding
            integer (c_int32_t), value :: convolutions
        end function

        integer (c_size_t) function grid_bin_count(grid) bind(c, name = 'pineappl_grid_bin_count')
            use iso_c_binding
            type (c_ptr), value :: grid
        end function

        integer (c_size_t) function grid_bin_dimensions(grid) bind(c, name = 'pineappl_grid_bin_dimensions')
            use iso_c_binding
            type (c_ptr), value :: grid
        end function

        subroutine grid_bin_limits_left(grid, dimension, left) bind(c, name = 'pineappl_grid_bin_limits_left')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: dimension
            real (c_double)           :: left(*)
        end subroutine

        subroutine grid_bin_limits_right(grid, dimension, right) bind(c, name = 'pineappl_grid_bin_limits_right')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: dimension
            real (c_double)           :: right(*)
        end subroutine

        subroutine grid_bin_normalizations(grid, bin_normalizations) bind(c, name = 'pineappl_grid_bin_normalizations')
            use iso_c_binding
            type (c_ptr), value :: grid
            real (c_double)     :: bin_normalizations(*)
        end subroutine

        type (c_ptr) function grid_clone(grid) bind(c, name = 'pineappl_grid_clone')
            use iso_c_binding
            type (c_ptr), value :: grid
        end function

        subroutine grid_convolve(grid, xfx, alphas, pdfs_state, alphas_state, order_mask, channel_mask, &
            bin_indices, nb_scales, mu_scales, results) &
            bind(c, name = 'pineappl_grid_convolve')

            use iso_c_binding
            type (c_ptr), value        :: grid
            type (c_funptr), value     :: xfx
            type (c_funptr), value     :: alphas
            type(c_ptr)                :: pdfs_state(*)
            type (c_ptr), value        :: alphas_state
            logical (c_bool)           :: order_mask(*), channel_mask(*)
            integer (c_size_t)         :: bin_indices(*)
            integer (c_size_t), value  :: nb_scales
            real (c_double)            :: mu_scales(*), results(*)
        end subroutine

        subroutine grid_dedup_channels(grid, ulps) bind(c, name = 'pineappl_grid_dedup_channels')
            use iso_c_binding
            type (c_ptr), value        :: grid
            integer (c_int64_t), value :: ulps
        end subroutine

        subroutine grid_delete(grid) bind(c, name = 'pineappl_grid_delete')
            use iso_c_binding
            type (c_ptr), value :: grid
        end subroutine

        subroutine grid_fill2(grid, order, observable, channel, ntuple, weight) bind(c, name = 'pineappl_grid_fill2')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: order, channel
            real (c_double), value    :: observable, weight
            real (c_double)           :: ntuple(*)
        end subroutine

        subroutine grid_fill_all2(grid, order, observable, ntuple, weights) bind(c, name = 'pineappl_grid_fill_all2')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: order
            real (c_double), value    :: observable
            real (c_double)           :: ntuple(*)
            real (c_double)           :: weights(*)
        end subroutine

        subroutine grid_fill_array2(grid, orders, observables, ntuples, lumis, weights, size) &
            bind(c, name = 'pineappl_grid_fill_array2')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t)        :: orders(*), lumis(*)
            real (c_double)           :: observables(*), ntuples(*), weights(*)
            integer (c_size_t), value :: size
        end subroutine

        function grid_key_value(grid, key) bind(c, name = 'pineappl_grid_key_value')
            use iso_c_binding
            type (c_ptr), value :: grid
            character (c_char)  :: key(*)
            type (c_ptr)        :: grid_key_value
        end function

        function grid_channels(grid) bind(c, name = 'pineappl_grid_channels')
            use iso_c_binding
            type (c_ptr), value :: grid
            type (c_ptr)        :: grid_channels
        end function

        subroutine grid_merge_and_delete(grid, other) bind(c, name = 'pineappl_grid_merge_and_delete')
            use iso_c_binding
            type (c_ptr), value :: grid, other
        end subroutine

        subroutine grid_merge_bins(grid, from, to) bind(c, name = 'pineappl_grid_merge_bins')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: from, to
        end subroutine

        type (c_ptr) function grid_new2(bins, bin_limits, orders, order_params, channels, &
            pid_basis, convolutions, interpolations, interp_info, &
            kinematics, mu_scales) bind(c, name = 'pineappl_grid_new2')
            use iso_c_binding
            import ! so we can use pineappl_kinematics and pineappl_interp

            integer (c_int32_t), value :: pid_basis
            type (c_ptr), value :: channels
            integer (c_size_t), value :: orders, bins, interpolations
            integer (c_int8_t) :: order_params(*)
            real (c_double) :: bin_limits(*)
            type (pineappl_conv) :: convolutions(*)
            type (pineappl_kinematics) :: kinematics(*)
            type (pineappl_interp) :: interp_info(*)
            type (pineappl_scale_func_form) :: mu_scales(*)
        end function

        subroutine grid_optimize(grid) bind(c, name = 'pineappl_grid_optimize')
            use iso_c_binding
            type (c_ptr), value :: grid
        end subroutine

        subroutine grid_optimize_using(grid, flags) bind(c, name = 'pineappl_grid_optimize_using')
            use iso_c_binding
            type (c_ptr), value        :: grid
            integer (c_int32_t), value :: flags
        end subroutine

        integer (c_size_t) function grid_order_count(grid) bind(c, name = 'pineappl_grid_order_count')
            use iso_c_binding
            type (c_ptr), value :: grid
        end function

        subroutine grid_order_params(grid, order_params) bind(c, name = 'pineappl_grid_order_params')
            use iso_c_binding
            type (c_ptr), value :: grid
            integer (c_int32_t) :: order_params(*)
        end subroutine

        type (c_ptr) function grid_read(filename) bind(c, name = 'pineappl_grid_read')
            use iso_c_binding
            character (c_char) :: filename(*)
        end function

        subroutine grid_scale(grid, factor) bind(c, name = 'pineappl_grid_scale')
            use iso_c_binding
            type (c_ptr), value    :: grid
            real (c_double), value :: factor
        end subroutine

        subroutine grid_scale_by_bin(grid, count, factors) bind(c, name = 'pineappl_grid_scale_by_bin')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: count
            real (c_double)           :: factors(*)
        end subroutine

        subroutine grid_scale_by_order(grid, alphas, alpha, logxir, logxif, global) &
            bind(c, name = 'pineappl_grid_scale_by_order')
            use iso_c_binding
            type (c_ptr), value :: grid
            real (c_double), value :: alphas, alpha, logxir, logxif, global
        end subroutine

        subroutine grid_set_key_value(grid, key, valju) bind(c, name = 'pineappl_grid_set_key_value')
            use iso_c_binding
            type (c_ptr), value :: grid
            character (c_char)  :: key(*), valju(*)
        end subroutine

        subroutine grid_set_remapper(grid, dimensions, normalizations, limits) bind(c, name = 'pineappl_grid_set_remapper')
            use iso_c_binding
            type (c_ptr), value       :: grid
            integer (c_size_t), value :: dimensions
            real (c_double)           :: normalizations(*), limits(*)
        end subroutine

        subroutine grid_split_channels(grid) bind(c, name = 'pineappl_grid_split_channels')
            use iso_c_binding
            type (c_ptr), value :: grid
        end subroutine

        subroutine grid_write(grid, filename) bind(c, name = 'pineappl_grid_write')
            use iso_c_binding
            type (c_ptr), value :: grid
            character (c_char)  :: filename(*)
        end subroutine

        integer (c_size_t) function channels_combinations(channels, entry) bind(c, name = 'pineappl_channels_combinations')
            use iso_c_binding
            type (c_ptr), value       :: channels
            integer (c_size_t), value :: entry
        end function

        integer (c_size_t) function channels_count(channels) bind(c, name = 'pineappl_channels_count')
            use iso_c_binding
            type (c_ptr), value :: channels
        end function

        subroutine channels_delete(channels) bind(c, name = 'pineappl_channels_delete')
            use iso_c_binding
            type (c_ptr), value :: channels
        end subroutine

        subroutine channels_entry(channels, entry, pdg_ids, factors) bind(c, name = 'pineappl_channels_entry')
            use iso_c_binding
            type (c_ptr), value       :: channels
            integer (c_size_t), value :: entry
            integer (c_int32_t)       :: pdg_ids(*)
            real (c_double)           :: factors(*)
        end subroutine

        subroutine string_delete(string) bind(c, name = 'pineappl_string_delete')
            use iso_c_binding
            character (c_char) :: string(*)
        end subroutine
    end interface

contains
    ! https://stackoverflow.com/a/20121335 and https://community.intel.com/t5/Intel-Fortran-Compiler/Converting-c-string-to-Fortran-string/m-p/959515
    function c_f_string(c_str) result(f_str)
        use :: iso_c_binding

        type(c_ptr), intent(in) :: c_str
        character(kind=c_char), dimension(:), pointer :: arr_f_ptr => null()
        character(len=:, kind=c_char), allocatable :: f_str
        integer(kind=c_size_t) :: i, length

        length = strlen(c_str)
        call c_f_pointer(c_str, arr_f_ptr, [length])

        if (.not.associated(arr_f_ptr)) then
            f_str = "NULL"
            return
        end if

        allocate(character(len=length)::f_str)

        do i = 1, length
            f_str(i:i) = arr_f_ptr(i)
        end do
    end function

    type (pineappl_channels) function pineappl_channels_new(convolutions)
        implicit none

        integer (c_int32_t), value :: convolutions

        pineappl_channels_new = pineappl_channels(channels_new(convolutions))
    end function

    integer function pineappl_grid_bin_count(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        pineappl_grid_bin_count = int(grid_bin_count(grid%ptr))
    end function

    integer function pineappl_grid_bin_dimensions(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        pineappl_grid_bin_dimensions = int(grid_bin_dimensions(grid%ptr))
    end function

    function pineappl_grid_bin_limits_left(grid, dimension) result(res)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, intent(in)              :: dimension
        real (dp), allocatable           :: res(:)

        allocate(res(pineappl_grid_bin_count(grid)))

        call grid_bin_limits_left(grid%ptr, int(dimension, c_size_t), res)
    end function

    function pineappl_grid_bin_limits_right(grid, dimension) result(res)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, intent(in)              :: dimension
        real (dp), allocatable           :: res(:)

        allocate(res(pineappl_grid_bin_count(grid)))

        call grid_bin_limits_right(grid%ptr, int(dimension, c_size_t), res)
    end function

    function pineappl_grid_bin_normalizations(grid) result(res)
        implicit none

        type (pineappl_grid), intent(in) :: grid
        real (dp), allocatable           :: res(:)

        allocate(res(pineappl_grid_bin_count(grid)))

        call grid_bin_normalizations(grid%ptr, res)
    end function

    type (pineappl_grid) function pineappl_grid_clone(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        pineappl_grid_clone = pineappl_grid(grid_clone(grid%ptr))
    end function

    function pineappl_grid_convolve(grid, xfx, alphas, pdfs_state, alphas_state, order_mask, &
        channel_mask, bin_indices, nb_scales, mu_scales) result(res)

        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in)   :: grid
        type (pineappl_xfx), value         :: xfx
        type (pineappl_alphas)             :: alphas
        type (c_ptr), intent(in)           :: pdfs_state(:)
        type (c_ptr), intent(in)           :: alphas_state
        logical, intent(in)                :: order_mask(:), channel_mask(:)
        integer, intent(in)                :: bin_indices(:), nb_scales
        real (dp), intent(in)              :: mu_scales(:)
        real (dp), allocatable             :: res(:)
        integer                            :: i

        allocate(res(size(bin_indices)))

        if (.not. c_associated(c_funloc(xfx%proc))) then
            error stop "xfx%proc is null"
        end if

        if (.not. c_associated(c_funloc(alphas%proc))) then
            error stop "alphas%proc is null"
        end if

        call grid_convolve( &
            grid%ptr, &
            c_funloc(xfx%proc), &
            c_funloc(alphas%proc), &
            pdfs_state, &
            alphas_state, &
            [(logical(order_mask(i), c_bool), i = 1, size(order_mask))], &
            [(logical(channel_mask(i), c_bool), i = 1, size(channel_mask))], &
            [(int(bin_indices, c_size_t), i = 1, size(bin_indices))], &
            int(nb_scales, c_size_t), &
            mu_scales, &
            res &
        )

    end function

    subroutine pineappl_grid_delete(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        call grid_delete(grid%ptr)
    end subroutine

    subroutine pineappl_grid_fill2(grid, order, observable, channel, ntuple, weight)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in)    :: grid
        real (dp), intent(in)               :: observable, ntuple(*), weight
        integer, intent(in)                 :: order, channel

        call grid_fill2(grid%ptr, int(order, c_size_t), observable, int(channel, c_size_t), ntuple, weight)
    end subroutine

    subroutine pineappl_grid_fill_all2(grid, order, observable, ntuple, weights)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in)    :: grid
        integer, intent(in)                 :: order
        real (dp), intent(in)               :: observable, weights(*)
        real (dp), dimension(*), intent(in) :: ntuple

        call grid_fill_all2(grid%ptr, int(order, c_size_t), observable, ntuple, weights)
    end subroutine

    subroutine pineappl_grid_fill_array2(grid, orders, observables, ntuples, lumis, weights)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in)    :: grid
        real (dp), intent(in)               :: observables(*), ntuples(*), weights(*)
        integer, intent(in)                 :: orders(:), lumis(:)
        integer (c_size_t)                  :: i

        call grid_fill_array2(grid%ptr, [(int(orders(i), c_size_t), i = 1, size(orders))], &
            observables, ntuples, [(int(lumis(i), c_size_t), i = 1, size(lumis))], weights, int(size(orders), c_size_t))
    end subroutine

    function pineappl_grid_key_value(grid, key) result(res)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        character (*), intent(in)        :: key
        character (:), allocatable       :: res

        res = c_f_string(grid_key_value(grid%ptr, key // c_null_char))
    end function

    type (pineappl_channels) function pineappl_grid_channels(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        pineappl_grid_channels = pineappl_channels(grid_channels(grid%ptr))
    end function

    subroutine pineappl_grid_merge_and_delete(grid, other)
        implicit none

        type (pineappl_grid), intent(in) :: grid
        type (pineappl_grid), intent(in) :: other

        call grid_merge_and_delete(grid%ptr, other%ptr)
    end subroutine

    subroutine pineappl_grid_merge_bins(grid, from, to)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, intent(in)              :: from, to

        call grid_merge_bins(grid%ptr, int(from, c_size_t), int(to, c_size_t))
    end subroutine

    type (pineappl_grid) function pineappl_grid_new2(bins, bin_limits, orders, order_params, &
        channels, pid_basis, convolutions, interpolations, interp_info, kinematics, mu_scales)
        implicit none

        integer(kind(pineappl_pid_basis)), intent(in)                             :: pid_basis
        type (pineappl_channels), intent(in)                                      :: channels
        integer, intent(in)                                                       :: orders, bins, interpolations
        integer(int8), dimension(5 * orders), intent(in)                          :: order_params
        real (dp), dimension(bins + 1), intent(in)                                :: bin_limits
        type (pineappl_conv), dimension(*), intent(in), target                    :: convolutions
        type (pineappl_kinematics), dimension(interpolations), intent(in), target :: kinematics
        type (pineappl_interp), dimension(interpolations), intent(in)             :: interp_info
        type (pineappl_scale_func_form), dimension(interpolations)                :: mu_scales

        pineappl_grid_new2 = pineappl_grid(grid_new2(&
            int(bins, c_size_t), &
            bin_limits, &
            int(orders, c_size_t), &
            order_params, &
            channels%ptr, &
            pid_basis, &
            convolutions, &
            int(interpolations, c_size_t), &
            interp_info, &
            kinematics, &
            mu_scales) &
        )
    end function

    subroutine pineappl_grid_optimize(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        call grid_optimize(grid%ptr)
    end subroutine

    subroutine pineappl_grid_optimize_using(grid, flags)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, intent(in)              :: flags

        call grid_optimize_using(grid%ptr, int(flags, c_int32_t))
    end subroutine

    integer function pineappl_grid_order_count(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        pineappl_grid_order_count = int(grid_order_count(grid%ptr))
    end function

    function pineappl_grid_order_params(grid) result(res)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, allocatable             :: res(:)

        allocate(res(4 * pineappl_grid_order_count(grid)))

        call grid_order_params(grid%ptr, res)
    end function

    type (pineappl_grid) function pineappl_grid_read(filename)
        use iso_c_binding

        implicit none

        character (*), intent(in) :: filename

        pineappl_grid_read = pineappl_grid(grid_read(filename // c_null_char))
    end function

    subroutine pineappl_grid_scale(grid, factor)
        implicit none

        type (pineappl_grid), intent(in) :: grid
        real (dp), intent(in)            :: factor

        call grid_scale(grid%ptr, factor)
    end subroutine

    subroutine pineappl_grid_scale_by_bin(grid, factors)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        real (dp), intent(in)            :: factors(:)

        call grid_scale_by_bin(grid%ptr, int(size(factors), c_size_t), factors)
    end subroutine

    subroutine pineappl_grid_scale_by_order(grid, alphas, alpha, logxir, logxif, global)
        implicit none

        type (pineappl_grid), intent(in) :: grid
        real (dp), intent(in)            :: alphas, alpha, logxir, logxif, global

        call grid_scale_by_order(grid%ptr, alphas, alpha, logxir, logxif, global)
    end subroutine

    subroutine pineappl_grid_set_key_value(grid, key, valju)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        character (*), intent(in)        :: key, valju

        call grid_set_key_value(grid%ptr, key // c_null_char, valju // c_null_char)
    end subroutine

    subroutine pineappl_grid_set_remapper(grid, dimensions, normalizations, limits)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        integer, intent(in)              :: dimensions
        real (dp), intent(in)            :: normalizations(*), limits(*)

        call grid_set_remapper(grid%ptr, int(dimensions, c_size_t), normalizations, limits)
    end subroutine

    subroutine pineappl_grid_split_channels(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        call grid_split_channels(grid%ptr)
    end subroutine

    subroutine pineappl_grid_write(grid, filename)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        character (*), intent(in)        :: filename

        call grid_write(grid%ptr, filename // c_null_char)
    end subroutine

    subroutine pineappl_channels_add(channels, combinations, pdg_id_combinations, factors)
        use iso_c_binding

        implicit none

        type (pineappl_channels), intent(in)             :: channels
        integer, intent(in)                              :: combinations
        integer, dimension(2 * combinations), intent(in) :: pdg_id_combinations
        real (dp), dimension(combinations), intent(in)   :: factors

        call channels_add(channels%ptr, int(combinations, c_size_t), pdg_id_combinations, factors)
    end subroutine

    integer function pineappl_channels_combinations(channels, entry)
        use iso_c_binding

        implicit none

        type (pineappl_channels), intent(in) :: channels
        integer, intent(in)                  :: entry

        pineappl_channels_combinations = int(channels_combinations(channels%ptr, int(entry, c_size_t)))
    end function

    integer function pineappl_channels_count(channels)
        use iso_c_binding

        implicit none

        type (pineappl_channels), intent(in) :: channels

        pineappl_channels_count = int(channels_count(channels%ptr))
    end function

    subroutine pineappl_channels_delete(channels)
        implicit none

        type (pineappl_channels), intent(in) :: channels

        call channels_delete(channels%ptr)
    end subroutine

    subroutine pineappl_channels_entry(channels, entry, pdg_ids, factors)
        use iso_c_binding

        implicit none

        type (pineappl_channels), intent(in) :: channels
        integer, intent(in)                  :: entry
        integer, intent(out)                 :: pdg_ids(*)
        real (dp), intent(out)               :: factors(*)

        call channels_entry(channels%ptr, int(entry, c_size_t), pdg_ids, factors)
    end subroutine
end module
