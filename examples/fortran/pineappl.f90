module pineappl
    use iso_c_binding, only: c_null_ptr, c_ptr

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

    interface
        subroutine grid_delete(grid) bind(c, name = 'pineappl_grid_delete')
            use iso_c_binding
            type (c_ptr), value :: grid
        end subroutine

        subroutine grid_fill(grid, x1, x2, q2, order, observable, lumi, weight) bind(c, name = 'pineappl_grid_fill')
            use iso_c_binding
            type (c_ptr), value       :: grid
            real (c_double), value    :: x1, x2, q2, observable, weight
            integer (c_size_t), value :: order, lumi
        end subroutine

        type (c_ptr) function grid_new(lumi, orders, order_params, bins, bin_limits, key_vals) &
                bind(c, name = 'pineappl_grid_new')
            use iso_c_binding
            type (c_ptr), value       :: lumi, key_vals
            integer (c_size_t), value :: orders, bins
            integer (c_int32_t)       :: order_params(*)
            real (c_double)           :: bin_limits(*)
        end function

        subroutine grid_optimize(grid) bind(c, name = 'pineappl_grid_optimize')
            use iso_c_binding
            type (c_ptr), value :: grid
        end subroutine

        subroutine grid_set_key_value(grid, key, valju) bind(c, name = 'pineappl_grid_set_key_value')
            use iso_c_binding
            type (c_ptr), value :: grid
            character (c_char)  :: key(*), valju(*)
        end subroutine

        subroutine grid_write(grid, filename) bind(c, name = 'pineappl_grid_write')
            use iso_c_binding
            type (c_ptr), value :: grid
            character (c_char)  :: filename(*)
        end subroutine

        subroutine keyval_delete(keyval) bind(c, name = 'pineappl_keyval_delete')
            use iso_c_binding
            type (c_ptr), value :: keyval
        end subroutine

        type (c_ptr) function keyval_new() bind(c, name = 'pineappl_keyval_new')
            use iso_c_binding
        end function

        subroutine keyval_set_string(keyval, key, valju) bind(c, name = 'pineappl_keyval_set_string')
            use iso_c_binding
            type (c_ptr), value :: keyval
            character (c_char)  :: key(*), valju(*)
        end subroutine

        subroutine lumi_add(lumi, combinations, pdg_id_pairs, factors) bind(c, name = 'pineappl_lumi_add')
            use iso_c_binding
            type (c_ptr), value       :: lumi
            integer (c_size_t), value :: combinations
            integer (c_int32_t)       :: pdg_id_pairs(*)
            real (c_double)           :: factors(*)
        end subroutine

        subroutine lumi_delete(lumi) bind(c, name = 'pineappl_lumi_delete')
            use iso_c_binding
            type (c_ptr), value :: lumi
        end subroutine

        type (c_ptr) function lumi_new() bind(c, name = 'pineappl_lumi_new')
            use iso_c_binding
        end function
    end interface

contains
    subroutine pineappl_grid_delete(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        call grid_delete(grid%ptr)
    end subroutine

    subroutine pineappl_grid_fill(grid, x1, x2, q2, order, observable, lumi, weight)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        real (dp), intent(in)            :: x1, x2, q2, observable, weight
        integer, intent(in)              :: order, lumi

        call grid_fill(grid%ptr, x1, x2, q2, int(order, c_size_t), observable, int(lumi, c_size_t), weight)
    end subroutine

    type (pineappl_grid) function pineappl_grid_new(lumi, orders, order_params, bins, bin_limits, key_vals)
        use iso_c_binding

        implicit none

        type (pineappl_lumi), intent(in)           :: lumi
        integer, intent(in)                        :: orders, bins
        integer, dimension(4 * orders), intent(in) :: order_params
        real (dp), dimension(bins + 1), intent(in) :: bin_limits
        type (pineappl_keyval), intent(in)         :: key_vals

        pineappl_grid_new = pineappl_grid(grid_new(lumi%ptr, int(orders, c_size_t), &
            order_params, int(bins, c_size_t), bin_limits, key_vals%ptr))
    end function

    subroutine pineappl_grid_optimize(grid)
        implicit none

        type (pineappl_grid), intent(in) :: grid

        call grid_optimize(grid%ptr)
    end subroutine

    subroutine pineappl_grid_set_key_value(grid, key, valju)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        character (*), intent(in)        :: key, valju

        call grid_set_key_value(grid%ptr, key // c_null_char, valju // c_null_char)
    end subroutine

    subroutine pineappl_grid_write(grid, filename)
        use iso_c_binding

        implicit none

        type (pineappl_grid), intent(in) :: grid
        character (*), intent(in)        :: filename

        call grid_write(grid%ptr, filename // c_null_char)
    end subroutine

    subroutine pineappl_keyval_delete(keyval)
        implicit none

        type (pineappl_keyval), intent(in) :: keyval

        call keyval_delete(keyval%ptr)
    end subroutine

    type (pineappl_keyval) function pineappl_keyval_new()
        implicit none

        pineappl_keyval_new = pineappl_keyval(keyval_new())
    end function

    subroutine pineappl_keyval_set_string(keyval, key, valju)
        use iso_c_binding

        implicit none
        type (pineappl_keyval), intent(in) :: keyval
        character (*), intent(in)          :: key, valju

        call keyval_set_string(keyval%ptr, key // c_null_char, valju // c_null_char)
    end subroutine

    subroutine pineappl_lumi_add(lumi, combinations, pdg_id_pairs, factors)
        use iso_c_binding

        implicit none

        type (pineappl_lumi), intent(in)                 :: lumi
        integer, intent(in)                              :: combinations
        integer, dimension(2 * combinations), intent(in) :: pdg_id_pairs
        real (dp), dimension(combinations), intent(in)   :: factors

        call lumi_add(lumi%ptr, int(combinations, c_size_t), pdg_id_pairs, factors)
    end subroutine

    subroutine pineappl_lumi_delete(lumi)
        implicit none

        type (pineappl_lumi), intent(in) :: lumi

        call lumi_delete(lumi%ptr)
    end subroutine

    type (pineappl_lumi) function pineappl_lumi_new()
        implicit none

        pineappl_lumi_new = pineappl_lumi(lumi_new())
    end function
end module
