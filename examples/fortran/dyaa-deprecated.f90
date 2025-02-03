program dyaa
    use pineappl

    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type psp2to2
        real (dp) s, t, u, x1, x2, jacobian
    end type

    type (pineappl_lumi)    :: lumi
    integer, dimension(2)   :: pdg_ids
    real (dp), dimension(1) :: ckm_factors

    type (pineappl_grid)     :: grid
    integer, dimension(4)    :: orders
    real (dp), dimension(25) :: bins

    type (pineappl_keyval) :: key_vals

    integer :: i

    integer              :: order_idx, lumi_idx, calls
    real (dp), parameter :: hbarc2 = 389379372.1_dp
    real (dp)            :: x1, x2, q2, weight, s, t, u, jacobian, ptl, mll, yll, ylp, ylm
    type (psp2to2)       :: tmp

    ! create a new luminosity function for the photon-photon initial state
    lumi = pineappl_lumi_new()
    pdg_ids = [ 22, 22 ]
    ckm_factors = [ 1.0_dp ]
    call pineappl_lumi_add(lumi, 1, pdg_ids, ckm_factors)

    ! only O(alphas^0 alpha^2 log^0(xiR^2) \log^0(xiF^2)
    orders = [ 0, 2, 0, 0 ]
    ! we bin in rapidity from 0 to 2.4 in steps of 0.1
    bins = [ (i * 0.1_dp, i = 0, 24) ]

    ! create the PineAPPL grid with default interpolation and binning parameters
    key_vals = pineappl_keyval_new()
    grid = pineappl_grid_new(lumi, 1, orders, 24, bins, key_vals)

    call pineappl_keyval_delete(key_vals)
    call pineappl_lumi_delete(lumi)

    ! number of phase-space points that are generated before cuts
    calls = 10000000

    ! fill the grid with phase-space points
    do i = 1, calls
        tmp = hadronic_pspgen(10.0_dp, 7000.0_dp)
        s = tmp%s
        t = tmp%t
        u = tmp%u
        x1 = tmp%x1
        x2 = tmp%x2
        jacobian = tmp%jacobian

        ptl = sqrt((t * u / s))
        mll = sqrt(s)
        yll = 0.5_dp * log(x1 / x2)
        ylp = abs(yll + acosh(0.5_dp * mll / ptl))
        ylm = abs(yll - acosh(0.5_dp * mll / ptl))

        jacobian = jacobian * hbarc2 / calls

        ! cuts for LO for the invariant-mass slice containing the Z-peak from the paper given in the
        ! metadata below
        if ((ptl < 14.0_dp) .or. (abs(yll) > 2.4_dp) .or. (ylp > 2.4_dp) .or. (ylm > 2.4_dp) .or. &
            (mll < 60.0_dp) .or. (mll > 120.0_dp)) then
            cycle
        end if

        weight = jacobian * int_photo(s, t, u)
        ! renormalisation and factorisation scale
        q2 = 90.0_dp**2
        order_idx = 0
        lumi_idx = 0

        ! fill
        ! - 'grid'
        ! - for PDF parameters 'x1, x2, q2',
        ! - for perturbative order O(alpha^2) ('order_idx = 0' corresponds to the first four powers
        !   given in 'orders' above)
        ! - for the bin of the differential distribution corresponding to 'abs(yll)'
        ! - for the first partonic channel ('lumi_idx = 0' corresponds to the lumi entry created
        !   above in 'lumi')
        ! - with the given 'weight'
        call pineappl_grid_fill(grid, x1, x2, q2, order_idx, abs(yll), lumi_idx, weight)
    end do

    ! set metadata - this isn't strictly needed, but usually useful (plot script, ...)
    call pineappl_grid_set_key_value(grid, 'arxiv', '1310.7291')
    call pineappl_grid_set_key_value(grid, 'description', 'CMS double-differential Drell—Yan cross section at 7 TeV')
    call pineappl_grid_set_key_value(grid, 'hepdata', '10.17182/hepdata.62207.v1/t8')
    call pineappl_grid_set_key_value(grid, 'x1_label', 'yll')
    call pineappl_grid_set_key_value(grid, 'x1_label_tex', '$y_{\ell\bar{\ell}}$')
    ! rapidity doesn't have a unit (other observables could be GeV, TeV, ...)
    call pineappl_grid_set_key_value(grid, 'x1_unit', '')
    call pineappl_grid_set_key_value(grid, 'y_unit', 'pb')

    ! the following are the default values
    !call pineappl_grid_set_key_value(grid, 'initial_state_1', '2212') ! proton
    !call pineappl_grid_set_key_value(grid, 'initial_state_2', '2212') ! proton

    ! optimize the grid representation (makes the file smaller)
    call pineappl_grid_optimize(grid)

    ! write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    call pineappl_grid_write(grid, 'DY-LO-AA-deprecated.pineappl.lz4')

    print *, 'Generated DY-LO-AA-deprecated.pineappl.lz4 containing a a -> l+ l-.'
    print *, ''
    print *, 'Try running (PDF sets must contain non-zero photon PDF):'
    print *, '  - pineappl convolve DY-LO-AA.pineappl.lz4 NNPDF31_nnlo_as_0118_luxqed'
    print *, '  - pineappl --silence-lhapdf plot DY-LO-AA-deprecated.pineappl.lz4 NNPDF31_nnlo_as_0118_luxqed MSHT20qed_nnlo > plot_script.py'
    print *, '  - pineappl --help'

    call pineappl_grid_delete(grid)

contains
    ! photon-photon initiated lepton-pair production
    real (dp) function int_photo(s, t, u)
        implicit none
        real (dp), intent(in) :: s, t, u
        real (dp), parameter :: alpha0 = 1.0_dp / 137.03599911_dp

        int_photo = alpha0 * alpha0 / 2.0_dp / s * (t / u + u / t)
    end function

    ! phase-space generator for 2->2 scattering
    type (psp2to2) function hadronic_pspgen(mmin, mmax)
        implicit none

        real (dp) :: mmin, mmax, smin, smax, r1, r2, r3, tau0, tau, y, x1, x2, s, t, u, jacobian, cos_theta

        smin = mmin * mmin
        smax = mmax * mmax

        call random_number(r1)
        call random_number(r2)
        call random_number(r3)

        tau0 = smin / smax
        tau = tau0**r1
        y = tau**(1.0_dp - r2)
        x1 = y
        x2 = tau / y
        s = tau * smax

        jacobian = tau * log(tau0) * log(tau0) * r1

        ! theta integration (in the CMS)
        cos_theta = 2.0_dp * r3 - 1.0_dp
        jacobian = jacobian * 2.0_dp

        t = -0.5_dp * s * (1.0_dp - cos_theta)
        u = -0.5_dp * s * (1.0_dp + cos_theta)

        ! phi integration
        jacobian = jacobian * 2.0_dp * acos(-1.0_dp)

        hadronic_pspgen = psp2to2(s, t, u, x1, x2, jacobian)
    end function
end program
