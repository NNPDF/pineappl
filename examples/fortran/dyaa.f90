program dyaa
    use pineappl

    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type psp2to2
        real (dp) s, t, u, x1, x2, jacobian
    end type

    type (pineappl_channels) :: channels
    integer, dimension(2)    :: pdg_ids
    real (dp), dimension(1)  :: factors

    type (pineappl_grid)     :: grid

    integer :: i

    integer              :: order_idx, channel_idx, calls
    real (dp), parameter :: hbarc2 = 389379372.1_dp
    real (dp)            :: x1, x2, q2, weight, s, t, u, jacobian, ptl, mll, yll, ylp, ylm
    type (psp2to2)       :: tmp

    channels = pineappl_channels_new(2)
    ! create a new channel for the photon-photon initial state
    pdg_ids = [ 22, 22 ]
    factors = [ 1.0_dp ]
    call pineappl_channels_add(channels, 1, pdg_ids, factors)


    grid = pineappl_grid_new2( &
        ! number of bins
        24, &
        ! one-dimensional fill limits: we bin in rapidity from 0 to 2.4 in steps of 0.1
        [ (i * 0.1_dp, i = 0, 24) ], &
        ! number of orders
        1, &
        ! perturbative orders: only O(alpha^2)
        [ 0_1, 2_1, 0_1, 0_1, 0_1 ], &
        channels, &
        pineappl_pid_basis_pdg, &
        [ &
            pineappl_conv(pineappl_conv_type_unpol_pdf, 2212), &
            pineappl_conv(pineappl_conv_type_unpol_pdf, 2212) &
        ], &
        3, &
        [ &
            pineappl_interp(1e2_dp, 1e8_dp, 40, 3, pineappl_reweight_meth_no_reweight, pineappl_map_applgrid_h0, pineappl_interp_meth_lagrange), &
            pineappl_interp(2e-7_dp, 1.0_dp, 50, 3, pineappl_reweight_meth_applgrid_x, pineappl_map_applgrid_f2, pineappl_interp_meth_lagrange), &
            pineappl_interp(2e-7_dp, 1.0_dp, 50, 3, pineappl_reweight_meth_applgrid_x, pineappl_map_applgrid_f2, pineappl_interp_meth_lagrange) &
        ], &
        [ &
            pineappl_kinematics(pineappl_kinematics_tag_scale, 0), &
            pineappl_kinematics(pineappl_kinematics_tag_x, 0), &
            pineappl_kinematics(pineappl_kinematics_tag_x, 1) &
        ], &
        [ &
            pineappl_scale_func_form(pineappl_scale_func_form_tag_scale, pineappl_scale_func_form_body(0, 0)), &
            pineappl_scale_func_form(pineappl_scale_func_form_tag_scale, pineappl_scale_func_form_body(0, 0)), &
            pineappl_scale_func_form(pineappl_scale_func_form_tag_no_scale, pineappl_scale_func_form_body(0, 0)) &
        ] &
    )

    ! The `pineappl_scale_func_form_body` objects have to defined with two fields - if not required, the value(s) will be ignored

    call pineappl_channels_delete(channels)

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
        channel_idx = 0

        ! fill
        ! - 'grid'
        ! - for perturbative order O(alpha^2) ('order_idx = 0' corresponds to the first four powers
        !   given in 'orders' above)
        ! - for the bin of the differential distribution corresponding to 'abs(yll)'
        ! - for the first partonic channel ('channel_idx = 0' corresponds to the channel created
        !   above in 'channels')
        ! - for PDF parameters 'x1, x2, q2',
        ! - with the given 'weight'
        call pineappl_grid_fill2(grid, order_idx, abs(yll), channel_idx, [ x1, x2, q2 ], weight)
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

    ! optimize the grid representation (makes the file smaller)
    call pineappl_grid_optimize(grid)

    ! write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    call pineappl_grid_write(grid, 'DY-LO-AA-deprecated.pineappl.lz4')

    print *, 'Generated DY-LO-AA-deprecated.pineappl.lz4 containing a a -> l+ l-.'
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
