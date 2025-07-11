module callbacks4
    implicit none
contains
    function wrap_xfx(pdg_id, x, q2, state) bind(c)
        use iso_c_binding

        implicit none

        integer(c_int32_t), value, intent(in) :: pdg_id
        real(c_double), value, intent(in)     :: x, q2
        type(c_ptr), value, intent(in)        :: state
        real(c_double)                        :: wrap_xfx
        integer, pointer                      :: state_array(:)

        call c_f_pointer(state, state_array, [2])
        call lhapdf_xfxq2(state_array(1), state_array(2), pdg_id, x, q2, wrap_xfx)
    end function

    function wrap_alphasq2(q2, state) bind(c)
        use iso_c_binding

        implicit none

        real(c_double), value, intent(in) :: q2
        type(c_ptr), value, intent(in)    :: state
        real(c_double)                    :: wrap_alphasq2
        integer, pointer                      :: state_array(:)

        call c_f_pointer(state, state_array, [2])
        call lhapdf_alphasq2(state_array(1), state_array(2), q2, wrap_alphasq2)
    end function
end module

program lhapdf_example
    use iso_c_binding
    use pineappl
    use callbacks4

    implicit none

    integer, parameter :: dp = kind(0.0d0)

    type(pineappl_grid)            :: grid
    type(pineappl_channels)        :: channels
    type(pineappl_kinematics)      :: kinematics(3)
    type(pineappl_scale_func_form) :: mu_scales_form(3)
    type(pineappl_interp)          :: interp_info(3)
    type(pineappl_conv)            :: convolutions(2)

    type(pineappl_xfx) :: xfx
    type(pineappl_alphas) :: alphas

    integer(kind(pineappl_reweight_meth)) :: q2_reweight
    integer(kind(pineappl_reweight_meth)) :: x_reweight
    integer(kind(pineappl_map)) :: q2_mapping
    integer(kind(pineappl_map)) :: x_mapping
    integer(kind(pineappl_interp_meth)) :: interpolation_meth

    integer, target        :: alphas_flags(2)
    type(c_ptr), target    :: pdfs_state(2)
    integer(c_int), target :: pdfs_array(2,2)
    character(len=30)      :: pdfset1, pdfset2

    channels = pineappl_channels_new(2) ! The argument is the number of convolutions
    call pineappl_channels_add(channels, 3, [0, 0, 1, -1, 2, -2], [1.0_dp, 1.0_dp, 1.0_dp])

    kinematics = [&
        pineappl_kinematics(pineappl_kinematics_tag_scale, 0), &
        pineappl_kinematics(pineappl_kinematics_tag_x, 0), &
        pineappl_kinematics(pineappl_kinematics_tag_x, 1) &
    ]

    q2_reweight = pineappl_reweight_meth_no_reweight
    x_reweight = pineappl_reweight_meth_applgrid_x
    q2_mapping = pineappl_map_applgrid_h0
    x_mapping = pineappl_map_applgrid_f2
    interpolation_meth = pineappl_interp_meth_lagrange
    interp_info = [ &
        pineappl_interp(1e2_dp, 1e8_dp, 40, 3, q2_reweight, q2_mapping, interpolation_meth), &
        pineappl_interp(2e-7_dp, 1.0_dp, 50, 3, x_reweight, x_mapping, interpolation_meth), &
        pineappl_interp(2e-7_dp, 1.0_dp, 50, 3, x_reweight, x_mapping, interpolation_meth) &
    ]

    ! The `pineappl_scale_func_form_body` objects have to defined with two fields - if not required, the value(s) will be ignored
    mu_scales_form = [ &
        pineappl_scale_func_form(pineappl_scale_func_form_tag_scale, pineappl_scale_func_form_body(0, 0)), &
        pineappl_scale_func_form(pineappl_scale_func_form_tag_scale, pineappl_scale_func_form_body(0, 0)), &
        pineappl_scale_func_form(pineappl_scale_func_form_tag_no_scale, pineappl_scale_func_form_body(0, 0)) &
    ]

    convolutions = [ &
        pineappl_conv(pineappl_conv_type_unpol_pdf, 2212), &
        pineappl_conv(pineappl_conv_type_unpol_pdf, 2212) &
    ]

    grid = pineappl_grid_new2(2, [0.0_dp, 1.0_dp, 2.0_dp], 1, [2_1, 0_1, 0_1, 0_1, 0_1], channels, pineappl_pid_basis_pdg, &
        convolutions, 3, interp_info, kinematics, mu_scales_form)

    call pineappl_grid_fill_all2(grid, 0, 0.5_dp, [100.0_dp, 0.5_dp, 0.5_dp], [0.5_dp, 0.5_dp, 0.5_dp])
    call pineappl_grid_fill_all2(grid, 0, 1.5_dp, [100.0_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 1.5_dp])

    call setlhaparm("SILENT")
    pdfset1 = "NNPDF31_nlo_as_0118_luxqed"
    pdfset2 = "MSHT20qed_nnlo"
    call lhapdf_initpdfset_byname(0, trim(pdfset1)) ! Init 1st PDF with ID=0
    call lhapdf_initpdfset_byname(1, trim(pdfset2)) ! Init 2nd PDF with ID=1

    ! Construct the callable to the function `xfx` and `alphasQ2`
    xfx = pineappl_xfx(wrap_xfx)
    alphas = pineappl_alphas(wrap_alphasq2)

    ! Define the array used to select the PDF and member ID.
    ! The array is of the form [[ISET, IMEMBER], ...] where the first element represents
    ! the 1st PDF set and is a tuple containing the set identification and the replica id.
    pdfs_array = reshape([0, 0, 0, 0], [2,2])
    pdfs_state(1) = c_loc(pdfs_array(1,1))
    pdfs_state(2) = c_loc(pdfs_array(1,2))

    ! [ISET, IMEMBER] for the computation of alphasQ2
    ! Here we first choose the 1st PDF to compute the alphasQ2
    alphas_flags = [0, 0]
    print *, "Computing predictions with the same PDF: ", trim(pdfset1)
    write(*, *) pineappl_grid_convolve(grid, xfx, alphas, pdfs_state, c_loc(alphas_flags(1)), &
        [.true.], [.true.], [0, 1], 1, [1.0_dp, 1.0_dp, 1.0_dp])

    pdfs_array = reshape([0, 0, 1, 0], [2,2])
    pdfs_state(1) = c_loc(pdfs_array(1,1))
    pdfs_state(2) = c_loc(pdfs_array(1,2))
    print *, "Computing predictions with different PDFs and alphasQ2(", trim(pdfset1), "):"
    write(*, *) pineappl_grid_convolve(grid, xfx, alphas, pdfs_state, c_loc(alphas_flags(1)), &
        [.true.], [.true.], [0, 1], 1, [1.0_dp, 1.0_dp, 1.0_dp])

    ! [ISET, IMEMBER] for the computation of alphasQ2
    ! Here we first choose the 1st PDF to compute the alphasQ2
    alphas_flags = [1, 0]
    print *, "Computing predictions with different PDFs and alphasQ2(", trim(pdfset2), "):"
    write(*, *) pineappl_grid_convolve(grid, xfx, alphas, pdfs_state, c_loc(alphas_flags(1)), &
        [.true.], [.true.], [0, 1], 1, [1.0_dp, 1.0_dp, 1.0_dp])

    ! call pineappl_grid_write(grid, 'test.pineappl.lz4')
    call pineappl_channels_delete(channels)
    call pineappl_grid_delete(grid)
end program lhapdf_example
