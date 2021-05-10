project(
    'appl2pine',
    'cpp',
    version : '0.1',
    license : 'GPL3',
    default_options : [ 'cpp_std=c++17', 'warning_level=3' ],
)

applgrid_config = find_program('applgrid-config')

meson.override_dependency('APPLgrid', declare_dependency(
    compile_args : run_command(applgrid_config, '--cxxflags').stdout().split(),
    link_args : run_command(applgrid_config, '--ldflags').stdout().split(),
    version : run_command(applgrid_config, '--version').stdout(),
))

applgrid_dep = dependency('APPLgrid', version : '=1.4.70')
pineappl_capi_dep = dependency('pineappl_capi', version : '>=0.3')
lhapdf_dep = dependency('lhapdf')

bin = executable(
    'appl2pine',
    'appl2pine.cpp',
    dependencies : [ applgrid_dep, lhapdf_dep, pineappl_capi_dep ],
    install : true
)
