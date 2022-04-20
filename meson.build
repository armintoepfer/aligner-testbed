project(
  'aligner-testbed',
  'cpp',
  version : '0.0.1',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++20',
    'c_std=c17',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.57.0')

############
# CXXFLAGS #
############

at_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wduplicated-cond',
  '-Wduplicated-branches',
  '-Wlogical-op',
  '-Wrestrict',
  '-Wold-style-cast',
  '-Wuseless-cast',
  '-Wdouble-promotion',
  '-Wunused-local-typedef',
  '-Wunused-parameter',
  '-Wno-sign-compare',
  '-Wshadow',
  '-Wformat=1']
    if cpp.has_argument(cflag)
      at_warning_flags += cflag
    endif
endforeach

################
# dependencies #
################

# threads
at_thread_dep = dependency('threads')

# boost
at_boost_dep = dependency('boost', include_type : 'system')

# pbcopper
at_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

# pbbam
at_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

# pancake
at_pancake_dep = dependency('pancake', fallback : ['pancake', 'pancake_dep'])

# wfa
at_wfa_dep = dependency('wfa', fallback : ['wfa', 'wfa_dep'])

# miniwfa
at_miniwfa_dep = dependency('miniwfa', fallback : ['miniwfa', 'miniwfa_dep'])

at_lib_deps = [
  at_thread_dep,
  at_boost_dep,
  at_pbcopper_dep,
  at_pbbam_dep,
  at_pancake_dep,
  at_wfa_dep,
  at_miniwfa_dep,
]

###########
# sources #
###########

executable(
  'at',
  files([
    'src/main.cpp',
  ]),
  dependencies : at_lib_deps,
  install : true,
  cpp_args : at_warning_flags
)
