project(
  'cas',
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

cas_warning_flags = []
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
      cas_warning_flags += cflag
    endif
endforeach

################
# dependencies #
################

# threads
cas_thread_dep = dependency('threads')

# boost
cas_boost_dep = dependency('boost', include_type : 'system')

# pbcopper
cas_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

# pbbam
cas_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

# pancake
cas_pancake_dep = dependency('pancake', fallback : ['pancake', 'pancake_dep'])

# wfa
cas_wfa_dep = dependency('wfa', fallback : ['wfa', 'wfa_dep'])

cas_lib_deps = [
  cas_thread_dep,
  cas_boost_dep,
  cas_pbcopper_dep,
  cas_pbbam_dep,
  cas_pancake_dep,
  cas_wfa_dep,
]

###########
# sources #
###########

executable(
  'cas',
  files([
    'src/main.cpp',
  ]),
  dependencies : cas_lib_deps,
  install : true,
  cpp_args : cas_warning_flags
)
