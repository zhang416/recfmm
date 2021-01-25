project('recFMM', ['c', 'fortran'],
        version : '0.1.1',
        license : 'GPL-3.0-or-later',
        meson_version : '>=0.49.0',
        default_options : ['c_std=gnu99',
                           'buildtype=release',
                           'default_library=both',
                           'prefix=' + meson.current_source_dir()])

vso = '0'

# ========== compiler options and dependencies

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')

cc_is_icc = cc.get_id().contains('intel')
use_opencilk = not cc_is_icc and cc.has_argument('-fopencilk')
use_cilkplus = not use_opencilk \
               and (cc_is_icc or cc.has_argument('-fcilkplus'))

# HACK GNU gcc 8.x and 9.x still accept `-fcilkplus` even though they do not
# actually support Cilk...
cc_is_gcc_with_cilk = cc.get_id() == 'gcc' and \
                      cc.version().version_compare('>=5.0.0') and \
                      cc.version().version_compare('<8.1.0')
if use_cilkplus and cc.get_id() == 'gcc' and not cc_is_gcc_with_cilk
  warning('The specified GCC version does not support Cilk' +
          ' (even though -fcilkplus does not trigger an error).')
  use_cilkplus = false
endif

if use_opencilk or use_cilkplus
  have_cilk_headers = cc.has_header('cilk/cilk.h') and \
                      cc.has_header('cilk/cilk_api.h') and \
                      cc.has_header('cilk/reducer_min.h') and \
                      cc.has_header('cilk/reducer_max.h')
  if not have_cilk_headers
    warning('Required Cilk headers not found. Switching to sequential build.')
  endif
endif

if cc_is_icc
  add_project_arguments(['-Dicc'], language : 'c') # enables SIMD reduction
endif

m_dep = cc.find_library('m', required : false)

# ========== parallel vs sequential implementation

flag_parallel = (use_opencilk or use_cilkplus) and have_cilk_headers

if flag_parallel                # ----- parallel build

  # Cilk flags
  if use_opencilk
    cilk_cflags  = ['-fopencilk']
    cilk_ldflags = ['-fopencilk']
  elif use_cilkplus and not cc_is_icc
    cilk_cflags  = ['-fcilkplus']
    cilk_ldflags = ['-lcilkrts']
  else
    cilk_cflags  = []
    cilk_ldflags = []
  endif
  add_project_arguments([cilk_cflags, '-DUSING_CILK'], language : 'c')
  add_project_link_arguments([cilk_ldflags], language : 'c')

  # Intel TBB scalable allocator
  tbb_dir = get_option('tbb_dir')
  if tbb_dir == ''
    tbb_dir = []
  endif
  tbbmalloc_dep = dependency('tbbmalloc', required : false)
  if not tbbmalloc_dep.found()
    tbbmalloc_dep = cc.find_library('tbbmalloc', required : true,
                                    dirs : tbb_dir)
  endif

else                            # ----- sequential build

  warning('Could not verify Cilk support with the specified compiler.' +
          ' Cilk keywords elided.')

  tbbmalloc_dep = dependency('', required : false)

endif

recfmm_dep = [tbbmalloc_dep, m_dep]

# ========== build targets

# headers
subdir('include')
install_headers(recfmm_hdrs, install_dir : 'include')

# source files
subdir('src')
subdir('test')
subdir('example')

def_lap = ['-DLAPLACE']
def_yuk = ['-DYUKAWA']

# recFMM libraries: Laplace & Yukawa
lib_lap = library('adap_laplace', [src_recfmm_laplace],
                  c_args : def_lap,
                  dependencies : recfmm_dep,
                  include_directories : [recfmm_include_dir],
                  install : true, install_dir : 'lib',
                  soversion : vso)
lib_yuk = library('adap_yukawa', [src_recfmm_yukawa],
                  c_args : def_yuk,
                  dependencies : recfmm_dep,
                  include_directories : [recfmm_include_dir],
                  install : true, install_dir : 'lib',
                  soversion : vso)

# test executables
exe_test_lap = executable('recfmm-test-laplace', [src_recfmm_test],
                          c_args : def_lap, link_with : lib_lap,
                          dependencies : recfmm_dep,
                          include_directories : [recfmm_include_dir],
                          install_rpath : get_option('prefix') / 'lib',
                          install : true, install_dir : 'test')
exe_test_yuk = executable('recfmm-test-yukawa', [src_recfmm_test],
                          c_args : def_yuk, link_with : lib_yuk,
                          dependencies : recfmm_dep,
                          include_directories : [recfmm_include_dir],
                          install_rpath : get_option('prefix') / 'lib',
                          install : true, install_dir : 'test')

# example executable
exe_example_lap = executable('recfmm-example', [src_recfmm_example],
                             c_args : def_lap, link_with : lib_lap,
                             dependencies : recfmm_dep,
                             include_directories : [recfmm_include_dir],
                             install_rpath : get_option('prefix') / 'lib',
                             install : true, install_dir : 'example')

# ========== tests

test_verifier = find_program('test/verify_test_pass.sh')

test('Laplace [box, acc=3]', test_verifier,
     args : [exe_test_lap, '-d1', '-a3'],
     suite : ['laplace', 'box', 'low_acc'])
test('Laplace [sphere, acc=3]', test_verifier,
     args : [exe_test_lap, '-d2', '-a3'],
     suite : ['laplace', 'sphere', 'low_acc'])
test('Laplace [box, acc=6]', test_verifier,
     args : [exe_test_lap, '-d1', '-a6'],
     suite : ['laplace', 'box', 'high_acc'])
test('Laplace [sphere, acc=6]', test_verifier,
     args : [exe_test_lap, '-d2', '-a6'],
     suite : ['laplace', 'sphere', 'high_acc'])

test('Yukawa [box, acc=3]', test_verifier,
     args : [exe_test_yuk, '-d1', '-a3'],
     suite : ['yukawa', 'box', 'low_acc'])
test('Yukawa [sphere, acc=3]', test_verifier,
     args : [exe_test_yuk, '-d2', '-a3'],
     suite : ['yukawa', 'sphere', 'low_acc'])
test('Yukawa [box, acc=6]', test_verifier,
     args : [exe_test_yuk, '-d1', '-a6'],
     suite : ['yukawa', 'box', 'high_acc'])
test('Yukawa [sphere, acc=6]', test_verifier,
     args : [exe_test_yuk, '-d2', '-a6'],
     suite : ['yukawa', 'sphere', 'high_acc'])
