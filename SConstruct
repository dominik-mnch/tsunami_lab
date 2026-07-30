import os
import SCons

print( '####################################' )
print( '### Tsunami Lab                  ###' )
print( '###                              ###' )
print( '### https://scalable.uni-jena.de ###' )
print( '####################################' )
print()
print('runnning build script')

# --- CUDA Detection ---

# 1. Try CUDA_PATH env variable
cuda_path = os.environ.get('CUDA_PATH','')

# 2. Fall back to common default location if not set
if not cuda_path:
        candidates = ['/usr/local/cuda', '/usr/cuda']
        for path in candidates:
                if os.path.isdir(path):
                        cuda_path = path
                        break

# 3. If still not found, disable CUDA support (graceful degradation for CI)
cuda_available = bool(cuda_path)
if not cuda_available:
        print("WARNING: CUDA not found. Building without GPU support.")
        print("  To enable GPU: set CUDA_PATH environment variable or install CUDA toolkit locally")
else:
        print(f"CUDA found at: {cuda_path}")

# --- CUDA Compiler Setup (only if available) ---

if cuda_available:
        # Path to nvcc binary
        nvcc = os.path.join(cuda_path, 'bin', 'nvcc')

        # GPU architecture. Defaults to sm_120 (Blackwell, e.g. RTX 50-series);
        # override via the CUDA_ARCH environment variable for other GPUs, e.g.
        # CUDA_ARCH=sm_89 for Ada (RTX 40-series) or sm_86 for Ampere.
        cuda_arch = os.environ.get('CUDA_ARCH', 'sm_120')
else:
        nvcc = None
        cuda_arch = None

# configuration
vars = Variables()

vars.AddVariables(
  EnumVariable( 'mode',
                'compile modes, option \'san\' enables address and undefined behavior sanitizers',
                'release',
                allowed_values=('release', 'debug', 'release+san', 'debug+san' )
              ),
  ('CXX', 'C++ compiler executable', 'g++'),
  EnumVariable( 'opt',
                'optimization level passed to the compiler as -O<opt>',
                '2',
                allowed_values=('0', '1', '2', '3', 'fast')
              ),
  BoolVariable( 'report',
                'enable compiler optimization remarks (vectorization/inlining reports)',
                False
              )
)

# exit in the case of unknown variables
if vars.UnknownVariables():
  print( "build configuration corrupted, don't know what to do with: " + str(vars.UnknownVariables().keys()) )
  exit(1)

# create environment
# ENV=os.environ forwards the surrounding shell environment (PATH,
# LD_LIBRARY_PATH, loaded modules, etc.) into SCons' subprocess environment,
# so e.g. `module load llvm` on the login/compute node is visible here too.
env = Environment( variables = vars, ENV = os.environ )
env.Append(CPPPATH=['src'])

# Task 1: generic compiler support via the CXX environment variable.
# Usage:   CXX=clang++ scons
#          CXX=g++     scons
# An explicit `scons CXX=...` command-line argument still wins over the
# environment variable, since it was already applied by vars.Update() above;
# we only fall back to the shell's CXX if the user did not pass it on the
# command line.
if 'CXX' in os.environ and 'CXX' not in ARGUMENTS:
  env['CXX'] = os.environ['CXX']

print(f"Using CXX = {env['CXX']}")

# Configure CUDA if available
if cuda_available:
  env.Append(CPPPATH=[os.path.join(cuda_path, 'include')])

  # Create an SCons Environment with CUDA settings
  env['NVCC'] = nvcc
  env['CUDA_ARCH'] = cuda_arch
  env['CUDA_PATH'] = cuda_path

  # Tell SCons how to compile .cu files using nvcc
  env['BUILDERS']['CUDAObject'] = Builder(
      action='$NVCC $NVCCFLAGS -c $SOURCE -o $TARGET',
      suffix='.o',
      src_suffix='.cu'
  )

  # CUDA compiler flags
  env['NVCCFLAGS'] = [
    '-std=c++17',
    '-O2',
    f'-arch={cuda_arch}',
    f'-I{os.path.join(cuda_path, "include")}',
    '-Isrc',
    '-Isrc/cuda',
    '-DCUDA_ENABLED',
    '--fmad=false',
    '--ftz=false',
    '--prec-div=true',
    '--prec-sqrt=true'
  ]

  # CUDA runtime library linking
  # The library directory differs between CUDA distributions: traditional
  # installs use 'lib64', while Nix-provided toolkits use 'lib'.
  cuda_libdirs = [d for d in ('lib64', 'lib')
                  if os.path.isdir(os.path.join(cuda_path, d))]
  if not cuda_libdirs:
    cuda_libdirs = ['lib64']
  env.Append( LIBS     = ['cudart'] )
  env.Append( LIBPATH  = [os.path.join(cuda_path, d) for d in cuda_libdirs] )

  # Define CUDA macro
  env.Append( CPPDEFINES = ['CUDA_ENABLED'] )

# env.CUDAObject('test_kernel.cu')

# The final programs are assembled from object files. Make the linker explicit so
# SCons uses the C++ compiler driver and links the C++ standard library reliably.
env.Replace( LINK = '$CXX' )

# generate help message
Help( vars.GenerateHelpText( env ) )

# add default flags
env.Append( CXXFLAGS = [ '-std=c++17',
                         '-Wall',
                         '-Wextra',
                         '-Wpedantic',
                         '-fopenmp',
                         '-ffp-contract=off' ] )
env.Append( LINKFLAGS = [ '-fopenmp' ] )

# add middle_states.csv location environment variable
env.Append(CPPDEFINES=[
  ('MIDDLE_STATES_CSV', '\\"' + Dir('#').abspath + '/res/middle_states.csv\\"'),
  ('BATHYMETRY_CSV', '\\"' + Dir('#').abspath + '/res/bathymetry.csv\\"'),
  ('ARTIFICIAL_TSUNAMI_BATHY_NC', '\\"' + Dir('#').abspath + '/res/artificialtsunami_bathymetry_1000.nc\\"'),
  ('ARTIFICIAL_TSUNAMI_DISP_NC', '\\"' + Dir('#').abspath + '/res/artificialtsunami_displ_1000.nc\\"'),
  ('STATIONS_CSV', '\\"' + Dir('#').abspath + '/stations/Stations.csv\\"'),
  ('STATIONS_OUTPUT_DIR', '\\"' + Dir('#').abspath + '/stations/output\\"')
])

# set optimization mode
# Task 3: the optimization level is now selectable via the `opt` build
# variable, e.g.  scons opt=fast   (-> -Ofast)   or   scons opt=3 (-> -O3)
# `mode=debug` still forces -O0 -g regardless of `opt`, matching the
# previous behaviour.
if 'debug' in env['mode']:
  env.Append( CXXFLAGS = [ '-g',
                           '-O0' ] )
else:
  env.Append( CXXFLAGS = [ '-O' + env['opt'] ] )

# add sanitizers
if 'san' in  env['mode']:
  env.Append( CXXFLAGS =  [ '-g',
                            '-fsanitize=float-divide-by-zero',
                            '-fsanitize=bounds',
                            '-fsanitize=address',
                            '-fsanitize=undefined',
                            '-fno-omit-frame-pointer' ] )
  env.Append( LINKFLAGS = [ '-g',
                            '-fsanitize=address',
                            '-fsanitize=undefined' ] )
else:
  env.Append( CXXFLAGS = [ '-Werror' ] )

# Task 4: optimization reports.
# Usage: scons report=1
# Clang uses -Rpass / -Rpass-missed / -Rpass-analysis remarks; GCC uses the
# -fopt-info family. We detect the compiler front-end from $CXX (this is a
# heuristic string check, not a compiler invocation) and emit the matching
# flags. Remarks go to stderr, so redirect the build to a file to inspect
# them, e.g.:  CXX=clang++ scons report=1 2> opt_report.txt
if env['report']:
  cxx_name = os.path.basename(env['CXX'])
  if 'clang' in cxx_name:
    env.Append( CXXFLAGS = [ '-Rpass=.*',
                             '-Rpass-missed=.*',
                             '-Rpass-analysis=.*' ] )
  else:
    # gcc (and g++-compatible front ends)
    env.Append( CXXFLAGS = [ '-fopt-info-all' ] )
  # Keep debug info so line numbers in the remarks are meaningful.
  env.Append( CXXFLAGS = [ '-gline-tables-only' ] if 'clang' in cxx_name else [ '-g' ] )

# add Catch2
env.Append( CXXFLAGS = [ '-isystem', 'submodules/Catch2/single_include' ] )

# add netCDF-C and netCDF-C++
env.ParseConfig('nc-config --libs --cflags')

# Draco's system library stack requires this explicitly for std::filesystem.
# Keep it after nc-config's libraries so the linker can resolve filesystem users.
env.Append( LIBS = [ 'stdc++fs' ] )

# get source files
VariantDir( variant_dir = 'build/src',
            src_dir     = 'src' )

env.sources = []
env.tests = []

Export('env')
SConscript( 'build/src/SConscript' )
Import('env')

env.Program( target = 'build/tsunami_lab',
             source = env.sources + env.standalone )

env.Program( target = 'build/benchmark',
             source = env.sources + env.benchmark )

env.Program( target = 'build/benchmark_grace',
             source = env.sources + env.benchmark_grace )

if hasattr( env, 'benchmark_cuda' ):
  env.Program( target = 'build/benchmark_cuda',
               source = env.sources + env.benchmark_cuda )

env.Program( target = 'build/tests',
             source = env.sources + env.tests )
