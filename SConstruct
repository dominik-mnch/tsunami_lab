##
# @author Alexander Breuer (alex.breuer AT uni-jena.de)
#
# @section DESCRIPTION
# Entry-point for builds.
##
import SCons
import re
import subprocess

print( '####################################' )
print( '### Tsunami Lab                  ###' )
print( '###                              ###' )
print( '### https://scalable.uni-jena.de ###' )
print( '####################################' )
print()
print('runnning build script')

# configuration
vars = Variables()

vars.AddVariables(
  EnumVariable( 'mode',
                'compile modes, option \'san\' enables address and undefined behavior sanitizers',
                'release',
                allowed_values=('release', 'debug', 'release+san', 'debug+san' )
              ),
  ('CXX', 'C++ compiler executable', 'g++')
)

# exit in the case of unknown variables
if vars.UnknownVariables():
  print( "build configuration corrupted, don't know what to do with: " + str(vars.UnknownVariables().keys()) )
  exit(1)

# create environment
env = Environment( variables = vars )
env.Append(CPPPATH=['src'])

# The final programs are assembled from object files. Make the linker explicit so
# SCons uses the C++ compiler driver and links the C++ standard library reliably.
env.Replace( LINK = '$CXX' )

# generate help message
Help( vars.GenerateHelpText( env ) )

# add default flags
env.Append( CXXFLAGS = [ '-std=c++17',
                         '-Wall',
                         '-Wextra',
                         '-Wpedantic' ] )

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
if 'debug' in env['mode']:
  env.Append( CXXFLAGS = [ '-g',
                           '-O0' ] )
else:
  env.Append( CXXFLAGS = [ '-O2' ] )

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

# add Catch2
env.Append( CXXFLAGS = [ '-isystem', 'submodules/Catch2/single_include' ] )

# add netCDF-C and netCDF-C++
env.ParseConfig('nc-config --libs --cflags')

# GCC versions before 9 keep std::filesystem in a separate library.
def add_filesystem_library_if_needed( i_env ):
  try:
    l_version = subprocess.check_output( [ i_env.subst('$CXX'), '-dumpfullversion', '-dumpversion' ],
                                         stderr = subprocess.DEVNULL,
                                         universal_newlines = True ).strip().splitlines()[0]
    l_match = re.match( r'^(\d+)', l_version )
    if l_match and int( l_match.group(1) ) < 9:
      print( 'detected GCC ' + l_version + ', adding stdc++fs for std::filesystem' )
      i_env.Append( LIBS = [ 'stdc++fs' ] )
  except ( OSError, subprocess.CalledProcessError, IndexError ):
    pass

add_filesystem_library_if_needed( env )

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

env.Program( target = 'build/tests',
             source = env.sources + env.tests )
