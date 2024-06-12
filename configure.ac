# Initialize the package and version
AC_INIT([ICESat2VegR],[0.1.33])

# Initialize autoconfig
AC_PREREQ([2.71])

# Check for the required version of pkg-config
PKG_PROG_PKG_CONFIG([0.15.0])

# Pick R specific compiler and flags
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_MSG_CHECKING([the operating system name])
os_name=`"${R_HOME}/bin/Rscript" -e "cat(Sys.info()[['sysname']])"`
AC_MSG_RESULT([$os_name])

# Check for the presence of a C compiler
AC_PROG_CXX()


#################
## GDAL
#################
# Define a function to give an error message if GDAL is not installed
AC_DEFUN([GDAL_INSTALLATION_ERROR], [
  AC_MSG_ERROR([
      Please make sure the GDAL development package is installed:
      
      - On Debian/Ubuntu: sudo apt-get install libgdal-dev
      - On Fedora: sudo dnf install gdal-devel
      - On Arch Linux: sudo pacman -S gdal
      - On macOS (Homebrew): brew install gdal
      - On macOS (MacPorts): sudo port install gdal
  ])
])

# Define a function to test if gdal test program compiles
# args: message, action if compiles, action if fails
AC_DEFUN([GDAL_TEST_PROGRAM_IFELSE], [
  AC_MSG_CHECKING([$1])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([
    AC_LANG_PROGRAM(
        [[#include <gdal.h>]],
        [[GDALAllRegister();]]
    )
  ]
  ,[
      AC_MSG_RESULT([yes])
      $2
    ], # Action if the program compiles correctly
    [
      AC_MSG_RESULT([no])
      $3
  ]   # Action if the program fails to compile
  )
  AC_LANG_POP([C++])
])


# Check for GDAL using pkg-config
PKG_CHECK_MODULES([GDAL], [gdal])
PKG_CHECK_MODULES_STATIC([GDAL_STATIC], [gdal])

GDAL_PKG_CONFIG_CFLAGS=$GDAL_CFLAGS
GDAL_PKG_CONFIG_LIBS=$GDAL_LIBS

# Check for GDAL gdal-config
# Allow the user to specify the path to gdal-config
AC_ARG_WITH([gdal-config],
  [AS_HELP_STRING([--with-gdal-config=PATH], [specify the path to gdal-config])],
  [gdal_config_path=$withval], [gdal_config_path=""])

# Check if the user has specified a custom path
AS_IF([test -n "$gdal_config_path"],
  [
    # Check if the gdal-config exists at the specified path
  AC_CHECK_PROG([GDAL_CONFIG], [$gdal_config_path], [$gdal_config_path])
  AS_IF([test -z "$GDAL_CONFIG"], [
      AC_MSG_ERROR([gdal-config not found at specified path: $gdal_config_path])])],
  [
    # Check for gdal-config in the default PATH
    AC_PATH_PROG([GDAL_CONFIG], [gdal-config])
  
  # If not found, error out with instructions
  AS_IF([test -z "$GDAL_CONFIG"], [
    AC_MSG_NOTICE([
      gdal-config not found.

      If you have installed the GDAL libraries, then make sure that
      gdal-config is in your path. Try typing gdal-config at a
      shell prompt and see if it runs. If not, use:
        --configure-args='--with-gdal-config=/usr/local/bin/gdal-config'
      with appropriate values for your installation.

      If you have not installed the GDAL libraries,
    ])
    GDAL_INSTALLATION_ERROR
  ]
  )])

# Output the found path of gdal-config
AC_MSG_NOTICE([Using gdal-config at $GDAL_CONFIG])

# Get compiler and linker flags from gdal-config
GDAL_CONFIG_CFLAGS=`${GDAL_CONFIG} --cflags`
GDAL_CONFIG_LIBS=$($GDAL_CONFIG --libs)
GDAL_CONFIG_DEP_LIBS=$($GDAL_CONFIG --dep-libs)


# Define CXXFLAGS which are used by default to compile the test program
CXXFLAGS="${CXXFLAGS} ${GDAL_CFLAGS}"

# Define CPPFLAGS which are used by default by AC_CHECK_HEADER
CPPFLAGS="${CPPFLAGS} ${GDAL_CFLAGS}" 

# Check if gdal.h header file is present
AC_CHECK_HEADER([gdal.h],
  [],
  [
    AC_MSG_NOTICE([gdal.h not found in given locations.])
    GDAL_INSTALLATION_ERROR
  ]
)

GDAL_LIBS=$GDAL_CONFIG_LIBS
LIBS="${GDAL_LIBS}"
gdal_can_compile=no

# Try to link GDAL with gdal-config --libs only
GDAL_LIBS=$GDAL_CONFIG_LIBS
LIBS="${GDAL_LIBS}"

GDAL_TEST_PROGRAM_IFELSE(
  [whether GDAL can be linked using only `gdal-config --libs`],
  [
    gdal_can_compile=yes
  ]
)

# Try to link GDAL with gdal-config --libs and --dep-libs
AS_IF([test "x$gdal_can_compile" = "xno"], [
  GDAL_LIBS="${GDAL_CONFIG_LIBS} ${GDAL_CONFIG_DEP_LIBS}"
  LIBS="${GDAL_LIBS}"
  
  GDAL_TEST_PROGRAM_IFELSE(
    [whether GDAL can be linked using `gdal-config --libs` and `--dep-libs`],
    [
      gdal_can_compile=yes
    ]
  )
])

# Try to link GDAL with `pkg-config --libs gdal` only
AS_IF([test "x$gdal_can_compile" = "xno"], [
  GDAL_LIBS=$GDAL_PKG_CONFIG_LIBS
  LIBS="${GDAL_LIBS}"


  GDAL_TEST_PROGRAM_IFELSE(
    [whether GDAL can be linked using only `pkg-config --libs gdal`],
    [
      gdal_can_compile=yes
    ]
  )
])

AS_IF([test "x$gdal_can_compile" = "xno"], [
  # Now try to link GDAL with pkg-config --libs --static gdal
  GDAL_LIBS="${$GDAL_STATIC_LIBS}"
  LIBS="${GDAL_LIBS}"
  
  GDAL_TEST_PROGRAM_IFELSE(
    [whether GDAL can be linked using `pkg-config --libs --static gdal`],
    [
      gdal_can_compile=yes
    ]
  )
])

AS_IF([test "x$gdal_can_compile" = "xno"], [
      AC_MSG_NOTICE([
        Could not link to GDAL using `pkg-config` nor `gdal-config`.
      ])
      GDAL_INSTALLATION_ERROR
    ]
)

AC_SUBST([GDAL_LIBS])

################
## PROJ
################
PKG_CHECK_MODULES([PROJ], [proj])

# Copy PROJ datadir
PKG_CHECK_VAR([PROJ_DATADIR], [proj], [datadir])

AC_MSG_CHECKING([proj datadir path])
AC_MSG_RESULT([$PROJ_DATADIR])

# Try to use the PROJ datadir
AC_MSG_CHECKING([whether PROJ can be ran using the datadir])

AC_LANG_PUSH([C++])
LIBS="${PROJ_LIBS}"
AC_RUN_IFELSE([
    AC_LANG_PROGRAM([[
        #include <proj.h>
    ]], [[
        const char* paths[] = {"$PROJ_DATADIR"};
        proj_context_set_search_paths(NULL, 1, paths);
        proj_create(NULL, "EPSG:4326");
        int ret = proj_context_errno(NULL);
        return ret;
    ]])
], [
  AC_MSG_RESULT([yes])
], [
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([PROJ test program failed to run])
])
AC_LANG_POP([C++])

AC_SUBST([PROJ_DATADIR])

##################
## END
##################

# Create CPPFLAGS
CPPFLAGS="${PROJ_CFLAGS} ${GDAL_CFLAGS} ${CPPFLAGS} -I."

# Configure the output files
output_files="src/Makevars"
AS_IF([test "x$os_name" = "xWindows"],
  [
    CPPFLAGS="${CPPFLAGS} -DDLL_EXPORTS"    
    output_files="src/Makevars.ucrt"
    rm -f src/Makevars.ucrt.in
    cp src/Makevars.in src/Makevars.ucrt.in
  ]
)

AC_SUBST([CPPFLAGS])
AC_CONFIG_FILES([$output_files])
AC_OUTPUT