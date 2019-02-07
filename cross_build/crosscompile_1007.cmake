# usage:  cmake -DCMAKE_INSTALL_PREFIX="C:/Users/gwe1si/AppData/Local/acados_install" -DCMAKE_TOOLCHAIN_FILE=crosscompile_1007.cmake -G "Unix Makefiles" ..


#cmake -DCMAKE_INSTALL_PREFIX="C:/Users/gwe1si/AppData/Local/acadosWindows_install" -DCMAKE_GENERATOR_PLATFORM=x64 ..
#call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\Tools\VsDevCmd.bat"
#call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\vcvars64.bat"
#msbuild ALL_BUILD.vcxproj
#msbuild INSTALL.vcxproj

# - when using cmake with Git Bash: -DCMAKE_SH="CMAKE_SH-NOTFOUND"
# - Add this file to acados/cmake
# - remove content from acados/cmake/Platform/dSpace.cmake 
#		- Issue: acados checks for variable "CMAKE_SYSTEM_NAME dSpace" to set generic dSpace settings
# 			- But: cmake requires a platform file dSpace.cmake which Robin adapted to the MicroAutoBox
# - Check which files make use of define _DS1007, DS_QNX_VERSION (compiler/1007, compiler/qnx) and what effects are
# Check if Robin's cmake file makes sense, given the dspace makefile for ds1104
# add_definitions(-D__DSPACE__) or custom __DSPACE__PPC; also, remove_definitions (see dSpace.cmake)
# custom implementation: timer.c/h, math.c/h fmax(), print.c printf, balsfeo i_aux_ext_dep_lib.c malloc_aligned()?
#	- check where these are implemented in ds1007 and ds1104

#set(CMAKE_SYSTEM_NAME Generic) # or 'Generic' if nothing special needs to be set/overwritten (like "-I" option)
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}) # look for custom platform file 'dSpace.cmake'

#QNX_HOST_TMP                  = $(CFG_COMPILER_ROOT_VARNAME)\host\win32\x86    # construct QNX_HOST relative to DSPACE_ROOT
#QNX_HOST_UNIX                 = $(QNX_HOST_TMP,S'\'/')                         # convert slashes into UNIX format
#%setenv QNX_HOST              = $(QNX_HOST_UNIX)                               # export QNX_HOST to DOS environment. It is needed by ntoppc-*.exe
#QNX_TARGET_TMP                = $(CFG_COMPILER_ROOT_VARNAME)\target\qnx6       # construct QNX_TARGET relative to DSPACE_ROOT
#QNX_TARGET_UNIX               = $(QNX_TARGET_TMP,S'\'/')                       # convert slashes into UNIX format
#%setenv QNX_TARGET            = $(QNX_TARGET_UNIX)                             # export QNX_HOST to DOS environment. It is needed by ntoppc-*.exe
#QNX_PATH                      = $(QNX_HOST_TMP)\usr\bin;$(DSPACE_ROOT)\exe     # add the QNX compiler's bin directory to the DOS search path
#%setenv PATH                  = $(QNX_PATH)                                    # (this is important because we need the cygwin1.dll)

##### Set compiler/linker on host system #####
# set(QNX_HOST "C:/ProgramData/dSPACE/071346EA-BFFA-4465-9551-2E48EDF35320/Compiler/QNX650_520/host/win32/x86")

# find_program(CMAKE_MAKE_PROGRAM  NAMES ${DSPACE_ROOT}/Exe/DSMAKE.exe)

set(CMAKE_MAKE_PROGRAM "C:/Program\ Files/dSPACE RCPHIL\ 2017-A/Exe/DSMAKE.exe")
set(CMAKE_C_COMPILER   "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-g++-5.2.0.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_CXX_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-g++-5.2.0.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_ASM_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-as.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_ASM_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-ar.exe") #  # ntoppc-gcc-5.2.0.exe

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x c" CACHE STRING "" FORCE) # needed when using ntoppcg++ as done in dspace makefile
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mpowerpc -mcpu=e500mc  -mcpu=e500mc -mtune=e500mc64 -mhard-float -mdouble-float -funsafe-math-optimizations -EB " CACHE STRING "" FORCE)
#set(CMAKE_RANLIB ":")  # tool to create index of archive
#set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <OBJECTS> -o <TARGET>") # syntax of archiver
##### Preprocessor directives
#set(ENV{QNX_HOST} "C:/ProgramData/dSPACE/071346EA-BFFA-4465-9551-2E48EDF35320/Compiler/QNX650_520/host/win32/x86")
#set(ENV{QNX_TARGET} "C:/ProgramData/dSPACE/071346EA-BFFA-4465-9551-2E48EDF35320/Compiler/QNX650_520/target/qnx6")
#set(ENV{PATH} "C:\\ProgramData\\dSPACE\\071346EA-BFFA-4465-9551-2E48EDF35320\\Compiler\\QNX650_520\\host\\usr\\bin;C:\\Program Files\\dSPACE RCPHIL 2018-A\\Exe")
add_definitions(-D_DS1007 -DDS_PLATFORM_PPC -DDS_PLATFORM_PPCE5500 -DDS_PLATFORM_BE -DDS_PLATFORM_QNX -DDS_PLATFORM_POSIX -DDS_PLATFORM_CN -DDS_PLATFORM_SMART -DDS_PLATFORM_SMARTRTK -DFEATURE_GIGALINK -DFEATURE_PHS_BUS -DFEATURE_ONBOARD_UART -DFEATURE_IO_ETHERNET -DFEATURE_DAQ -DFEATURE_MOTION_MANAGER) # CACHE STRING "" FORCE
#set(CMAKE_COMMAND "C://Program1 Files//CMake//bin//cmake.exe" CACHE STRING "" FORCE)

### Set make application ###

#find_program(CMAKE_MAKE_PROGRAM ${DSPACE_ROOT}/DS1007/WIN32/dsmake.exe)
#find_program(CMAKE_MAKE_PROGRAM ${DSPACE_ROOT}/exe/dsmake.exe)
#find_program(CMAKE_MAKE_PROGRAM NAMES DSMAKE.exe PATHS ${DSPACE_ROOT}/Exe/)
set(CMAKE_C_COMPILER_WORKS 1 CACHE INTERNAL "")
set(CMAKE_CXX_COMPILER_WORKS 1 CACHE INTERNAL "")

### Tell cmake where to find headers/libs of target system 

set(CMAKE_FIND_ROOT_PATH "C:/ProgramData/dSPACE/071346EA-BFFA-4465-9551-2E48EDF35320/Compiler/QNX650_520/target/qnx6")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

