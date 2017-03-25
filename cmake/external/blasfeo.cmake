include(ExternalProject)

if(NOT DEFINED BLASFEO_TARGET)
#   set(BLASFEO_TARGET X64_INTEL_HASWELL)
#   set(BLASFEO_TARGET X64_INTEL_SANDY_BRIDGE)
    set(BLASFEO_TARGET GENERIC)
endif()
if(NOT DEFINED BLASFEO_LA)
#   set(BLASFEO_LA HIGH_PERFORMANCE)
    set(BLASFEO_LA REFERENCE)
#   set(BLASFEO_LA BLAS)
endif()

ExternalProject_Add(
    blasfeo_project
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make clean shared_library -j 2 TARGET=${BLASFEO_TARGET} LA=${BLASFEO_LA}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    STEP_TARGETS build install
    EXCLUDE_FROM_ALL TRUE
	#	INSTALL_COMMAND make install PREFIX=${CMAKE_INSTALL_PREFIX}
  # INSTALL_COMMAND ""
)

set_target_properties(blasfeo_project PROPERTIES VERSION 0.22)
set_target_properties(blasfeo_project PROPERTIES SOVERSION 0.22)
ExternalProject_Get_Property(blasfeo_project source_dir)
add_library(blasfeo SHARED IMPORTED)
add_dependencies(blasfeo blasfeo_project-build)
set_property(TARGET blasfeo PROPERTY IMPORTED_LOCATION "${source_dir}/libblasfeo.so")
