include(ExternalProject)

ExternalProject_Add(
    hpmpc_project
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/hpmpc"
    BINARY_DIR "${PROJECT_SOURCE_DIR}/external/hpmpc/build"
    CONFIGURE_COMMAND cmake ..
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(hpmpc_project source_dir)

add_library(hpmpc SHARED IMPORTED)
add_dependencies(hpmpc hpmpc_project blasfeo_project)
set_property(TARGET hpmpc PROPERTY IMPORTED_LOCATION "${source_dir}/build/libhpmpc.so")
