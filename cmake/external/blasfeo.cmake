include(ExternalProject)

ExternalProject_Add(
    blasfeo_project
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo"
    BINARY_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo/build"
    CONFIGURE_COMMAND cmake ..
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(blasfeo_project source_dir)

add_library(blasfeo SHARED IMPORTED)
add_dependencies(blasfeo blasfeo_project)
set_property(TARGET blasfeo PROPERTY IMPORTED_LOCATION "${source_dir}/build/libblasfeo.so")
