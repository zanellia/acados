cmake_minimum_required(VERSION 2.8.11)
project(acados)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror -fPIC")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -fPIC")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message("\n\n --- The CMake build type is ${CMAKE_BUILD_TYPE} ---\n\n")

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fdiagnostics-show-option")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wno-unused-function -Wall -Wextra -Werror -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-unused-function -Wno-unused-parameter")
endif()

if(CMAKE_BUILD_TYPE MATCHES Profile)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use")
endif()

if(CMAKE_BUILD_TYPE MATCHES Test)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    enable_testing(true)
endif()


# Define include folders
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/external/")
include_directories("${PROJECT_SOURCE_DIR}/external/qpOASES/include/")
if(CMAKE_BUILD_TYPE MATCHES Test)
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/Eigen")
        message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but Eigen can not be found")
    endif()
    include_directories("${PROJECT_SOURCE_DIR}/external/eigen/")

    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.1.1")
        message(FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but casadi-octave-v3.1.1 can not be found")
    endif()
endif()

# external ocps
include_directories("${PROJECT_SOURCE_DIR}/experimental/andrea/modular_interface/ocps/")
include_directories("${PROJECT_SOURCE_DIR}/experimental/andrea/modular_interface/models/")
include_directories("${PROJECT_SOURCE_DIR}/experimental/andrea/modular_interface/autogen/")
include_directories("${PROJECT_SOURCE_DIR}/experimental/andrea/modular_interface/utils/")

# Define source files
include(acados_sources)
add_definitions(-DUSE_ACADOS_TYPES)

# Process external libraries
include(external/openblas)
include(external/blasfeo)
include(external/hpmpc)
#include(external/qpoases)
#include(external/qpdunes)
#include_directories("${PROJECT_SOURCE_DIR}/external/qpDUNES-dev/include/")

if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    add_definitions(-DOOQP)
    include(external/ma27)
    include(external/ooqp)
endif ()

# Define acados library
#add_library(acados ${ACADOS_SRC})
#add_dependencies(acados hpmpc blasfeo)

if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    add_dependencies(ooqpbase ma27)
    add_dependencies(acados ooqpbase)
endif ()

#add_subdirectory(swig)

# Define examples
add_library(kitty_acados SHARED ${ACADOS_CONTROLLER_MEM_SRC} ${ACADOS_SRC})

# force the compiler to add hpmpc to the shared library
SET (hpmpc -Wl,--whole-archive)
SET (blasfeo -Wl,--whole-archive)

#add_dependencises(kitty_acados hpmpc)
target_link_libraries(kitty_acados hpmpc blasfeo m)

set_target_properties(kitty_acados PROPERTIES VERSION 0.22)
set_target_properties(kitty_acados PROPERTIES SOVERSION 0.22)

# Automatic generation of unit testing data
if(CMAKE_BUILD_TYPE MATCHES Test)
    include(test_sources)
    add_custom_command(OUTPUT ${UNIT_TESTS_SRC_CASADI}
        COMMAND ${CMAKE_COMMAND} -P ../cmake/generate_test_data.cmake
        COMMENT "Generating unit testing files")
    add_custom_target(generate_test_data
        DEPENDS ${UNIT_TESTS_SRC_CASADI}
        COMMENT "Checking if regeneration of test data is needed")
    add_executable(unit_tests ${UNIT_TESTS_SRC} ${UNIT_TESTS_SRC_SIM} ${UNIT_TESTS_SRC_OCP_QP} ${UNIT_TESTS_SRC_OCP_NLP} ${UNIT_TESTS_SRC_TEST_UTILS})
    add_dependencies(unit_tests generate_test_data)

    target_link_libraries(unit_tests acados qpoases hpmpc blasfeo qpdunes m)

    if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
        find_library(FLIB gfortran HINTS /usr/lib/gcc/x86_64-linux-gnu/* /usr/local/lib/gcc/*)
        target_link_libraries(unit_tests ooqpgensparse ooqpsparse ooqpgondzio ooqpbase ma27 blas ${FLIB} m)
    endif ()

    set_target_properties(unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)
endif()

# use, i.e. don't skip the full RPATH for the build tree
#install(TARGETS acados acados_controller_mem DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

install(DIRECTORY acados DESTINATION include/kitty_acados)
install(DIRECTORY acados/utils/ DESTINATION include/kitty_acados/acados)
install(DIRECTORY acados/sim/ DESTINATION include/kitty_acados/acados)
install(DIRECTORY external/hpmpc/include DESTINATION include/kitty_acados/external/hpmpc)
install(DIRECTORY experimental/andrea/modular_interface/models DESTINATION include/kitty_acados/experimental/andrea/modular_interface)
install(DIRECTORY experimental/andrea/modular_interface/ocps DESTINATION include/kitty_acados/experimental/andrea/modular_interface)
install(DIRECTORY experimental/andrea/modular_interface/utils DESTINATION include/kitty_acados/experimental/andrea/modular_interface)

#install(external/hpmpc/libhpmpc.a DESTINATION share/lib)
#install(external/blasfeo/libblasfeo.a DESTINATION share/lib)

target_include_directories(kitty_acados PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/kitty_acados>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/acados/utils>
  $<INSTALL_INTERFACE:include/kitty_acados/acados/utils>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/acados/sim>
  $<INSTALL_INTERFACE:include/kitty_acados/acados/sim>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/hpmpc/include>
  $<INSTALL_INTERFACE:include/kitty_acados/external/hpmpc/include>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/experimental/andrea/modular_interface/models>
  $<INSTALL_INTERFACE:include/kitty_acados/experimental/andrea/modular_interface/models>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/experimental/andrea/modular_interface/ocps>
  $<INSTALL_INTERFACE:include/kitty_acados/experimental/andrea/modular_interface/ocps>

  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/experimental/andrea/modular_interface/utils>
  $<INSTALL_INTERFACE:include/kitty_acados/experimental/andrea/modular_interface/utils>
  )

install(TARGETS kitty_acados EXPORT KittyAcadosConfig
  ARCHIVE  DESTINATION lib
  LIBRARY  DESTINATION lib
  RUNTIME  DESTINATION bin)  # This is for Windows
#install(DIRECTORY include/ DESTINATION include)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT KittyAcadosConfig DESTINATION share/KittyAcados/cmake)

## This makes the project importable from the build directory
#export(TARGETS kitty_acados FILE KittyAcadosConfig.cmake)
