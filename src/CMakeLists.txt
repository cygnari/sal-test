add_library(sal-test
    fast_sum_funcs.cpp
    general_utils.cpp
    green_func.cpp
    icos_funcs.cpp
    interp_utils.cpp
    io_utils.cpp)

target_link_libraries(sal-test PUBLIC MPI::MPI_C LAPACK::LAPACK)
target_include_directories(sal-test PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

install(FILES
    # all FastBVE header/include files go here
    ${CMAKE_BINARY_DIR}/sal-test-config.h
    fast_sum_funcs.hpp
    general_utils.hpp
    green_func.hpp
    interp_utils.hpp
    io_utils.hpp
    icos_funcs.hpp
    structs.hpp
    DESTINATION ${PROJECT_SOURCE_DIR}/include)
