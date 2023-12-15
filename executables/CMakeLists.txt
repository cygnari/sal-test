add_executable(driver driver.cpp)
target_link_libraries(driver PUBLIC sal-test)
target_include_directories(driver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)
