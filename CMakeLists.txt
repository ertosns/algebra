cmake_minimum_required(VERSION 3.2)

set(project_name algebra)

project(${project_name})


include_directories(include /usr/include/eigen3)

include_directories(include /usr/include/gtest)
find_library(gtest gtest)
set(test_name algebra_test)

find_library(pthread pthread)

#######################
# build libalgebra.a
#######################
add_library(${project_name} STATIC
  ./source/algebra.cpp
  )

################
# build tests
################

add_executable(${test_name} ./tests/algebra_test.cpp)
target_link_libraries(${test_name} algebra gtest pthread)

####################
# build docs
####################
execute_process(COMMAND doxygen
  WORKING_DIRECTORY docs
  )
execute_process(COMMAND make
  WORKING_DIRECTORY docs/latex
  )
