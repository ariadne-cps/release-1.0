cmake_minimum_required (VERSION 2.8.12)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Buddy REQUIRED)

add_executable(tutorial tutorial.cc)
target_link_libraries(tutorial ariadne ${Buddy_LIBRARIES})