cmake_minimum_required (VERSION 2.8.12)
if(POLICY CMP0045)
  cmake_policy(SET CMP0045 OLD)
endif()

if(POLICY CMP0058)
  cmake_policy(SET CMP0058 NEW)
endif()

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-Wno-pointer-bool-conversion")

find_package(Boost REQUIRED COMPONENTS system serialization thread)
find_package(GTK2 REQUIRED)
find_package(Cairo REQUIRED)
find_package(Buddy REQUIRED)

add_executable(tutorial tutorial.cc)
target_link_libraries(tutorial ariadne ${Boost_LIBRARIES} ${GTK2_LIBRARIES} ${CAIRO_LIBRARIES} ${Buddy_LIBRARIES})