cmake_minimum_required (VERSION 2.8.12)

project(Ariadne)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-Wno-pointer-bool-conversion")

find_package(Boost REQUIRED COMPONENTS system serialization thread)
find_package(Cairo REQUIRED)
find_package(Buddy REQUIRED)
include(UseDoxygen OPTIONAL)
include(cotire)

configure_file ("include/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h")
set(INCLUDE_TARGET "include")
configure_file ("include/ariadne.h.in"
        "${PROJECT_BINARY_DIR}/ariadne.h")
set(INCLUDE_TARGET "ariadne")
configure_file ("include/ariadne_install.h.in"
        "${PROJECT_BINARY_DIR}/ariadne_install.h")

include_directories ("include" ${PROJECT_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${CAIRO_INCLUDE_DIRS})

add_subdirectory(src)

add_library(ariadne SHARED $<TARGET_OBJECTS:ariadne-src>)

target_link_libraries(ariadne ${Boost_LIBRARIES} ${CAIRO_LIBRARIES} ${Buddy_LIBRARIES} --coverage)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
	target_compile_options(ariadne INTERFACE
			-O0        # no optimization
			-g         # generate debug info
			--coverage # sets all required flags
			)
endif()

cotire(ariadne CONFIGURATIONS None)

install (TARGETS ariadne DESTINATION lib)
install (FILES ${PROJECT_BINARY_DIR}/ariadne_install.h RENAME ariadne.h DESTINATION include)
install (FILES ${PROJECT_BINARY_DIR}/config.h DESTINATION include/ariadne)
file (GLOB library_includes "include/*.h")
install (FILES ${library_includes} DESTINATION include/ariadne)

enable_testing()

add_subdirectory(test)
add_subdirectory(examples)

