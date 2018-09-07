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

include_directories ("include" ${PROJECT_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${GTK2_INCLUDE_DIRS} ${CAIRO_INCLUDE_DIRS})

add_library (ariadne SHARED "src/numeric.cc" "src/vector.cc" "src/matrix.cc" "src/linear_programming.cc" 
							"src/differential.cc" "src/taylor_model.cc" "src/expression.cc" "src/affine.cc" "src/polynomial.cc" 
							"src/function.cc" "src/taylor_function.cc" "src/simulation_toolbox.cc" "src/taylor_calculus.cc"
							"src/solver.cc" "src/integrator.cc" "src/point.cc" "src/box.cc" "src/zonotope.cc" "src/polytope.cc" 	
							"src/polyhedron.cc" "src/curve.cc" "src/function_set.cc" "src/taylor_set.cc" "src/grid.cc" 
							"src/denotable_set.cc" "src/hybrid_set.cc" "src/bdd_set.cc" "src/geometry.cc"
							"src/hybrid_automaton.cc" "src/hybrid_io_automaton.cc" "src/orbit.cc" "src/hybrid_evolver-image.cc" 
							"src/set_checker.cc" "src/reachability_restriction.cc" "src/reachability_analyser.cc" "src/logging.cc" 
							"src/graphics.cc" "src/textplot.cc" "src/parametric.cc" "src/verifier.cc" "src/verification_input.cc")

target_link_libraries(ariadne ${Boost_LIBRARIES} ${GTK2_LIBRARIES} ${CAIRO_LIBRARIES} ${Buddy_LIBRARIES})

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    include(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    set(COVERAGE_LCOV_EXCLUDES '/usr/*' 'examples/*' 'test/*')
    SETUP_TARGET_FOR_COVERAGE_LCOV(NAME lcov_coverage EXECUTABLE ctest DEPENDENCIES tests)
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

