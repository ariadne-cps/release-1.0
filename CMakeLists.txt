cmake_minimum_required (VERSION 2.6)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Boost REQUIRED COMPONENTS system serialization thread)
find_package(Gtk2 REQUIRED)
find_package(Cairo REQUIRED)
find_package(Buddy REQUIRED)
find_package(Doxygen)

configure_file ("include/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h")

include_directories ("include" ${PROJECT_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${GTK2_INCLUDE_DIR} ${CAIRO_INCLUDE_DIRS})

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

enable_testing()

set(TESTS
	#bdd_set
	differential
	discretised_evolution
	expansion
	float
	function_set
	function
	graphics
	hybrid_evolution
	hybrid_io_automaton
	hybrid_set
	interval
	matrix
	misc
	multi_index
	polynomial
	reachability_analysis
	reachability_restriction
	taylor_function
	taylor_model
	taylor_set
	textplot
	vector
)

foreach(test ${TESTS})
    add_executable(test_${test} test/test_${test}.cc)
    target_link_libraries(test_${test} ariadne)
    add_test(test_${test} test_${test})
endforeach()

