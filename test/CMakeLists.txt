set(TESTS
	bdd_set
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
    add_executable(test_${test} test_${test}.cc)
    target_link_libraries(test_${test} ariadne)
    add_test(test_${test} test_${test})
endforeach()