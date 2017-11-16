set(TEST_SUITE
	test_bdd_set
	test_differential
	test_discretised_evolution
	test_expansion
	test_float
	test_function_set
	test_function
	test_graphics
	test_hybrid_evolution
	test_hybrid_io_automaton
	test_hybrid_set
	test_interval
	test_matrix
	test_misc
	test_multi_index
	test_polynomial
	test_reachability_analysis
	test_reachability_restriction
	test_taylor_function
	test_taylor_model
	test_taylor_set
	test_integrator
	test_textplot
	test_vector
)

foreach(test ${TEST_SUITE})
    add_executable(${test} ${test}.cc)
    target_link_libraries(${test} ariadne)
    add_test(${test} ${test})
endforeach()

add_custom_target(tests DEPENDS ${TEST_SUITE})
