add_library(ariadne-src OBJECT
        numeric.cpp
        vector.cpp
        matrix.cpp
        linear_programming.cpp
        differential.cpp
        taylor_model.cpp
        expression.cpp
        affine.cpp
        polynomial.cpp
        function.cpp
        taylor_function.cpp
        simulation_toolbox.cpp
        taylor_calculus.cpp
        solver.cpp
        integrator.cpp
        point.cpp
        box.cpp
        zonotope.cpp
        polytope.cpp
        polyhedron.cpp
        curve.cpp
        function_set.cpp
        taylor_set.cpp
        grid.cpp
        denotable_set.cpp
        hybrid_set.cpp
        bdd_set.cpp
        geometry.cpp
        hybrid_automaton.cpp
        hybrid_io_automaton.cpp
        orbit.cpp
        hybrid_evolver-image.cpp
        set_checker.cpp
        reachability_restriction.cpp
        reachability_analyser.cpp
        logging.cpp
        graphics.cpp
        textplot.cpp
        parametric.cpp
        verifier.cpp
        verification_input.cpp
        )

cotire(ariadne-src CONFIGURATIONS None)

target_link_libraries(ariadne-src PUBLIC coverage_config)