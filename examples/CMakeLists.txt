set(EXAMPLES
		 watertank-monolithic-hysteresis 
		 watertank-compositional-hysteresis 
	 	 watertank-proportional 
		 watertank-nonlinear-monolithic-hysteresis 
		 watertank-nonlinear-monolithic-proportional
		 watertank-monolithic-dominance 
		 watertank-nonlinear-monolithic-dominance 
		 bouncingball
		 springs
		 boost-converter
		 home
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(experimental)
