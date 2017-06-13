set(EXAMPLES
		 watertank-hysteresis 
	 	 watertank-proportional
		 watertank-dominance
		 bouncingball
		 springs
		 vanderpol
		 boost-converter
		 home
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(experimental)
