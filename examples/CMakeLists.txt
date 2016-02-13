set(EXAMPLES
		 watertank-monolithic-hysteresis 
		 watertank-compositional-hysteresis 
	 	 watertank-monolithic-proportional 
		 watertank-nonlinear-monolithic-hysteresis 
		 watertank-nonlinear-monolithic-proportional 
		 watertank-aasap 
		 watertank-monolithic-dominance 
		 watertank-nonlinear-monolithic-dominance 
		 bouncingball 
		 3bouncingballs 
		 ecoli 
		 sinusoid-sin 
		 rectifier 
		 springs
		 cmos_inverter_analoginput 
		 robotarm
		 fb_refinement 
		 traffic-arbitrator4 
		 boost-converter
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()