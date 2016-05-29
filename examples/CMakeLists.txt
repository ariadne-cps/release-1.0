set(EXAMPLES
		 watertank-monolithic-hysteresis 
		 watertank-compositional-hysteresis 
	 	 watertank-proportional 
		 watertank-nonlinear-monolithic-hysteresis 
		 watertank-nonlinear-monolithic-proportional 
		 watertank-aasap 
		 watertank-monolithic-dominance 
		 watertank-nonlinear-monolithic-dominance 
		 bouncingball 
		 rectifier 
		 springs
		 cmos_inverter_analoginput 
		 robotarm
		 fb_refinement 
		 traffic-arbitrator4 
		 boost-converter
		 exponential
		 circle
		 thermostat
)

foreach(example ${EXAMPLES})
    add_executable(${example} ${example}.cc)
    target_link_libraries(${example} ariadne)
endforeach()
