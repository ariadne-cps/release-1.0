set(EXPERIMENTAL
		 watertank-aasap
		 rectifier
		 cmos_inverter_analoginput 
		 robotarm
		 fb_refinement 
		 traffic-arbitrator4 
		 exponential
		 circle
		 circle-unstable
)

foreach(experiment ${EXPERIMENTAL})
    add_executable(${experiment} ${experiment}.cc)
    target_link_libraries(${experiment} ariadne)
endforeach()
