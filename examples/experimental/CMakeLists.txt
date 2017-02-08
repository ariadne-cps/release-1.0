set(EXPERIMENTAL
		 watertank-aasap
		 rectifier
		 cmos_inverter_analoginput 
		 robotarm
		 fb_refinement 
		 traffic-arbitrator4 
		 exponential
		 double-exponential
		 circle
		 range-check
		 sin
		 composition
         vanderpol
)

foreach(experiment ${EXPERIMENTAL})
    add_executable(${experiment} ${experiment}.cc)
    target_link_libraries(${experiment} ariadne)
endforeach()

add_subdirectory(laser)
add_subdirectory(painting)
