add_executable(laser laser.cc)
target_link_libraries(laser ariadne)
add_executable(laser-period laser-period.cc)
target_link_libraries(laser-period ariadne)

add_executable(laser-zi-check laser-zi-check.cc)
target_link_libraries(laser-zi-check ariadne)