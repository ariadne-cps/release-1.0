add_executable(laser-generic laser.cc)
target_link_libraries(laser-generic ariadne)
add_executable(laser-zi-check laser-zi-check.cc)
target_link_libraries(laser-zi-check ariadne)