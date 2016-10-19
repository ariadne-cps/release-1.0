add_executable(painting painting.cc)
target_link_libraries(painting ariadne)

add_executable(painting-spatial painting-spatial.cc)
target_link_libraries(painting-spatial ariadne)

add_executable(painting-spatial-parametric painting-spatial-parametric.cc)
target_link_libraries(painting-spatial-parametric ariadne)