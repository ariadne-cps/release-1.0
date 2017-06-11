cmake_minimum_required (VERSION 2.8.12)

add_executable(tutorial tutorial.cc)
target_link_libraries(tutorial ariadne bdd)