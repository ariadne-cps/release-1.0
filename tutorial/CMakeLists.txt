# Preferable to avoid warnings during execution
cmake_minimum_required (VERSION 3.1)

# Set the executable along with the required source files
add_executable(tutorial tutorial.cpp)

# Instruct to link against the ariadne library, and the bdd library
target_link_libraries(tutorial ariadne bdd)