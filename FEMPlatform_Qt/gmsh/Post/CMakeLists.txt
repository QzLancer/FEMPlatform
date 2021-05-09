# message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SRC RELATIVE "*.cpp" "*.c" "*.cxx")
file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
append_gmsh_src(Post "${SRC};${HDR}")
