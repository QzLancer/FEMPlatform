# message(${CMAKE_CURRENT_SOURCE_DIR})

set(SRC
  Gmsh.tab.cpp
  Gmsh.yy.cpp
  FunctionManager.cpp
)
file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)


append_gmsh_src(Parser "${SRC};${HDR}")
