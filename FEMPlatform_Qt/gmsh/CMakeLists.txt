# PROJECT (gmsh)

# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# set(CMAKE_BUILD_TYPE "Release")


# if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wno-unused-variable -pthread")
# else( CMAKE_BUILD_TYPE STREQUAL "Debug" )
#         set(CMAKE_BUILD_TYPE "Release")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -pthread -fopenmp")
# endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
# message("Build Type:${CMAKE_BUILD_TYPE} ${CMAKE_CXX_FLAGS}")


macro(append_gmsh_src DIRNAME FILES)
foreach(FILE ${FILES})
  list(APPEND LIST ${DIRNAME}/${FILE})
endforeach()
set(GMSH_SRC ${GMSH_SRC};${LIST} PARENT_SCOPE)
set(GMSH_DIRS ${GMSH_DIRS};${DIRNAME} PARENT_SCOPE)
endmacro()

add_subdirectory(Common)
add_subdirectory(Numeric)
add_subdirectory(Geo)
add_subdirectory(Mesh)
add_subdirectory(Post)
add_subdirectory(Parser)

# message(${GMSH_DIRS})
include_directories(Common Geo Mesh Numeric Parser Post )

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/contrib/Netgen)
    add_subdirectory(contrib/Netgen)
    include_directories(contrib/Netgen contrib/Netgen/libsrc/include
                        contrib/Netgen/nglib)
    #set_config_option(HAVE_NETGEN "Netgen")
    #add_definitions(-DNO_PARALLEL_THREADS -DNOTCL)
endif()

add_library(gmsh STATIC ${GMSH_SRC})

