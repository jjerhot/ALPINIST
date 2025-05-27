cmake_minimum_required(VERSION 3.5)
project(DecayMC)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

### Add ROOT
find_package(ROOT REQUIRED EG)
if( ${ROOT_FOUND} )
    include(${ROOT_USE_FILE})
    message (STATUS "Found ROOT ${ROOT_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
else()
    message (FATAL_ERROR "ROOT not found")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# root_generate_dictionary(ExpParameters_Dict ExpParameters.h LINKDEF include/ExpParameters_LinkDef.h)

# add_library(ExpParameters SHARED src/ExpParameters.C ExpParameters_Dict.cxx)
# set_target_properties(ExpParameters PROPERTIES PUBLIC_HEADER include/ExpParameters.h)
# set_target_properties(ExpParameters PROPERTIES LINKER_LANGUAGE CXX)
# target_link_libraries(ExpParameters ${ROOT_LIBRARIES})

add_executable(DecayMC src/DecayMC.C src/ExpParameters.C src/Particle.C)
# target_link_libraries(DecayMC ExpParameters ${ROOT_LIBRARIES})
target_link_libraries(DecayMC ${ROOT_LIBRARIES})