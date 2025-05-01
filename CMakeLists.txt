cmake_minimum_required(VERSION 3.15)
project(ray_tracer CXX)

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(OpenMP REQUIRED)

set(RT_CORE_SOURCES
    Sphere.cpp
    TriangleMesh.cpp
    stb_impl.cpp
)

add_library(rt_core STATIC ${RT_CORE_SOURCES})

target_include_directories(rt_core PUBLIC
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(rt_core PUBLIC OpenMP::OpenMP_CXX)

add_executable(spheres spheres.cpp)
target_link_libraries(spheres PRIVATE rt_core)

add_executable(cat cat.cpp)
target_link_libraries(cat PRIVATE rt_core)

install(TARGETS spheres cat
        RUNTIME DESTINATION bin)
