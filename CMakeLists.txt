cmake_minimum_required(VERSION 3.1)
project(clion_carpio)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/)

set(SOURCE_FILES src/main.cpp)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(clion_carpio ${SOURCE_FILES})

if(VTK_LIBRARIES)
    target_link_libraries(clion_carpio ${VTK_LIBRARIES})
else()
    target_link_libraries(clion_carpio vtkHybrid vtkWidgets)
endif()