set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

cmake_minimum_required(VERSION 3.14.5)
project(main)
file(GLOB SOURCES "src/*.cpp")

add_executable(sim ${SOURCES})
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})
	
target_link_libraries(sim ${YAML_CPP_LIBRARIES})
