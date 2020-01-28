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

# if everything breaks with the testing, just comment everything below
# add packages in for catch2 unit testing
find_package(Catch2 REQUIRED)
include_directories(${CATCH2_INCLUDE_DIR})

add_executable(test_sim 
	catch_testing/main_config.cpp 
	catch_testing/properties_test.hpp
	src/Particle.cpp
	src/Properties.cpp
	src/Parameters.cpp)

target_link_libraries(test_sim Catch2::Catch2 ${YAML_CPP_LIBRARIES})

include(CTest)
include(Catch)
catch_discover_tests(test_sim)
