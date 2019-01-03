cmake_minimum_required (VERSION 2.8.11)
project (vector-commitment)
set(NAME "vector-commitment")

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} cmake/)
include_directories(${CMAKE_SOURCE_DIR})
include_directories("./ate-pairing/include/")
include_directories("./ate-pairing/test/")
link_directories("./ate-pairing/lib/")




set (CMAKE_C_FLAGS "-pthread  -march=native -O3 -maes")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -ggdb")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11")


add_executable(test test.cpp vcs.cpp)
target_link_libraries(test gmp zm gmpxx)
