cmake_minimum_required (VERSION 3.16)
project (assignment3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Debug")

set(${PROJECT_NAME}_VERSION_MAJOR 1)
set(${PROJECT_NAME}_VERSION_MINOR 0)

include_directories( ${CMAKE_SOURCE_DIR}/include/ )

# CGAL
find_package( CGAL QUIET COMPONENTS )
if ( CGAL_FOUND )
  message(STATUS "CGAL found")
  message(STATUS ${CGAL_LIBRARIES})
  message(STATUS ${CGAL_3RD_PARTY_LIBRARIES})
else()
  message(SEND_ERROR "CGAL library is required")
  return()  
endif()

find_package(Eigen3 REQUIRED)
include(CGAL_Eigen3_support)

add_executable(assignment3 src/main.cpp
        include/structs.h)

target_compile_options(assignment3 PRIVATE -Wall)

target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::Eigen3_support)

if(TBB_FOUND)
  target_link_libraries(${PROJECT_NAME} TBB::tbb)
  target_compile_definitions(${PROJECT_NAME} PRIVATE CGAL_LINKED_WITH_TBB)
endif()