cmake_minimum_required (VERSION 2.6)
CMAKE_POLICY(SET CMP0012 NEW)

project (lmmin)


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) 

set (NMR_VERSION_MAJOR 1)
set (NMR_VERSION_MINOR 0)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h"
  "${PROJECT_BINARY_DIR}/config.h"
  )

SET(BUILD_SHARED_LIBS  false)
SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/myprogram ) 
SET(CMAKE_PREFIX_PATH $ENV{HOME}/myprogram )

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -vec-report3")

set(USE_OPENMP ON CACHE BOOL "Enable openmp")
if(${USE_OPENMP})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
    set (USE_OPENMP off)
endif(OPENMP_FOUND)
endif(${USE_OPENMP})


include_directories("${PROJECT_BINARY_DIR}")

add_library(lmmin STATIC lmmin.c lmcurve.c ann.cpp ) 
install (TARGETS lmmin DESTINATION lib)
install (FILES  lmstruct.h lmcurve.h lmmin.h ann.h DESTINATION include)

