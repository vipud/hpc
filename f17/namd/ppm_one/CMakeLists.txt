cmake_minimum_required (VERSION 2.6)
CMAKE_POLICY(SET CMP0012 NEW)
project (ppm)

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

set(USE_OPENMP ON CACHE BOOL "Enable openmp")

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -vec-report3")

if(${USE_OPENMP})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
    set (USE_OPENMP off)
endif(OPENMP_FOUND)
endif(${USE_OPENMP})


find_package (common)
IF(${LIBCOMMON_FOUND})
   MESSAGE(STATUS "Found common")
ELSE()
   MESSAGE(STATUS "Could not locate common.")
   MESSAGE("You can use cmake -DLIBCOMMON_DIR=common_main_dir to solve this")
ENDIF()

find_package (lmmin)
IF(${LIBLMMIN_FOUND})
   MESSAGE(STATUS "Found lmmin")
ELSE()
   MESSAGE(STATUS "Could not locate lmmin.")
   MESSAGE("You can use cmake -DLIBLMMIN_DIR=lmmin_main_dir to solve this")
ENDIF()

include_directories("${PROJECT_BINARY_DIR}")
include_directories(${LIBCOMMON_INCLUDE_DIRS})
include_directories(${LIBLMMIN_INCLUDE_DIRS})
set (EXTRA_LIBS ${EXTRA_LIBS} ${LIBCOMMON_LIBRARIES} ${LIBLMMIN_LIBRARIES}) 

add_subdirectory (ppm)
add_subdirectory (compare_cs)
