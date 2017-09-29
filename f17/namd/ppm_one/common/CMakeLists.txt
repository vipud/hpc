cmake_minimum_required (VERSION 2.6)

project (common)

set (NMR_VERSION_MAJOR 1)
set (NMR_VERSION_MINOR 0)


configure_file (
  "${PROJECT_SOURCE_DIR}/config.h"
  "${PROJECT_BINARY_DIR}/config.h"
  )

SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/myprogram ) 
SET(CMAKE_PREFIX_PATH $ENV{HOME}/myprogram )

include_directories("${PROJECT_BINARY_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/common")
include_directories ("${PROJECT_SOURCE_DIR}/munkres")

set (EXTRA_LIBS ${EXTRA_LIBS} common xdrfile munkres )

add_subdirectory (common)
add_subdirectory (munkres)
