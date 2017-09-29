
find_path(LIBLMMIN_INCLUDE_DIR lmmin.h 
          HINTS "${LIBLMMIN_DIR}"
	  PATH_SUFFIXES include )

find_library(LIBLMMIN_LIBRARY NAMES liblmmin.a
          HINTS "${LIBLMMIN_DIR}"
	  PATH_SUFFIXES lib  )

set(LIBLMMIN_LIBRARIES ${LIBLMMIN_LIBRARY} )
set(LIBLMMIN_INCLUDE_DIRS ${LIBLMMIN_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBLMMIN DEFAULT_MSG LIBLMMIN_LIBRARY LIBLMMIN_INCLUDE_DIR)
