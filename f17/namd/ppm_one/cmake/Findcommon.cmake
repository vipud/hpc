#set( LIBCOMMON_DIR "/panfs/storage.local/bruschweiler/home/dli2/xdr-files" )

find_path(LIBCOMMON_INCLUDE_DIR aa.h 
          HINTS "${LIBCOMMON_DIR}"
	  PATH_SUFFIXES include )

find_library(LIBCOMMON_LIBRARY NAMES libcommon.a
          HINTS "${LIBCOMMON_DIR}"
	  PATH_SUFFIXES lib  )

set(LIBCOMMON_LIBRARIES ${LIBCOMMON_LIBRARY} )
set(LIBCOMMON_INCLUDE_DIRS ${LIBCOMMON_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBCOMMON DEFAULT_MSG LIBCOMMON_LIBRARY LIBCOMMON_INCLUDE_DIR)
