project(munkres)
add_library(munkres STATIC  munkres.cpp) 
install (TARGETS munkres DESTINATION lib)
install (FILES munkres.h matrix.h matrix.cpp  DESTINATION include)

