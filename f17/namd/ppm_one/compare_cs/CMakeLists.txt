cmake_minimum_required (VERSION 2.6)
project (compare_cs)



# add the executable
add_executable (compare_cs main.cpp )
target_link_libraries (compare_cs  ${EXTRA_LIBS})

# add the install targets
install (TARGETS compare_cs DESTINATION bin)


