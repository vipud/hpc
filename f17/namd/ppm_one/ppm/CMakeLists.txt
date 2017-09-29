cmake_minimum_required (VERSION 2.6)
project (ppm_stable)


# add the executable
add_executable (ppm_stable main.cpp mainbody.cpp)
target_link_libraries (ppm_stable  ${EXTRA_LIBS})

# add the install targets
install (TARGETS ppm_stable DESTINATION bin)


