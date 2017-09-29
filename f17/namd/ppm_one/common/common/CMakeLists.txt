project(common)
add_library(common STATIC aa.cpp bmrb.cpp pdb.cpp  rdc.cpp  supply.cpp traj.cpp ) 
install (TARGETS common DESTINATION lib)
install (FILES aa.h  bmrb.h  pdb.h  rdc.h  supply.h  traj.h  DESTINATION include)

