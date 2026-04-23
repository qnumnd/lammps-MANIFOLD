if(PKG_MANIFOLD)
  find_package(Boost REQUIRED)
  target_link_libraries(lammps PRIVATE Boost::headers)
endif()