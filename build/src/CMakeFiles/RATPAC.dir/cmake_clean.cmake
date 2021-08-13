file(REMOVE_RECURSE
  "../lib/libRATPAC.a"
  "../lib/libRATPAC.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX Fortran)
  include(CMakeFiles/RATPAC.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
