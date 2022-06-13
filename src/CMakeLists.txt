#
#
#

#
# ========================================================
# make_random_sources
# ========================================================
#
add_executable(make_random_sources make_random_sources.cpp)
set_target_properties(make_random_sources       
PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
install(TARGETS make_random_sources DESTINATION "${MY_INSTALL_BINDIR}")

#
# ========================================================
# make_random_targets
# ========================================================
#
add_executable(make_random_targets make_random_targets.cpp)
set_target_properties(make_random_targets       
PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
install(TARGETS make_random_targets DESTINATION "${MY_INSTALL_BINDIR}")


#
# ========================================================
# cgil2
# ========================================================
#
add_executable(cgil2 cgil2.cpp)

set_target_properties(cgil2       
PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_compile_options(cgil2 
        PUBLIC
                "$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
                "-O3"
                "-DNDEBUG"
)

target_include_directories(cgil2 
   PUBLIC
      ${CGAL_INCLUDE_DIR}
)

target_link_libraries(cgil2 
        Boost::headers
)

install(TARGETS cgil2 DESTINATION "${MY_INSTALL_BINDIR}")

#
# ========================================================
# cgil3
# ========================================================
#
add_executable(cgil3 cgil3.cpp)

set_target_properties(cgil3     
PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_compile_options(cgil3 
        PUBLIC
                "$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
                "-O3"
                "-DNDEBUG"
)

target_include_directories(cgil3 
   PUBLIC
      ${CGAL_INCLUDE_DIR}
)

target_link_libraries(cgil3 
        Boost::headers
)

install(TARGETS cgil3 DESTINATION "${MY_INSTALL_BINDIR}")
