#
# Search OR Download+Install CGL
#
include(ExternalProject)
ExternalProject_Add(CGAL_Download

    # Directory Options
    #PREFIX       "${CMAKE_BINARY_DIR}/third_party/CGAL"
    TMP_DIR      "${CMAKE_BINARY_DIR}/third_party/tmp/CGAL"
    STAMP_DIR    "${CMAKE_BINARY_DIR}/third_party/Stamp/CGAL"
    LOG_DIR      "${CMAKE_BINARY_DIR}/third_party/Stamp/CGAL"
    DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/third_party/Download/CGAL"
    SOURCE_DIR   "${CMAKE_BINARY_DIR}/third_party/Source/CGAL"
    BINARY_DIR   "${CMAKE_BINARY_DIR}/third_party/Build/CGAL"
    INSTALL_DIR  "${CMAKE_BINARY_DIR}/third_party/Install/CGAL"

    # Download Options
    GIT_REPOSITORY "https://github.com/CGAL/cgal"
    GIT_TAG "v5.4.1"
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE

    # Configure Step Options
    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=\"Release\"" 
               "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/third_party/Install/CGAL"

    # Build Step Options
    # - <default>

    # Install Step Options
    # - <default>
)

#
# Get Location if Installation
#
ExternalProject_Get_Property(CGAL_Download install_dir)

#
# Set Location to now find this stuff
# 
set(CGAL_FOUND       TRUE PARENT_SCOPE)
set(CGAL_ROOT        "${install_dir}" PARENT_SCOPE)
set(CGAL_INCLUDE_DIR "${install_dir}/include" PARENT_SCOPE)
# set(CGAL_LIB_DIR     "${install_dir}/lib64" PARENT_SCOPE)

# include(FindPackageHandleStandardArgs) 
# find_package_handle_standard_args(CGAL DEFAULT_MSG   
#     CGAL_ROOT
# 	CGAL_INCLUDE_DIR
#     CGAL_LIB_DIR
# )

# #
# # Export Modern Library
# #
# if(CGAL_FOUND AND NOT TARGET CGAL::CGAL)
#     add_library(CGAL::CGAL INTERFACE)
#     set_target_properties(CGAL::CGAL PROPERTIES
#         INTERFACE_INCLUDE_DIRECTORIES "${CGAL_INCLUDE_DIR}"
#         INTERFACE_LINK_LIBRARIES "${CGAL_LIB_DIR}"
#     )
# endif()

# #
# # Display Paths
# #
# if(CGAL_FOUND)
# 	if(NOT CGAL_FIND_QUIETLY)
#         message(STATUS "Future Locations of CGAL:")
#         message(STATUS "CGAL_ROOT        = ${CGAL_ROOT}")
#         message(STATUS "CGAL_INCLUDE_DIR = ${CGAL_INCLUDE_DIR}")
#         message(STATUS "CGAL_LIB_DIR     = ${CGAL_LIB_DIR}")
#     endif()
# endif()
