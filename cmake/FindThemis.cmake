# Author : Simon Kallweit simon.kallweit@gmail.com

# This module will define the following variables:
#  Themis_INCLUDE_DIRS - Location of the ilmbase includes
#  Themis_LIBRARIES - [TODO] Required libraries for all requested bindings
#  Themis_FOUND - true if ILMBASE was found on the system
#  Themis_LIBRARY_DIRS - the full set of library directories

# find_path(Themis_Base_Dir include/themis/themis.h ENV Themis_ROOT)

# if(Themis_Base_Dir)

#   set(Themis_INCLUDE_DIRS ${Themis_Base_Dir}/include CACHE STRING "Themis include directories")
#   set(Themis_LIBRARY_DIRS ${Themis_Base_Dir}/lib CACHE STRING "Themis library directories")
#   set(Themis_FOUND TRUE)

# endif(Themis_Base_Dir)

find_path(Themis_INCLUDE_DIR NAMES themis/themis.h)
find_library(Themis_LIBRARY NAMES themis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Themis REQUIRED_VARS Themis_INCLUDE_DIR Themis_LIBRARY)

set(Themis_LIBRARIES ${Themis_LIBRARY})
set(Themis_INCLUDE_DIRS ${Themis_INCLUDE_DIR})

mark_as_advanced(Themis_INCLUDE_DIR Themis_LIBRARY)