# - Try to find the buddy library
# Once done this will define
#
#  Buddy_FOUND - system has bdd
#  Buddy_INCLUDE_DIRS - the bdd include directory
#  Buddy_LIBRARIES - Link these to use bdd
#
# Define Buddy_MIN_VERSION for which version desired.
#

INCLUDE(FindPkgConfig)

IF(Buddy_FIND_REQUIRED)
	SET(_pkgconfig_REQUIRED "REQUIRED")
ELSE(Buddy_FIND_REQUIRED)
	SET(_pkgconfig_REQUIRED "")
ENDIF(Buddy_FIND_REQUIRED)

#IF(Buddy_MIN_VERSION)
#	PKG_SEARCH_MODULE(Buddy ${_pkgconfig_REQUIRED} bdd>=${Buddy_MIN_VERSION})
#ELSE(Buddy_MIN_VERSION)
#	PKG_SEARCH_MODULE(Buddy ${_pkgconfig_REQUIRED} bdd)
#ENDIF(Buddy_MIN_VERSION)

#IF(NOT Buddy_FOUND AND NOT PKG_CONFIG_FOUND)
	FIND_PATH(Buddy_INCLUDE_DIRS bdd.h)
	FIND_LIBRARY(Buddy_LIBRARIES bdd)

	# Report results
	IF(Buddy_LIBRARIES AND Buddy_INCLUDE_DIRS)
		SET(Buddy_FOUND 1)
		#IF(NOT Buddy_FIND_QUIETLY)
			MESSAGE(STATUS "Found Buddy: ${Buddy_LIBRARIES}")
		#ENDIF(NOT Buddy_FIND_QUIETLY)
	ELSE(Buddy_LIBRARIES AND Buddy_INCLUDE_DIRS)	
		IF(Buddy_FIND_REQUIRED)
			MESSAGE(SEND_ERROR "Could not find Buddy")
		ELSE(Buddy_FIND_REQUIRED)
			#IF(NOT Buddy_FIND_QUIETLY)
				MESSAGE(STATUS "Could not find Buddy")	
			#ENDIF(NOT Buddy_FIND_QUIETLY)
		ENDIF(Buddy_FIND_REQUIRED)
	ENDIF(Buddy_LIBRARIES AND Buddy_INCLUDE_DIRS)
#ENDIF(NOT Buddy_FOUND AND NOT PKG_CONFIG_FOUND)

# Hide advanced variables from CMake GUIs
MARK_AS_ADVANCED(Buddy_LIBRARIES Buddy_INCLUDE_DIRS)