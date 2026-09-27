cmake_minimum_required( VERSION 3.16 FATAL_ERROR )

project( criticalFlocking )

# Location of the Cinder checkout. Defaults to a "Cinder" directory next to this repo
# (which is where scripts/setup_cinder.sh puts it). Override with -DCINDER_PATH=/path/to/Cinder
set( CINDER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../Cinder" CACHE PATH "Path to the root of the Cinder checkout" )
get_filename_component( CINDER_PATH "${CINDER_PATH}" ABSOLUTE )

if( NOT EXISTS "${CINDER_PATH}/proj/cmake/modules/cinderMakeApp.cmake" )
	message( FATAL_ERROR "Cinder not found at '${CINDER_PATH}'. Run scripts/setup_cinder.sh or pass -DCINDER_PATH=..." )
endif()

# Cinder's configure step defaults to Debug; default to Release instead since the simulation is O(N^2)
if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug or Release)" FORCE )
endif()

include( "${CINDER_PATH}/proj/cmake/modules/cinderMakeApp.cmake" )

# Simulation core (no Cinder dependency), shared by the app and the headless tool
add_library( flocksim STATIC
	${CMAKE_CURRENT_SOURCE_DIR}/src/sim/Flock.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sim/Measure.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sim/Params.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/sim/Predator.cpp
)
target_include_directories( flocksim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src )
target_compile_features( flocksim PUBLIC cxx_std_17 )

ci_make_app(
	APP_NAME    criticalFlocking
	CINDER_PATH ${CINDER_PATH}
	SOURCES     ${CMAKE_CURRENT_SOURCE_DIR}/src/FlockingApp.cpp
	LIBRARIES   flocksim
)

# Headless driver for calibration and validation runs (see tools/flocksim_cli.cpp)
add_executable( flocksim_cli ${CMAKE_CURRENT_SOURCE_DIR}/tools/flocksim_cli.cpp )
target_link_libraries( flocksim_cli PRIVATE flocksim )
