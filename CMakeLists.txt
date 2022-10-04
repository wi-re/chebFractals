cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
# CMake policy settings
cmake_policy(VERSION 3.9)												# required for boost 1.65.1

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/cmake") # add local CMake files 
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON) 						
SET(CMAKE_INCLUDE_CURRENT_DIR ON)	# needed for QT 
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
     
# Project settings
SET(NAME "Chebyshev Fractals")													# name of the Project, variable for convenience
project(${NAME} LANGUAGES CXX) 

FIND_PACKAGE( Boost 		REQUIRED COMPONENTS iostreams system thread program_options filesystem)
FIND_PACKAGE( yaml-cpp		REQUIRED )
FIND_PACKAGE( glfw3			REQUIRED )
FIND_PACKAGE( glad			REQUIRED ) 
FIND_PACKAGE( LAPACK 		REQUIRED ) 
FIND_PACKAGE( Armadillo 	REQUIRED ) 


# Configure file to store the binary and source path for resolving file paths
configure_file(config.h.in config/config.h) 
# Set Include paths, due to everything needing almost everything this is done here
SET(INCLUDE_DIRS
	${Boost_INCLUDE_DIR} 		
	${CMAKE_SOURCE_DIR} 
	${CMAKE_SOURCE_DIR} 		
	${CMAKE_BINARY_DIR}   
	${YAML_CPP_INCLUDE_DIR}
	${GLAD_INCLUDE_DIR}
	${GLFW3_INCLUDE_DIR}
	${ARMADILLO_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})   
  
IF(WIN32)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /arch:AVX2 /openmp /permissive-")
ELSE()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp -march=native")
ENDIF()
add_compile_definitions(NOMINMAX _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)


LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} )			
   
# Macro useful to gather all sources for cuda, c++ and qt code and sort them in the correct source group
macro(_GLOB_SOURCES LIBRARY_NAME)
	file(GLOB_RECURSE CPP_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.cpp)
	file(GLOB_RECURSE CPP_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.h ${LIBRARY_NAME}/*.hpp)
	SET(source_list  ${CPP_SOURCES} ${CPP_HEADERS})
	   
	foreach(source IN LISTS source_list)
		get_filename_component(source_path "${source}" PATH) 
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()  
	SET(source_list ${CPP_SOURCES} ${CPP_HEADERS})
endmacro(_GLOB_SOURCES) 
     
# Macro used to link all system libraries, e.g.  boost and TBB but NOT cuda
macro(_SETUP_DEPENDENCIES MODULE_NAME)
	target_link_libraries(${MODULE_NAME}
  		${Boost_SYSTEM_LIBRARY_RELEASE}
  		${Boost_THREAD_LIBRARY_RELEASE} 
		${Boost_LIBRARIES}  
		Boost::filesystem
		Boost::program_options
		${BLOSC_blosc_LIBRARY}
		yaml-cpp
		${ARGN}
		${YAML_CPP_LIBRARIES}
		${ARMADILLO_LIBRARIES}
		${LAPACK_LIBRARIES}
		Boost::boost
	) 
	IF(WIN32)
	ELSE()	
	target_link_libraries(${MODULE_NAME} stdc++)
	ENDIF()
	target_include_directories(${MODULE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})
endmacro(_SETUP_DEPENDENCIES)  
  
function(buildLibrary LIBRARY_NAME TARGETLIST PREFIX ) 
	_GLOB_SOURCES("${PREFIX}/${LIBRARY_NAME}")
	add_library( ${LIBRARY_NAME} ${CPP_SOURCES} ${CPP_HEADERS}) 
	_SETUP_DEPENDENCIES(${LIBRARY_NAME} ${ARGN})
endfunction(buildLibrary) 

_GLOB_SOURCES("cheb")
SET(chebSources "${source_list}")
_GLOB_SOURCES("gui")
SET(guiSources "${source_list}")
_GLOB_SOURCES("imgui")
SET(imguiSources "${source_list}")
_GLOB_SOURCES("fsVis")
SET(visSources "${source_list}")
_GLOB_SOURCES("tools")
SET(toolSources "${source_list}")


add_executable(online Source.cpp ${chebSources} ${guiSources} ${imguiSources} ${visSources} ${toolSources})

add_executable(offline offline.cpp ${chebSources} ${guiSources} ${imguiSources} ${visSources} ${toolSources})

add_executable(viewer viewer.cpp ${chebSources} ${guiSources} ${imguiSources} ${visSources} ${toolSources})

_SETUP_DEPENDENCIES(online glfw glad::glad ) 
_SETUP_DEPENDENCIES(viewer glfw glad::glad ) 
_SETUP_DEPENDENCIES(offline glfw glad::glad ) 
 