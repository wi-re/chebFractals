cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
# CMake policy settings
cmake_policy(VERSION 3.9)												# required for boost 1.65.1
cmake_policy(SET CMP0071 NEW)											# automoc setting policy introduced in CMake 3.10
   
# CMake settings to build cuda projects properly and to work well with other depdencies, e.g. boost
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
SET(NAME DemoProject)													# name of the Project, variable for convenience
project(${NAME} LANGUAGES CXX CUDA) 
enable_language(CUDA)
#find_package(CUDAToolkit REQUIRED)

# steps to determine the GPU architecture of the system to only build one architecture
SET(OUTPUTFILE ${CMAKE_CURRENT_BINARY_DIR}/cuda_script) # No suffix required
SET(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_cuda.cu)
execute_process(COMMAND nvcc -lcuda ${CUDAFILE} -o ${OUTPUTFILE})
execute_process(COMMAND nvcc --version)
execute_process(COMMAND ${OUTPUTFILE} RESULT_VARIABLE CUDA_RETURN_CODE OUTPUT_VARIABLE ARCH)
if(${CUDA_RETURN_CODE} EQUAL 0)
	SET(CUDA_SUCCESS "TRUE")
else() 
	SET(CUDA_SUCCESS "FALSE")
endif() 
if (${CUDA_SUCCESS})
	message(STATUS "CUDA Architecture: ${ARCH}")
else()
    message(WARNING ${ARCH})
endif()
  
SET(ARCH 75)
SET(CUDA_ARCHITECTURES 75)

add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

FIND_PACKAGE( Boost 		REQUIRED COMPONENTS iostreams system thread program_options filesystem)
#FIND_PACKAGE( Bullet		REQUIRED )											# Rigid body physics
#FIND_PACKAGE( ALEMBIC 		REQUIRED )									# Required for particle based I/O
FIND_PACKAGE( BLOSC			REQUIRED )									
FIND_PACKAGE( ZLIB			REQUIRED )									
# Package finding for all components required should be build as per the readme.md file
FIND_PACKAGE( yaml-cpp		REQUIRED )
FIND_PACKAGE( glfw3			REQUIRED )
FIND_PACKAGE( glad			REQUIRED ) 
FIND_PACKAGE( LAPACK 			REQUIRED ) 
FIND_PACKAGE( OpenBLAS 			REQUIRED ) 
FIND_PACKAGE( Armadillo 			REQUIRED ) 
FIND_PACKAGE(FFTW3 REQUIRED)
find_package(benchmark REQUIRED)
find_package(GSL REQUIRED)
find_library(GMP_LIBRARY NAMES gmp mpir)
#find_package(mpfr REQUIRED)


# Configure file to store the binary and source path for resolving file paths
configure_file(cmake/config.h.in config/config.h) 
# Set Include paths, due to everything needing almost everything this is done here
SET(INCLUDE_DIRS
	#${ALEMBIC_INCLUDE_DIR} 
	#${ALEMBIC_INCLUDE_DIRS} 
	#${OPENEXR_INCLUDE_DIR} 	
	#${OPENEXR_INCLUDE_DIR}/OpenEXR
	${Boost_INCLUDE_DIR} 		
	#${ILMBASE_INCLUDE_DIR} 	
	#${ILMBASE_INCLUDE_DIR}/OpenEXR
	${CMAKE_SOURCE_DIR} 	
	${TBB_INCLUDE_DIR} 
	${ZLIB_INCLUDE_DIR} 
	${BLOSC_INCLUDE_DIR}  
	${CMAKE_SOURCE_DIR} 		
	${CMAKE_BINARY_DIR}   
	${YAML_CPP_INCLUDE_DIR}
	${GLAD_INCLUDE_DIR}
	${GLFW3_INCLUDE_DIR}
	${ARMADILLO_INCLUDE_DIR}
	${FFTW_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})   
  
SET(CUDA_SEPARABLE_COMPILATION ON)								 		
SET(CUDA_HOST_COMPILATION_CPP ON)	
SET(CUDA_ATTACH_VS_BUILD_RULE_TO_FILE ON) 	
INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
LINK_DIRECTORIES(${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
LINK_DIRECTORIES(${CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES})
#CMAKE is broken with respect to CMAKE_CXX_STANDARD and clang-tidy as it generates -std=c++17, which is unknown, instead of /std:c++17
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp")
add_compile_definitions(NOMINMAX _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
SET(CMAKE_CUDA_FLAGS_RELEASE 
"${CMAKE_CUDA_FLAGS} -O3 --std c++17 
-gencode=arch=compute_${ARCH},code=sm_${ARCH} -gencode=arch=compute_${ARCH},code=compute_${ARCH} 
--expt-extended-lambda --expt-relaxed-constexpr -ftz=true -keep -use_fast_math -src-in-ptx -lineinfo -Xptxas \"-v -warn-double-usage\"")
SET(CMAKE_CUDA_FLAGS_RELWITHDEBINFO
"${CMAKE_CUDA_FLAGS} -O3 --std c++17 
-gencode=arch=compute_${ARCH},code=sm_${ARCH} -gencode=arch=compute_${ARCH},code=compute_${ARCH} 
--expt-extended-lambda --expt-relaxed-constexpr -ftz=true -keep -use_fast_math -src-in-ptx -lineinfo -Xptxas \"-v -warn-double-usage\"")
SET(CMAKE_CUDA_FLAGS_DEBUG
"${CMAKE_CUDA_FLAGS} -g --std c++17 
-gencode=arch=compute_${ARCH},code=sm_${ARCH} -gencode=arch=compute_${ARCH},code=compute_${ARCH} 
--expt-extended-lambda --expt-relaxed-constexpr -ftz=true -keep -use_fast_math -src-in-ptx -lineinfo -Xptxas \"-warn-double-usage\"")
#SET(CMAKE_CUDA_FLAGS "${CUDA_FLAGS_FOR_RELEASE}")
 
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} )					 		
LINK_DIRECTORIES(${TBB_LIBRARY_DIRS} )
   
# Macro useful to gather all sources for cuda, c++ and qt code and sort them in the correct source group
macro(_GLOB_SOURCES LIBRARY_NAME)
	file(GLOB_RECURSE CUDA_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.cu)
	file(GLOB_RECURSE CUDA_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.cuh)
	file(GLOB_RECURSE CPP_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.cpp)
	file(GLOB_RECURSE CPP_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.h ${LIBRARY_NAME}/*.hpp)
	SET(source_list ${CUDA_SOURCES} ${CUDA_HEADERS} ${CPP_SOURCES} ${CPP_HEADERS})
	   
	foreach(source IN LISTS source_list)
		get_filename_component(source_path "${source}" PATH) 
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()  
	SET(source_list ${CUDA_SOURCES} ${CUDA_HEADERS} ${CPP_SOURCES} ${CPP_HEADERS})
	# MESSAGE(WARNING ${CUDA_SOURCES} ${CUDA_HEADERS} ${CPP_SOURCES} ${CPP_HEADERS})
endmacro(_GLOB_SOURCES) 
     
# Macro used to link all system libraries, e.g.  boost and TBB but NOT cuda
macro(_SETUP_DEPENDENCIES MODULE_NAME)
	target_link_libraries(${MODULE_NAME}
		#Alembic::Alembic
		#${OPENEXR_ILMIMF_LIBRARY_RELEASE}
		#${ALEMBIC_ILMBASE_IMATH_LIB} 
		#${ALEMBIC_ILMBASE_ILMTHREAD_LIB}
		#${ALEMBIC_ILMBASE_IEX_LIB}
		#${ALEMBIC_ILMBASE_IEXMATH_LIB}
		#${ALEMBIC_ILMBASE_HALF_LIB}
		#${ILMBASE_HALF_LIBRARY_RELEASE}  
  		${Boost_SYSTEM_LIBRARY_RELEASE}
  		${Boost_THREAD_LIBRARY_RELEASE} 
		${Boost_LIBRARIES}  
		Boost::filesystem
		Boost::program_options
		${GVDB_LIB_DIR}/${GVDB_LIB}
		${BLOSC_blosc_LIBRARY}
		${TBB_LIBRARY}   
		yaml-cpp
		${ARGN}
		${YAML_CPP_LIBRARIES}
		${ARMADILLO_LIBRARIES}
		${LAPACK_LIBRARIES}
		${BLAS_openblas_LIBRARY}
		FFTW3::fftw3
		benchmark::benchmark
		GSL::gsl GSL::gslcblas
		Boost::boost
		#C:/vcpkg/installed/x64-windows/lib/mpfr.lib

	)  
	target_link_libraries(${MODULE_NAME}
		${CUDA_CUDART_LIBRARY}
		${CUDA_LIBRARIES}
		${CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES}
		${CMAKE_CUDA_HOST_IMPLICIT_LINK_LIBRARIES}	
	)
	target_include_directories(${MODULE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})
endmacro(_SETUP_DEPENDENCIES)  
  
function(buildLibrary LIBRARY_NAME TARGETLIST PREFIX ) 
	_GLOB_SOURCES("${PREFIX}/${LIBRARY_NAME}")
	add_library( ${LIBRARY_NAME} ${CPP_SOURCES} ${CPP_HEADERS} ${CUDA_SOURCES} ${CUDA_HEADERS}) 
	_SETUP_DEPENDENCIES(${LIBRARY_NAME} ${ARGN})
endfunction(buildLibrary) 

_GLOB_SOURCES("cheb")
SET(chebSources "${source_list}")
# MESSAGE(WARNING ${chebSources})
# MESSAGE(WARNING ${source_list})
_GLOB_SOURCES("gui")
SET(guiSources "${source_list}")
# MESSAGE(WARNING ${chebSources})
# MESSAGE(WARNING ${source_list})
_GLOB_SOURCES("imgui")
SET(imguiSources "${source_list}")
_GLOB_SOURCES("simulation")
SET(simSources "${source_list}")
_GLOB_SOURCES("tools")
SET(toolSources "${source_list}")



add_executable(gui Source.cpp ${chebSources} ${guiSources} ${imguiSources} ${simSources} ${toolSources})
_SETUP_DEPENDENCIES(gui glfw glad::glad ) 
# target_compile_options(rayTracingCore PRIVATE /arch:AVX2) 
