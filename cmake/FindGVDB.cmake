#####################################################################################
# Find GVDB
#
unset(GVDB_FOUND CACHE)
unset(GVDB_INCLUDE_DIR CACHE)

if ( NOT DEFINED GVDB_ROOT_DIR )
  if (WIN32)
    get_filename_component ( BASEDIR "${CMAKE_MODULE_PATH}/../../_output" REALPATH )
  else()
    get_filename_component ( BASEDIR "/usr/local/gvdb/" REALPATH )
  endif()
  set ( GVDB_ROOT_DIR ${BASEDIR} CACHE PATH "Location of GVDB library" FORCE)
endif()
message ( STATUS "Searching for GVDB at.. ${GVDB_ROOT_DIR}")
set( GVDB_FOUND "YES" )

#----------------------------------------------- CROSS-PLATFORM FIND FILES
# Find one or more of a specific file in the given folder
# Returns the file name w/o path

macro(_FIND_FILE targetVar searchDir nameWin64 nameLnx cnt)
  unset ( fileList )  
  unset ( nameFind )
  unset ( targetVar )  
  if ( WIN32 ) 
     SET ( nameFind ${nameWin64} )
  else()
     SET ( nameFind ${nameLnx} )
  endif()
  if ( "${nameFind}" STREQUAL ""  )
    MATH (EXPR ${cnt} "${${cnt}}+1" )	
  else()
    file(GLOB fileList "${${searchDir}}/${nameFind}")  
    list(LENGTH fileList NUMLIST)  
    if (NUMLIST GREATER 0)	
       MATH (EXPR ${cnt} "${${cnt}}+1" )	
       list(APPEND ${targetVar} ${nameFind} )
    endif() 
  endif()
endmacro()

macro(debug msg)
    message(STATUS "DEBUG ${msg}")
endmacro()
macro(debugValue variableName)
    debug("${variableName}=\${${variableName}}")
endmacro()

#----------------------------------------------- CROSS-PLATFORM FIND MULTIPLE
# Find all files in specified folder with the given extension.
# This creates a file list, where each entry is only the filename w/o path
# Return the count of files
macro(_FIND_MULTIPLE targetVar searchDir extWin64 extLnx cnt)    
  unset ( fileList )    
  unset ( targetVar ) 
  unset ( ${cnt} )
  set ( ${cnt} "0" )
  if ( WIN32 ) 
     SET ( extFind ${extWin64} )
  else()
     SET ( extFind ${extLnx} )
  endif()
  file( GLOB fileList "${${searchDir}}/*.${extFind}")  
  list( LENGTH fileList NUMLIST)
  math( EXPR ${cnt} "${${cnt}}+${NUMLIST}" )  
  foreach ( _file ${fileList} )
     get_filename_component ( fname ${_file} NAME )
     list( APPEND ${targetVar} ${fname} )
  endforeach()
endmacro()
if ( GVDB_ROOT_DIR )

    #-- Paths to GVDB Library (cached so user can modify)
	set ( GVDB_INCLUDE_DIR "${GVDB_ROOT_DIR}/include" CACHE PATH "Path to include files" FORCE)
	set ( GVDB_LIB_DIR "${GVDB_ROOT_DIR}/lib" CACHE PATH "Path to libraries" FORCE)	
	set ( GVDB_SHARE_DIR "${GVDB_ROOT_DIR}/lib" CACHE PATH "Path to share files" FORCE)	
    set ( GVDB_BIN_DIR "${GVDB_ROOT_DIR}/bin" CACHE PATH "Path to binary files" FORCE)	

	#-------- Locate Header files
        set ( OK_H "0" )
	_FIND_FILE ( GVDB_HEADERS GVDB_INCLUDE_DIR "gvdb.h" "gvdb.h" OK_H )
	if ( OK_H EQUAL 1 ) 
	    message ( STATUS "  Found. GVDB Header files. ${GVDB_INCLUDE_DIR}" )
	else()
	    message ( "  NOT FOUND. GVDB Header files" )
		set ( GVDB_FOUND "NO" )
	endif ()

    #-------- Locate Library	
        set ( OK_DLL 0 )	
        set ( OK_LIB 0 )	
	_FIND_FILE ( LIST_DLL GVDB_LIB_DIR "gvdb.dll" "libgvdb.so" OK_DLL )	
  	_FIND_FILE ( LIST_LIB GVDB_LIB_DIR "gvdb.lib" "libgvdb.so" OK_LIB )
	if ( (${OK_LIB} EQUAL 1) ) 
	   message ( STATUS "  Found. GVDB Library. ${GVDB_LIB_DIR}" )	   
	else()
	   set ( GVDB_FOUND "NO" )	   
	   message ( "  NOT FOUND. GVDB Library. (so/dll or lib missing)" )	   
	endif()

	#-------- Locate PTX/GLSL
        set ( OK_PTX 0 )	
        set ( OK_GLSL 0 )	
	_FIND_MULTIPLE( LIST_PTX GVDB_BIN_DIR "ptx" "ptx" OK_PTX )  
         _FIND_MULTIPLE( LIST_GLSL GVDB_BIN_DIR "glsl" "glsl" OK_GLSL )    
	if ( (${OK_PTX} EQUAL 3) AND (${OK_GLSL} GREATER 2) ) 
	   message ( STATUS "  Found. GVDB Ptx/Glsl. ${GVDB_SHARE_DIR}" )	   
	else()
	   set ( GVDB_FOUND "NO" )	   
	   message ( "  NOT FOUND. GVDB Ptx/Glsl. (ptx or glsl missing)" )	   	   
	endif()	
	
	#-------- Locate Support DLLs
        set ( OK_EXTRA "0" )		
	_FIND_MULTIPLE ( LIST_EXTRA GVDB_LIB_DIR "dll" "so" OK_EXTRA )	

endif()
 
if ( ${GVDB_FOUND} STREQUAL "NO" )
   message( FATAL_ERROR "
      Please set GVDB_ROOT_DIR to the root location 
      of installed GVDB library containing /include and /lib.
      Not found at GVDB_ROOT_DIR: ${GVDB_ROOT_DIR}\n"
   )
endif()

set ( GVDB_DLL ${LIST_DLL} CACHE INTERNAL "" FORCE)
set ( GVDB_LIB ${LIST_LIB} CACHE INTERNAL "" FORCE)
set ( GVDB_PTX ${LIST_PTX} CACHE INTERNAL "" FORCE)
set ( GVDB_GLSL ${LIST_GLSL} CACHE INTERNAL "" FORCE)
set ( GVDB_EXTRA ${LIST_EXTRA} CACHE INTERNAL "" FORCE)

#-- Create a list of all binary files needed for exes
unset ( LIST_FULL )	
LIST ( APPEND LIST_FULL ${LIST_DLL} )
LIST ( APPEND LIST_FULL ${LIST_PTX} )	
LIST ( APPEND LIST_FULL ${LIST_GLSL} )	
LIST ( APPEND LIST_FULL ${LIST_EXTRA} )	
set ( GVDB_LIST ${LIST_FULL} CACHE INTERNAL "" )

#-- We do not want user to modified these vars, but helpful to show them
message ( STATUS "  GVDB_ROOT_DIR: ${GVDB_ROOT_DIR}" )
message ( STATUS "  GVDB_DLL:  ${GVDB_DLL}" )
message ( STATUS "  GVDB_LIB:  ${GVDB_LIB}" )
message ( STATUS "  GVDB_PTX:  ${GVDB_PTX}" )
message ( STATUS "  GVDB_GLSL: ${GVDB_GLSL}" )
message ( STATUS "  GVDB_EXTRA:${GVDB_EXTRA}" )

mark_as_advanced(GVDB_FOUND)