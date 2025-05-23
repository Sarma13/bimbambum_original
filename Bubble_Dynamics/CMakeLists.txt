# -- Set the project --
#

CMAKE_MINIMUM_REQUIRED (VERSION 3.18 FATAL_ERROR)
IF (POLICY CMP0042)
  # policies not known to CMake 2.8 ...
  CMAKE_POLICY (SET CMP0042 NEW) # Mac OS -rpath behavior
ENDIF (POLICY CMP0042)

# Compiler environment
SET (CMAKE_CXX_STANDARD 14)
SET (CMAKE_CXX_FLAGS -Wall)

# Project information
SET (${PROJECT_NAME}_VERSION_MAJOR 0)
SET (${PROJECT_NAME}_VERSION_MINOR 0)
SET (${PROJECT_NAME}_VERSION_PATCH 1)
PROJECT (BimBamBum)
SET (PROJECT_VERSION 0.0.1)

# Set output path
SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -- Include headers and sources for the project --
#

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/include/stream
  ${PROJECT_SOURCE_DIR}/include/core)
FILE (GLOB SOURCES_STREAM ${PROJECT_SOURCE_DIR}/src/stream/*.cpp)
FILE (GLOB SOURCES_CORE ${PROJECT_SOURCE_DIR}/src/core/*.cpp)
FILE (GLOB HEADERS_STREAM ${PROJECT_SOURCE_DIR}/include/stream/*.hpp)
FILE (GLOB HEADERS_CORE ${PROJECT_SOURCE_DIR}/include/core/*.hpp)

# -- Include the examples from the examples folder --
#

FILE (GLOB SOURCES_EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/solver/*.cpp)
SET (EXAMPLES_CONFIG_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_free_surface_gamma056.json
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_free_surface_gamma075.json
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_Plesset_unbounded_epsilon100.json
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_Plesset_unbounded_epsilon1000.json
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_rigid_boundary.json
    ${CMAKE_CURRENT_SOURCE_DIR}/config_examples/inputs_Rayleigh_unbounded.json)

# -- Find libraries --
#

# GSL
FIND_PACKAGE(GSL REQUIRED)
IF (GSL_FOUND)
  MESSAGE ("-- GSL found")
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
ELSEIF (NOT GSL_FOUND)
  MESSAGE (FATAL_ERROR "-- GSL not found")
ENDIF()

# Boost
FIND_PACKAGE(Boost REQUIRED)
IF (Boost_FOUND)
  MESSAGE ("-- Boost found")
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSEIF (NOT Boost_FOUND)
  MESSAGE (FATAL_ERROR "-- Boost not found")
ENDIF()

# Armadillo
FIND_PACKAGE(Armadillo REQUIRED)
IF (Armadillo_FOUND)
  MESSAGE ("-- Armadillo found")
  INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
ELSEIF (NOT Armadillo_FOUND)
  MESSAGE (FATAL_ERROR "-- Armadillo not found")
ENDIF()

# Open MP

IF (BUILD_OPENMP)


    FIND_PACKAGE(OpenMP)
    IF(OpenMP_FOUND)
      MESSAGE ("-- OpenMP found")
      SET (OPENMP_LIBRARIES OpenMP::OpenMP_CXX)
    ELSEIF (NOT OpenMP_FOUND)
      MESSAGE (FATAL_ERROR "-- OpenMP not found")
    ENDIF()

  
ENDIF()



# -- Doxygen documentation --
#

IF (BUILD_DOXYGEN)
  FIND_PACKAGE(Doxygen REQUIRED)

  if( DOXYGEN_FOUND)
    SET(DOXYGEN_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)
    ADD_CUSTOM_TARGET(doxydoc
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
      COMMENT "Generating Doxygen documentation")
    ADD_DEPENDENCIES (doxydoc ${PROJECT_NAME})
  ENDIF()

ENDIF()

# -- Link the library and install --
#

ADD_LIBRARY(${PROJECT_NAME} SHARED ${SOURCES_STREAM} ${SOURCES_CORE})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME}
  SYSTEM
  PUBLIC ${GSL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${OPENMP_INCLUDE_DIRS}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/stream>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/core>
  PUBLIC $<INSTALL_INTERFACE:include/BimBamBum>
  PUBLIC $<INSTALL_INTERFACE:include/BimBamBum/stream>
  PUBLIC $<INSTALL_INTERFACE:include/BimBamBum/core>
  )
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
  PUBLIC ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} ${OPENMP_LIBRARIES})
SET_TARGET_PROPERTIES(${PROJECT_NAME}
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  LIBRARY_OUTPUT_DIRECTORY lib/
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib
  )

INSTALL (TARGETS
  ${PROJECT_NAME}
  EXPORT BimBamBumConfig
  ARCHIVE DESTINATION lib/
  LIBRARY DESTINATION lib/
  RUNTIME DESTINATION bin/
  )

INSTALL( DIRECTORY include/ # source directory
         DESTINATION include/BimBamBum # target directory
         FILES_MATCHING PATTERN "*.h" PATTERN "*.hh" PATTERN "*.hpp" # select header files
       )

EXPORT (TARGETS
  ${PROJECT_NAME}
  NAMESPACE BimBamBum::
  FILE share/cmake/BimBamBumConfig.cmake
  )

INSTALL (EXPORT
  BimBamBumConfig
  DESTINATION lib/cmake/${PROJECT_NAME}-${PROJECT_VERSION} #< installed by pkg's prefix, see https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html#package-configuration-file
  NAMESPACE BimBamBum::
  )

# -- Examples --
#

# Add the build targets for each example executable
FOREACH( SOURCE_EXAMPLE ${SOURCES_EXAMPLES} )
  GET_FILENAME_COMPONENT ( EXAMPLE_NAME ${SOURCE_EXAMPLE} NAME_WE )
  ADD_EXECUTABLE ( ${EXAMPLE_NAME} ${SOURCE_EXAMPLE} )
  TARGET_LINK_LIBRARIES ( ${EXAMPLE_NAME} ${PROJECT_NAME} )
  SET_TARGET_PROPERTIES ( ${EXAMPLE_NAME}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    EXCLUDE_FROM_ALL TRUE)

  # Copy the necessary configuration files
  FOREACH (CONFIG_FILE ${EXAMPLES_CONFIG_FILES})
    GET_FILENAME_COMPONENT (CONFIG_NAME ${CONFIG_FILE} NAME)
    CONFIGURE_FILE (${CONFIG_FILE}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CONFIG_NAME}
      COPYONLY)
  ENDFOREACH()

  LIST ( APPEND EXAMPLE_NAMES ${EXAMPLE_NAME} )
ENDFOREACH()

# Add a custom target to build all examples
ADD_CUSTOM_TARGET ( solver )
FOREACH (EXAMPLE_NAME ${EXAMPLE_NAMES} )
  ADD_DEPENDENCIES ( solver ${EXAMPLE_NAME} )
ENDFOREACH()
