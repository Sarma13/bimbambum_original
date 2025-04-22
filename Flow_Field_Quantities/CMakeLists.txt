cmake_minimum_required(VERSION 3.18)

project(Flow_Field_Quantities)

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

SET (PYTHON_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/solver_codes/python/inputs.py
		${CMAKE_CURRENT_SOURCE_DIR}/solver_codes/python/main.py
		${CMAKE_CURRENT_SOURCE_DIR}/solver_codes/python/flow_field.py
		${CMAKE_CURRENT_SOURCE_DIR}/solver_codes/python/create_domain.py
		${CMAKE_CURRENT_SOURCE_DIR}/solver_codes/python/classify_domain.py)

SET (EXAMPLE_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/simulation_examples/rigid_boundary_Rayleigh_gamma095.txt
		${CMAKE_CURRENT_SOURCE_DIR}/simulation_examples/unbounded_Rayleigh.txt)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)

pybind11_add_module(flow_potential
			solver_codes/cpp/velocity_potential.cpp
			solver_codes/cpp/cubic_spline.cpp)

# -- Check for Python 3 --
#

FIND_PACKAGE(Python 3 COMPONENTS Interpreter Development REQUIRED)
IF(Python_FOUND)
 MESSAGE("-- Python 3 found")
ELSEIF(NOT Python_FOUND)
 MESSAGE(FATAL_ERROR "-- Python 3 not found")
ENDIF()

# -- Find dependencies --
#

# Armadillo
FIND_PACKAGE(Armadillo REQUIRED)
IF (Armadillo_FOUND)
	MESSAGE ("-- Armadillo found")
	INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
	TARGET_LINK_LIBRARIES(flow_potential PRIVATE ${ARMADILLO_LIBRARIES})
ELSEIF (NOT Armadillo_FOUND)
	MESSAGE (FATAL_ERROR "-- Armadillo not found")
ENDIF()

# GSL
FIND_PACKAGE(GSL REQUIRED)
IF (GSL_FOUND)
	MESSAGE ("-- GSL found")
	INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
	TARGET_LINK_LIBRARIES(flow_potential PRIVATE ${GSL_LIBRARIES})
ELSEIF (NOT GSL_FOUND)
	MESSAGE (FATAL_ERROR "-- GSL not found")
ENDIF()



# -- Copy the necessary python files --
# 

FOREACH (PYTHON_SCRIPT ${PYTHON_FILES})
	GET_FILENAME_COMPONENT (PYTHON_NAME ${PYTHON_SCRIPT} NAME)
	CONFIGURE_FILE (${PYTHON_SCRIPT}
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PYTHON_NAME}
			COPYONLY)
ENDFOREACH()

# -- Copy the necessary example files --
# 

FOREACH (EXAMPLE_INPUTS ${EXAMPLE_FILES})
	GET_FILENAME_COMPONENT (EXAMPLE_NAME ${EXAMPLE_INPUTS} NAME)
	CONFIGURE_FILE (${EXAMPLE_INPUTS}
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXAMPLE_NAME}
			COPYONLY)
ENDFOREACH()
