cmake_minimum_required(VERSION 3.5)

# Project Name
PROJECT(phase1)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_executable(${PROJECT_NAME}
    engine/main.cpp
    engine/transformations.cpp
    engine/tinyxml2.cpp
)

include_directories(
    engine
    TinyXML
)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})

if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)


find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})

if(WIN32)
    message(STATUS "Toolkits_DIR set to: " ${TOOLKITS_FOLDER})
    set(TOOLKITS_FOLDER "" CACHE PATH "Path to Toolkits folder")
    
    if(NOT EXISTS "${TOOLKITS_FOLDER}/glut/GL/glut.h" OR NOT EXISTS "${TOOLKITS_FOLDER}/glut/glut32.lib")
        message(ERROR ": GLUT not found")
    endif()
    
    include_directories(${TOOLKITS_FOLDER}/glut)
    target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} 
                                      ${TOOLKITS_FOLDER}/glut/glut32.lib
                                      glew32)  # Add GLEW for Windows
    
    if(EXISTS "${TOOLKITS_FOLDER}/glut/glut32.dll")
        file(COPY ${TOOLKITS_FOLDER}/glut/glut32.dll DESTINATION ${CMAKE_BINARY_DIR})
    endif()
    
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
    
else() #Linux and Mac
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
    
    find_package(GLUT REQUIRED)
    find_package(GLEW REQUIRED)  # Add this line to find GLEW
    
    include_directories(${GLUT_INCLUDE_DIR} ${GLEW_INCLUDE_DIRS})
    link_directories(${GLUT_LIBRARY_DIRS} ${GLEW_LIBRARY_DIRS})
    add_definitions(${GLUT_DEFINITIONS})

    find_package(DevIL REQUIRED)
	link_libraries(${IL_LIBRARIES} )
	include_directories(ยง{IL_INCLUDE_DIR})
    
    target_link_libraries(${PROJECT_NAME} 
        ${OPENGL_LIBRARIES} 
        ${GLUT_LIBRARIES}
        ${GLEW_LIBRARIES}  # Link GLEW
        IL
    )
    
    if(NOT GLUT_FOUND)
        message(ERROR ": GLUT not found!")
    endif()
    
    if(NOT GLEW_FOUND)
        message(ERROR ": GLEW not found!")
    endif()
    
endif()