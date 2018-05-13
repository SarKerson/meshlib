cmake_minimum_required(VERSION 2.6)
add_compile_options(-fpermissive -std=c++11 -g)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/")

project(MyDemo)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIR})

find_package(FREEGLUT REQUIRED)
include_directories(${FREEGLUT_INCLUDE_DIR})

# find_package( OpenCV REQUIRED )
# include_directories( ${OpenCV_INCLUDE_DIRS} )

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

include_directories( ${PROJECT_BINARY_DIR}../../MeshLib/algriothm/Conformal/HarmonicMapper )
include_directories( ${PROJECT_BINARY_DIR}../../MeshLib/core )
add_executable(${PROJECT_NAME} ${PROJECT_BINARY_DIR}/../../MeshLib/core/bmp/RgbImage.cpp main.cpp )

target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARY} ${FREEGLUT_LIBRARY} ${OPENGL_LIBRARY} ${OpenCV_LIBS} ${Eigen_LIBS})

