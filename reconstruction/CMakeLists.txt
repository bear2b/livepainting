# Build the realtime deformable surface 3D tracking and reconstruction
#
# Author: Ngo Tien Dat, datngotien@gmail.com
# 25 Sep 2014

cmake_minimum_required (VERSION 2.6)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_COMPILER /usr/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

set(LIB_NAME deformable)

# Find opencv
find_package( OpenCV REQUIRED )

# Find armadillo
find_library( armadillo NAMES "armadillo" PATHS /opt/local/lib/ )
message(STATUS "Found armadillo in ${armadillo}")

# Find Eigen
find_library( eigen NAMES "Eigen" PATHS /opt/local/lib/ )
message(STATUS "Found armadillo in ${Eigen}")

# Find ferns
#find_library( ferns NAMES "ferns" PATHS "${PROJECT_SOURCE_DIR}/3rdparty/ferns/" )
#message(STATUS "Found ferns in ${ferns}")

# Find timer
#find_library( timer NAMES "timer" PATHS "${PROJECT_SOURCE_DIR}/3rdparty/timer/" )
#message(STATUS "Found timer in ${timer}")

# Inlude directories
include_directories( ${OpenCV_INCLUDE_DIRS} )           # for OpenCV
include_directories( "${PROJECT_SOURCE_DIR}/src" )
#include_directories( "${PROJECT_SOURCE_DIR}/3rdparty" ) # for ferns and timer
include_directories( "/opt/local/include/" )            # for armadillo

set(SOURCE_FILES 

  ./data/im_corners.txt
  ./src/maillage.cpp
   ./src/maillage.hpp
  ./src/RealtimeDemo.cpp
  ./src/RealtimeDemo.hpp
./src/LaplacianMesh.cpp
   ./src/LaplacianMesh.hpp
./src/LinearAlgebraUtils.cpp
./src/LinearAlgebraUtils.hpp
./src/Camera.cpp
./src/Camera.hpp
  
#  ./src/Camera.cpp
  # ./src/Linear/EqualConstrFunction.cpp
  #  ./src/Linear/EqualConstrOptimize.cpp
  #  ./src/Linear/IneqConstrFunction.cpp
  #  ./src/Linear/IneqConstrOptimize.cpp
  #  ./src/Linear/LinearAlgebraUtils.cpp
  #  ./src/Linear/ObjectiveFunction.cpp
  #  ./src/MatchTracker/InlierMatchTracker.cpp
  #  ./src/MatchTracker/LKPointTracker.cpp
  #  ./src/PointMatching/FernKeypointMatcher3D2D.cpp
   ./src/KeypointMatcher3D2D.cpp
   ./src/KeypointMatcher3D2D.hpp
  #  ./src/RealtimeDemo/RealtimeDemo.cpp
  #  ./src/Reconstruction.cpp
    ./src/DUtils.cpp
    ./src/DUtils.hpp
  ./src/Visualization.hpp
  ./src/Visualization.cpp
  ./src/DefinedMacros.hpp
    ./src/Timer.cpp
    ./src/Timer.hpp
)

set( LIBS ${armadillo} ${OpenCV_LIBS} ${Eigen}) #${ferns} ${timer})

# Build library
add_library( ${LIB_NAME} SHARED ${SOURCE_FILES} )
target_link_libraries( ${LIB_NAME} ${LIBS})

# Build deformabletracker
set(APP_NAME firstversion)
add_executable( ${APP_NAME} ./main.cpp )
target_link_libraries( ${APP_NAME} ${LIB_NAME} ${LIBS} )



