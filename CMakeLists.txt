cmake_minimum_required(VERSION 3.15)
project(Graduation)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SRC_LIST
    ChiralityMathTools.cpp
    ClassA_Bezier.cpp
    Cone_Surface.cpp
    Cornu_Spiral.cpp
    DiscreteSpiral.cpp
    EulerBezier2D.cpp
    EulerBezier3D.cpp
    EulerBspline2D.cpp
    EulerBspline3D.cpp
    EulerPolygon3D.cpp
    Fillet_using_EB3d.cpp
    FilletSurface.cpp
    FixAxisBezier3D.cpp
    ParameterSurface.cpp
    read3dm.cpp
    Spiral.cpp
    write3dm.cpp
    OBJGenerator.cpp
    LogAestheticBezier.cpp
    EditingCurvatureCurve.cpp
    SphereSkinning.cpp
    BoundaryRecorder.cpp
    )
include_directories(thirdparty/opennurbs)
include_directories(thirdparty/eigen)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Building for windows")
    find_library(OpenNurbs opennurbsStatic.lib ./thirdparty/opennurbs/build/Debug)
    find_library(Z zlib.lib ./thirdparty/opennurbs/build/zlib/Debug)
    MESSAGE(STATUS ${OpenNurbs})
    MESSAGE(STATUS ${Z})
    link_libraries(${OpenNurbs})
    link_libraries(${Z})
    add_library(Chirality_LIB STATIC ${SRC_LIST})
    add_executable(Test main.cpp)
    target_link_libraries(Test PUBLIC Chirality_LIB)
    add_executable(Cup Cup_main.cpp)
    target_link_libraries(Cup PUBLIC Chirality_LIB)
    add_executable(LogAesthetic LogAesthetic_main.cpp)
    target_link_libraries(LogAesthetic PUBLIC Chirality_LIB)
    add_executable(ElectricFan ElectricFan_main.cpp)
    target_link_libraries(ElectricFan PUBLIC Chirality_LIB)
    add_executable(SphereSkin sphere_skin_main.cpp)
    target_link_libraries(SphereSkin PUBLIC Chirality_LIB)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    MESSAGE(STATUS "Building for linux")
    find_library(OpenNurbs libOpenNURBS.so ./thirdparty/opennurbs/build)
    #find_library(Z zlib.lib ./thirdparty/opennurbs/build/zlib/Debug)
    MESSAGE(STATUS ${OpenNurbs})
    #MESSAGE(STATUS ${Z})
    link_libraries(${OpenNurbs})
    link_libraries(pthread)
    #link_libraries(${Z})
    add_executable(${CMAKE_PROJECT_NAME} ${SRC_LIST})
endif()