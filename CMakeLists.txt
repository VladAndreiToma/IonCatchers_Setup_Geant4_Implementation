#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(IonCatchers)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
find_package(ROOT REQUIRED)

# Find CADMesh
find_package(cadmesh REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
		    ${cadmesh_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
add_executable(IonCatchers IonCatchers.cc ${sources} ${headers})
target_link_libraries(IonCatchers ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}  ${cadmesh_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build director. This is so that we can run the executable directly
# because it relies on these scripts being in the current working directory.
set(IonCatchers_SCRIPTS
  vis.mac
  JobsRun.pl
  JobsProcess.pl
  )

foreach(_script ${IonCatchers_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
