cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(sd1d LANGUAGES CXX)

# Extra CMake scripts in cmake/ subdirectory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Update submodules
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message(STATUS "Submodule update")
  execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT)
  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
  endif()
endif()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC SD1D_REVISION)
if(SD1D_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(SD1D_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${SD1D_REVISION}")

# BOUT++ is a dependency
set(BOUT_BUILD_EXAMPLES, OFF) # Don't create example makefiles
add_subdirectory(external/BOUT-dev)

set(SD1D_SOURCES
    sd1d.cxx
    div_ops.cxx
    loadmetric.cxx
    radiation.cxx
    atomicpp/ImpuritySpecies.cxx
    atomicpp/Prad.cxx
    atomicpp/RateCoefficient.cxx
    atomicpp/sharedFunctions.cxx
    div_ops.hxx
    loadmetric.hxx
    radiation.hxx
    atomicpp/ImpuritySpecies.hxx
    atomicpp/json.hxx
    atomicpp/Prad.hxx
    atomicpp/RateCoefficient.hxx
    atomicpp/sharedFunctions.hxx
    )

# The main executable target
add_executable(sd1d
               ${SD1D_SOURCES}
               ${CMAKE_CURRENT_BINARY_DIR}/revision.hxx)

target_link_libraries(sd1d PRIVATE bout++::bout++)

target_include_directories(sd1d PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

# Build the file containing just the commit hash
# This will be rebuilt on every commit!
configure_file(
  "${PROJECT_SOURCE_DIR}/revision.hxx.in"
  "${PROJECT_BINARY_DIR}/revision.hxx")

# Once built, copy the data and test directories
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/case-01 $<TARGET_FILE_DIR:${PROJECT_NAME}>/case-01)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/case-02 $<TARGET_FILE_DIR:${PROJECT_NAME}>/case-02)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/case-03 $<TARGET_FILE_DIR:${PROJECT_NAME}>/case-03)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/case-04 $<TARGET_FILE_DIR:${PROJECT_NAME}>/case-04)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/case-05 $<TARGET_FILE_DIR:${PROJECT_NAME}>/case-05)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/MAST-U $<TARGET_FILE_DIR:${PROJECT_NAME}>/MAST-U)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/json_database $<TARGET_FILE_DIR:${PROJECT_NAME}>/json_database)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                   ${CMAKE_SOURCE_DIR}/impurity_user_input.json $<TARGET_FILE_DIR:${PROJECT_NAME}>/impurity_user_input.json)
