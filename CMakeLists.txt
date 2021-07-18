cmake_minimum_required(VERSION 3.20)


# Project name and description
project(ODElibrary VERSION 1.0
    DESCRIPTION "Library of routines for ODE system integration"
    LANGUAGES C
)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})


add_library(odelib SHARED src/ode_multistep.c src/ode_singlestep.c)
target_include_directories(odelib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_link_libraries(linalglib m)


add_executable(test apps/test.c)
target_link_libraries(test PUBLIC odelib)
set_target_properties(test PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)


install(TARGETS odelib DESTINATION lib)
install(TARGETS test DESTINATION bin)


# To install in non-standard directory a shared lib see the Zaffy answer
# https://stackoverflow.com/questions/58997230/cmake-project-fails-to-find-shared-library
# from where I get the line `set_target_properties`, Another place to look at is
# https://github.com/pablospe/cmake-example-library/blob/master/cmake/SetEnv.cmake
# where the following variable is set
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Any problem search for more info about runtime path (cmake rpath)
