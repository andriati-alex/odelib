cmake_minimum_required(VERSION 3.20)


# Project name and description
project(ODElibrary VERSION 1.0
    DESCRIPTION "Library of routines for ODE system integration"
    LANGUAGES C
)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})


add_library(odelib SHARED src/ode_multistep.c src/ode_singlestep.c src/vector.c)
target_include_directories(odelib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


add_executable(quinney_examples apps/quinney_examples.c)
target_link_libraries(quinney_examples PUBLIC odelib)
set_target_properties(quinney_examples PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)


add_executable(quinney_corrector_iteration apps/quinney_corrector_iteration.c)
target_link_libraries(quinney_corrector_iteration PUBLIC odelib)
set_target_properties(quinney_corrector_iteration PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)


add_executable(adams4order_demo apps/adams4order_demo.c)
target_link_libraries(adams4order_demo PUBLIC odelib)
set_target_properties(adams4order_demo PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)


install(TARGETS odelib DESTINATION lib)
install(TARGETS quinney_examples DESTINATION bin)
install(TARGETS quinney_corrector_iteration DESTINATION bin)
install(TARGETS adams4order_demo DESTINATION bin)
