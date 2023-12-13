cmake_minimum_required(VERSION "3.26.4")

project(PGBAR)

add_executable(${PROJECT_NAME} src/analyze_data.cxx)

add_library(constants src/constants.cpp)
add_library(reparam src/reparam.cpp)
add_library(param src/param.cpp)
add_library(utils src/utils.cpp)
add_library(particle src/particle.cpp)
target_link_libraries(particle PUBLIC utils param reparam constants)

target_link_libraries(${PROJECT_NAME} PUBLIC utils constants param reparam particle)
target_link_libraries(${PROJECT_NAME} PUBLIC gsl jsoncpp openblas armadillo)
target_include_directories(${PROJECT_NAME} PUBLIC ../include/)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR}/bin)





