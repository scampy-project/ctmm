cmake_minimum_required(VERSION 3.7.2)
project(ctmm)

SET(CTEST_OUTPUT_ON_FAILURE TRUE) 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

add_executable(ctmm_test ctmm_test.c)
add_executable(ctmm_mathematics_test ctmm_mathematics_test.c)
add_executable(ctmm_fresnel_test ctmm_fresnel_test.c)


target_link_libraries(ctmm_test ctmm)
target_link_libraries(ctmm_mathematics_test ctmm)
target_link_libraries(ctmm_fresnel_test ctmm)

add_test(NAME ctmm_mathematics_test COMMAND ctmm_mathematics_test)
add_test(NAME ctmm_test COMMAND ctmm_test)

add_custom_command(TARGET ctmm_test ctmm_mathematics_test
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../lib/$<TARGET_FILE_NAME:ctmm>" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    COMMENT "Copy shared library file to ${CMAKE_CURRENT_BINARY_DIR} directory" VERBATIM
)