cmake_minimum_required(VERSION 3.7.2)
project(ctmm)

enable_testing()

if(MSVC)
    add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
    add_compile_options("-Wall" "-Wextra" "-Werror" "$<$<CONFIG:RELEASE>:-O3>")
endif()

add_subdirectory(src)
add_subdirectory(tests)