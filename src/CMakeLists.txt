cmake_minimum_required(VERSION 3.7.2)
project(ctmm)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(${PROJECT_NAME} SHARED ctmm.c)

if (NOT WIN32)
    target_link_libraries(${PROJECT_NAME} m)
endif ()

add_custom_command(TARGET ctmm
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/../python/pyctmm/
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../lib/$<TARGET_FILE_NAME:ctmm>" "${CMAKE_CURRENT_SOURCE_DIR}/../python/pyctmm/"
    COMMENT "Copy shared library file to ${CMAKE_CURRENT_SOURCE_DIR}/../python/pyctmm/ directory" VERBATIM
)

add_custom_command(TARGET ctmm
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/../python/include/
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../include/ctmm.h" "${CMAKE_CURRENT_SOURCE_DIR}/../python/include/"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../include/DLL_defines.h" "${CMAKE_CURRENT_SOURCE_DIR}/../python/include/"
    COMMENT "Copy ctmm header files to ${CMAKE_CURRENT_SOURCE_DIR}/../python/include/ directory" VERBATIM
)

if (WIN32)
    add_custom_command(TARGET ctmm
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../lib/ctmm.lib" "${CMAKE_CURRENT_SOURCE_DIR}/../python/pyctmm/"
		COMMENT "Copy import library file to ${CMAKE_CURRENT_SOURCE_DIR}/../python/pyctmm/ directory" VERBATIM
	)
endif ()