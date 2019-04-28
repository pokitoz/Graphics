cmake_minimum_required (VERSION 2.6)
project(rawtracing)

add_executable(${PROJECT_NAME} skeleton.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package (SDL)

if(NOT SDL_FOUND)
   message(FATAL_ERROR "SDL not found!")
endif()

include_directories(
	${SDL_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/glm
)

target_link_libraries(${PROJECT_NAME} ${SDL_LIBRARY})

