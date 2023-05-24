project(driven_QME)

set(OPENBLAS_INCLUDE_DIRS /opt/homebrew/Cellar/openblas/0.3.21/include/)
set(OPENBLAS_LIBRARIES /opt/homebrew/Cellar/openblas/0.3.21/lib/libopenblasp-r0.3.21.dylib)
set(ARMADILLO_INCLUDE_DIRS /opt/homebrew/Cellar/armadillo/11.4.2/include/)
set(ARMADILLO_LIBRARIES /opt/homebrew/Cellar/armadillo/11.4.2/lib/libarmadillo.11.4.2.dylib)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${OPENBLAS_INCLUDE_DIRS})

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(driven_QME main.cpp constant.h constant.cpp ioput.h ioput.cpp initialize.h initialize.cpp evolution.h evolution.cpp)
target_link_libraries(driven_QME PRIVATE ${ARMADILLO_LIBRARIES})
enable_testing()

add_test(NAME driven_QME COMMAND driven_QME)
