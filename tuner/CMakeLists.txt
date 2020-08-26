add_subdirectory(alglib)
add_executable(tuner tuner.cpp)
target_link_libraries(tuner alglib pgmindexlib)