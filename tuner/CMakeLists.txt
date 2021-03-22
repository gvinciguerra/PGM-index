include_directories(../benchmark)
add_executable(tuner tuner.cpp)
target_link_libraries(tuner pgmindexlib)