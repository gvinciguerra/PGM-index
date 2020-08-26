add_executable(tests tests.cpp)
target_link_libraries(tests pgmindexlib)
add_test(NAME test_all COMMAND tests)