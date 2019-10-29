add_executable(tests ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../lib/catch.hpp)
target_link_libraries(tests pgmindexlib)
