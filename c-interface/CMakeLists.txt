add_library(cpgmindexlib STATIC cpgm.cpp)
target_link_libraries(cpgmindexlib pgmindexlib)
target_include_directories(cpgmindexlib PUBLIC ${CMAKE_CURRENT_LIST_DIR})

if (${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif ()