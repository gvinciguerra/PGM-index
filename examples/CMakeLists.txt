set(examples
        mapped
        multidimensional
        multiset
        simple
        updates)

foreach (example ${examples})
    add_executable(example_${example} ${example}.cpp)
    target_link_libraries(example_${example} pgmindexlib)
endforeach (example)