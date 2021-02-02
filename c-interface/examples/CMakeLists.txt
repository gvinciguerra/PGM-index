set(examples
        simple
        updates)

foreach (example ${examples})
    add_executable(example_${example}_c ${example}.c)
    target_link_libraries(example_${example}_c cpgmindexlib)
endforeach (example)