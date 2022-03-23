set(target
        build
        benchmark)

foreach (target ${target})
    add_executable(kv_${target} ${target}.cpp)
    target_link_libraries(kv_${target} pgmindexlib)
endforeach (target)
