add_subdirectory(en/source)
add_subdirectory(ja/source)
add_subdirectory(tutorial/en/source)

add_custom_target(doc DEPENDS doc-ja doc-en)
add_custom_target(tutorial DEPENDS tutorial-en)
