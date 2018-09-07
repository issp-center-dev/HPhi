add_subdirectory(en/source)
add_subdirectory(ja/source)

add_custom_target(doc DEPENDS doc-ja doc-en)
