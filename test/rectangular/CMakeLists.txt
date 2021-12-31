dune_symlink_to_source_files(FILES "params.input" "channel.msh")
dumux_add_test(NAME test_channel SOURCES main.cc)
