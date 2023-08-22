cmake_minimum_required(VERSION 3.5)

project(macro_dr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#find_package(Catch2 REQUIRED)
#target_link_libraries(tests Catch2::Catch2)
add_executable(macro_dr main.cpp
   # CLI_base.h
   Maybe_error.h
   CLI_regular_types.h
   qmodel.h
   run_script.txt
   CLI_macro_dr.h
   CLI_grammar.h
   grammar_untyped.h
    grammar_Identifier.h
    lexer_untyped.h
    test.txt
#    tst_untyped_grammar.cpp
)

include(GNUInstallDirs)
install(TARGETS macro_dr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
