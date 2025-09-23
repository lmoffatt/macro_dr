# ──────────────────────────────────────────────────────────────────────
#  Macro DR – unified top-level CMakeLists.txt
# ──────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.18)
project(macro_dr LANGUAGES CXX)

# -------------------------------------------------------
# 1) Global compile rules              (kept from legacy)
# -------------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optional: integrate clang-tidy if available
option(MACRODR_ENABLE_CLANG_TIDY "Run clang-tidy during build" OFF)
if (MACRODR_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy)
    if (CLANG_TIDY_EXE)
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
        message(STATUS "clang-tidy enabled: ${CLANG_TIDY_EXE}")
    else()
        message(WARNING "clang-tidy requested but not found")
    endif()
endif()

# Common warnings
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Werror=return-type
)

# Per-compiler flags (Clang vs GCC)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(
        -ferror-limit=5
        # -ftime-trace  # enable locally if needed; produces JSON per object file
    )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        -fmax-errors=5
        -ftime-report
        -Wno-ctad-maybe-unsupported
    )
endif()

set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")


# -------------------------------------------------------
# 2) Obtain Git commit hash (unchanged)
# -------------------------------------------------------
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(GIT_COMMIT_HASH "unknown")
endif()

# -------------------------------------------------------
# 3) Third-party single-header interface target ### NEW
# -------------------------------------------------------
# 3) Third-party single-header interface target  ───────────────
add_library(macrodr_thirdparty INTERFACE)

# > BUILD_INTERFACE  – during build  → point at source tree
# > INSTALL_INTERFACE – after install → point inside include dir
target_include_directories(macrodr_thirdparty SYSTEM INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/macrodr/third_party>)

# Copy the header(s) so the install-side path is valid
install(DIRECTORY ${PROJECT_SOURCE_DIR}/third_party/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/macrodr/third_party
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# -------------------------------------------------------
# 4) Find heavy dependencies (OpenMP, BLAS/LAPACK, GSL)
# -------------------------------------------------------
find_package(OpenMP REQUIRED)                                    # 
find_package(BLAS)
find_package(LAPACK)
set(lapackblas_libraries "")
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    message(WARNING "BLAS/LAPACK not found – install libblas-dev/liblapack-dev or set MACRODR_BLAS_DIR")
    if(DEFINED ENV{MACRODR_BLAS_DIR})
        set(_macrodr_blas_root $ENV{MACRODR_BLAS_DIR})
        set(_macrodr_blas_lib "${_macrodr_blas_root}/lib/libcblas.so")
        set(_macrodr_lapack_lib "${_macrodr_blas_root}/lib/liblapack.so")
        if(EXISTS ${_macrodr_blas_lib} AND EXISTS ${_macrodr_lapack_lib})
            set(lapackblas_libraries ${_macrodr_blas_lib} ${_macrodr_lapack_lib})
        else()
            message(FATAL_ERROR "BLAS/LAPACK not found at ${_macrodr_blas_root}")
        endif()
    elseif(EXISTS "/nfs/home/lmoffatt/.nix-profile/lib/libcblas.so" AND EXISTS "/nfs/home/lmoffatt/.nix-profile/lib/liblapack.so")
        message(STATUS "Using fallback Nix BLAS/LAPACK")
        set(lapackblas_libraries
            /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so
            /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so)
    else()
        message(FATAL_ERROR "BLAS/LAPACK not found")
    endif()
endif()

find_package(GSL)
set(gsl_library "")
set(gsl_dir "")
if(GSL_FOUND)
    set(gsl_library ${GSL_LIBRARY})
    set(gsl_dir     ${GSL_INCLUDE_DIR})
else()
    message(WARNING "GSL not found – install libgsl-dev or set MACRODR_GSL_DIR")
    if(DEFINED ENV{MACRODR_GSL_DIR})
        set(_macrodr_gsl_root $ENV{MACRODR_GSL_DIR})
        set(_macrodr_gsl_lib "${_macrodr_gsl_root}/lib/libgsl.so")
        if(EXISTS ${_macrodr_gsl_lib})
            set(gsl_library ${_macrodr_gsl_lib})
            set(gsl_dir     "${_macrodr_gsl_root}/include")
        else()
            message(FATAL_ERROR "GSL not found at ${_macrodr_gsl_root}")
        endif()
    elseif(EXISTS "/nfs/home/lmoffatt/.nix-profile/lib/libgsl.so")
        message(STATUS "Using fallback Nix GSL")
        set(gsl_library /nfs/home/lmoffatt/.nix-profile/lib/libgsl.so)
        set(gsl_dir     /nix/store/snnl5y7hbk7b49cvfi8mhvgzrknjwz1p-gsl-2.7.1-dev/include/)
    else()
        message(FATAL_ERROR "GSL not found")
    endif()
endif()

# -------------------------------------------------------
# 5) Public include directory ### NEW
# -------------------------------------------------------
set(MACRODR_PUBLIC_INC ${PROJECT_SOURCE_DIR}/include)

include(CTest)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

# -------------------------------------------------------
## ---------- 6) CORE library ----------------------------------------
file(GLOB_RECURSE MACRODR_CORE_HEADERS CONFIGURE_DEPENDS
     ${PROJECT_SOURCE_DIR}/include/macrodr/core/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/core/*.hpp
     ${PROJECT_SOURCE_DIR}/include/macrodr/cmd/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/cmd/*.hpp
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/*.hpp
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/json/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/json/*.hpp
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/csv/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/io/csv/*.hpp
     ${PROJECT_SOURCE_DIR}/legacy/*.h
     ${PROJECT_SOURCE_DIR}/legacy/*.hpp
 )

file(GLOB_RECURSE MACRODR_CORE_SRC CONFIGURE_DEPENDS  src/core/*.cpp)

add_library(macrodr_core STATIC
    ${MACRODR_PUBLIC_HEADERS}
    ${MACRODR_CORE_HEADERS}
    ${MACRODR_CORE_SRC}
    ${MACRODR_LEGACY_HEADERS})

target_include_directories(macrodr_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_include_directories(macrodr_core SYSTEM PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/legacy>
    $<INSTALL_INTERFACE:legacy>)

target_compile_definitions(macrodr_core
    PUBLIC GIT_COMMIT_HASH="${GIT_COMMIT_HASH}")

# -------------------------------------------------------
# 7) Legacy object library (kept compiling) ### NEW
# -------------------------------------------------------
# immediately after the block that defines macrodr_core
target_link_libraries(macrodr_core
    PUBLIC
        OpenMP::OpenMP_CXX          # from find_package(OpenMP)
        ${lapackblas_libraries}     # BLAS/LAPACK combo you built above
        ${gsl_library}              # GSL (imported target or path)
)

 
 
# -------------------------------------------------------
# ---------- 8) CLI executable --------------------------------------
file(GLOB_RECURSE MACRODR_CLI_SRC      CONFIGURE_DEPENDS  src/cli/*.cpp)
file(GLOB_RECURSE MACRODR_CLI_HEADERS  CONFIGURE_DEPENDS
     ${PROJECT_SOURCE_DIR}/include/macrodr/cli/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/dsl/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/interface/*.h
     ${PROJECT_SOURCE_DIR}/include/macrodr/cmd/*.h)

list(REMOVE_ITEM MACRODR_CLI_SRC ${PROJECT_SOURCE_DIR}/src/cli/main_helpers.cpp)

add_library(macrodr_cli_helpers STATIC
    src/cli/main_helpers.cpp)

target_link_libraries(macrodr_cli_helpers PUBLIC macrodr_core)
target_include_directories(macrodr_cli_helpers PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

 add_executable(macrodr_cli
     ${MACRODR_CLI_SRC}
     ${MACRODR_CLI_HEADERS})

target_include_directories(macrodr_cli PRIVATE
    ${PROJECT_SOURCE_DIR}/include/macrodr/cli
    ${PROJECT_SOURCE_DIR}/include/macrodr/dsl
    ${PROJECT_SOURCE_DIR}/include/macrodr/interface
    ${PROJECT_SOURCE_DIR}/include/macrodr/cmd)
target_include_directories(macrodr_cli SYSTEM PRIVATE
    ${PROJECT_SOURCE_DIR}/legacy)

target_link_libraries(macrodr_cli PRIVATE macrodr_core macrodr_cli_helpers)

# -------------------------------------------------------
# 9) Data-file groups for Qt Creator (unchanged)
# -------------------------------------------------------
set(DATA_DIRS scripts experiments models slurm clusters docs multi_task spec MacroIR)
set(PROJECT_DATA_FILES)
foreach(dir IN LISTS DATA_DIRS)
    file(GLOB_RECURSE _tmp CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/${dir}/*.macroir"
        "${CMAKE_SOURCE_DIR}/${dir}/*.md"
         "${CMAKE_SOURCE_DIR}/${dir}/*.txt"
         "${CMAKE_SOURCE_DIR}/${dir}/*.sh"
         "${CMAKE_SOURCE_DIR}/${dir}/*.csv")
    list(APPEND PROJECT_DATA_FILES ${_tmp})
endforeach()
source_group("Headers" TREE ${CMAKE_SOURCE_DIR} FILES
    ${MACRODR_CORE_HEADERS} ${MACRODR_CLI_HEADERS})
source_group("Sources"     TREE ${CMAKE_SOURCE_DIR} FILES ${MACRODR_CORE_SRC} ${MACRODR_LEGACY_SRC})
source_group("Data Files"  TREE ${CMAKE_SOURCE_DIR} FILES ${PROJECT_DATA_FILES})
add_custom_target(data_files_group SOURCES ${PROJECT_DATA_FILES})

# -------------------------------------------------------
# 10) Install rules (binary + headers)
# -------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS macrodr_core ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS macrodr_cli  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${MACRODR_PUBLIC_INC}/macrodr
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
message(STATUS "PROJECT_SOURCE_DIR is: ${PROJECT_SOURCE_DIR}")
## Packaging/export of CMake config removed for now (was untested).

