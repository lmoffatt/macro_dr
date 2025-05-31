# -------------------------------------------------------
# CMakeLists.txt for macro_dr (with folder hierarchy in Qt Creator)
# -------------------------------------------------------

cmake_minimum_required(VERSION 3.10)
project(macro_dr LANGUAGES CXX)

# -------------------------------------------------------
# 1) Set C++ standard and common compile options
# -------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    -ftemplate-backtrace-limit=20
    -fmax-errors=5
    -fconcepts-diagnostics-depth=4
    -Werror=return-type
    -ftime-report
)

# Configuration-specific flags
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# -------------------------------------------------------
# 2) Print CMAKE_SOURCE_DIR for verification (optional)
# -------------------------------------------------------
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")

# -------------------------------------------------------
# 3) Obtain GIT commit hash if .git exists
# -------------------------------------------------------
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
else()
    set(GIT_COMMIT_HASH "unknown")
endif()

# -------------------------------------------------------
# 4) Find dependencies: OpenMP, BLAS/LAPACK, GSL
#    If find_package fails, fall back to hardcoded Nix paths with warnings.
# -------------------------------------------------------

# 4.1) OpenMP
find_package(OpenMP REQUIRED)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP not found. Load a compiler/module with -fopenmp support.")
endif()

# 4.2) BLAS and LAPACK
find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    message(WARNING
        "find_package(BLAS/LAPACK) failed. "
        "Falling back to hardcoded Nix paths for BLAS and LAPACK."
    )
    set(lapackblas_libraries
        /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so
        /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so
    )
endif()

# 4.3) GSL
find_package(GSL)
if(GSL_FOUND)
    set(gsl_library ${GSL_LIBRARY})
    set(gsl_dir     ${GSL_INCLUDE_DIR})
else()
    message(WARNING
        "find_package(GSL) failed. "
        "Falling back to hardcoded Nix path for GSL."
    )
    set(gsl_library /nfs/home/lmoffatt/.nix-profile/lib/libgsl.so)
    set(gsl_dir     /nix/store/snnl5y7hbk7b49cvfi8mhvgzrknjwz1p-gsl-2.7.1-dev/include/)
endif()

# -------------------------------------------------------
# 5) Collect only the C++ source/header files (*.cpp, *.h)
#    from the entire project tree.
# -------------------------------------------------------
file(GLOB_RECURSE PROJECT_SOURCES
    "${CMAKE_SOURCE_DIR}/*.cpp"
    "${CMAKE_SOURCE_DIR}/*.h"
)

# Exclude any file that lives under test_macro_dr/ (so they won't be compiled)
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*/test_macro_dr/.*")
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*/tst_.*\\.cpp$")
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*/tst_.*\\.h$")


# (Optional) Print how many sources remain after filtering
message(STATUS "After filtering, PROJECT_SOURCES contains:")
foreach(src IN LISTS PROJECT_SOURCES)
    message(STATUS "  ${src}")
endforeach()

# Create a source group so IDEs show the folder hierarchy under "Sources"
source_group(
    "Sources"
    TREE ${CMAKE_SOURCE_DIR}
    FILES ${PROJECT_SOURCES}
)

# -------------------------------------------------------
# Define the executable using only PROJECT_SOURCES (excluding test_macro_dr)
# -------------------------------------------------------
add_executable(macro_dr
    ${PROJECT_SOURCES}
)




# -------------------------------------------------------
# 6) Collect “data” files (scripts, .txt, .csv, .sh) only from specific subdirectories
#    so that they appear in the IDE but are NOT compiled.
# -------------------------------------------------------
set(DATA_DIRS
    scripts
    experiments
    models
    slurm
    clusters
    docs
)

set(PROJECT_DATA_FILES)

foreach(dir IN LISTS DATA_DIRS)
    file(GLOB_RECURSE _txt_files
        "${CMAKE_SOURCE_DIR}/${dir}/*.txt"
    )
    file(GLOB_RECURSE _sh_files
        "${CMAKE_SOURCE_DIR}/${dir}/*.sh"
    )
    file(GLOB_RECURSE _csv_files
        "${CMAKE_SOURCE_DIR}/${dir}/*.csv"
    )
    file(GLOB_RECURSE _md_files
        "${CMAKE_SOURCE_DIR}/${dir}/*.md"
    )
    list(APPEND PROJECT_DATA_FILES
        ${_txt_files}
        ${_sh_files}
        ${_csv_files}
        ${_md_files}
    )
endforeach()

# Also include any .txt/.sh/.csv directly in the project root (optional)
##file(GLOB_RECURSE _root_txt_files
#    "${CMAKE_SOURCE_DIR}/*.txt"
#)
#file(GLOB_RECURSE _root_sh_files
#    "${CMAKE_SOURCE_DIR}/*.sh"
#)
#file(GLOB_RECURSE _root_csv_files
#    "${CMAKE_SOURCE_DIR}/*.csv"
#)
list(APPEND PROJECT_DATA_FILES
    ${_root_txt_files}
    ${_root_sh_files}
    ${_root_csv_files}
)

# -------------------------------------------------------
# 7) Create source groups with TREE so Qt Creator shows a folder hierarchy
#    ----------------------------------------------------------------
#  - “Sources” group will mirror folders containing your .cpp/.h files.
#  - “Data Files” group will mirror folders containing .txt/.sh/.csv.
# -------------------------------------------------------
# For all source files:
source_group(
    "Sources"
    TREE ${CMAKE_SOURCE_DIR}
    FILES ${PROJECT_SOURCES}
)

# For all data files:
source_group(
    "Data Files"
    TREE ${CMAKE_SOURCE_DIR}
    PREFIX "Data Files"
    FILES ${PROJECT_DATA_FILES}
)

# -------------------------------------------------------
# 8) Add a custom target so Qt Creator (and other IDEs) register the data files
#    This target has no build commands; it only exists so the IDE lists those files.
# -------------------------------------------------------
add_custom_target(data_files_group
    SOURCES ${PROJECT_DATA_FILES}
)


# Attach the GIT_COMMIT_HASH define to the target
target_compile_definitions(macro_dr
    PUBLIC GIT_COMMIT_HASH="${GIT_COMMIT_HASH}"
)

# -------------------------------------------------------
# 10) Include directories and link libraries for macro_dr
# -------------------------------------------------------
target_include_directories(macro_dr
    PRIVATE ${gsl_dir}
)

target_link_libraries(macro_dr
    PUBLIC ${lapackblas_libraries}
           ${gsl_library}
           OpenMP::OpenMP_CXX
)

# -------------------------------------------------------
# 11) Installation rules (optional)
# -------------------------------------------------------
include(GNUInstallDirs)

# Install the executable
install(TARGETS macro_dr
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDDIR}
)

# If you want `make install` to also copy scripts/experiments/models/slurm/clusters,
# uncomment these blocks:

#install(DIRECTORY scripts/
#    DESTINATION ${CMAKE_INSTALL_BINDDIR}/scripts
#    FILES_MATCHING PATTERN "*.txt" PATTERN "*.sh"
#)

#install(DIRECTORY experiments/
#    DESTINATION ${CMAKE_INSTALL_BINDDIR}/experiments
#    FILES_MATCHING PATTERN "*.txt" PATTERN "*.csv"
#)

#install(DIRECTORY models/
#    DESTINATION ${CMAKE_INSTALL_BINDDIR}/models
#    FILES_MATCHING PATTERN "*.txt" PATTERN "*.csv"
#)

#install(DIRECTORY slurm/
#    DESTINATION ${CMAKE_INSTALL_BINDDIR}/slurm
#    FILES_MATCHING PATTERN "*.sh"
#)

#install(DIRECTORY clusters/
#    DESTINATION ${CMAKE_INSTALL_BINDDIR}/clusters
#    FILES_MATCHING PATTERN "*.sh"
#)
