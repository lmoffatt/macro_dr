# ──────────────────────────────────────────────────────────────────────
#  Macro DR – unified top-level CMakeLists.txt
# ──────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.18)
project(macro_dr LANGUAGES CXX)

# -------------------------------------------------------
# 1) Global compile rules              (kept from legacy)
# -------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(
    -ftemplate-backtrace-limit=20     # 
    -fmax-errors=5
    -fconcepts-diagnostics-depth=4
    -Werror=return-type
    -ftime-report
)
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

install(TARGETS macrodr_thirdparty EXPORT macrodrTargets)

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
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
    message(WARNING "BLAS/LAPACK not found – falling back to Nix store paths")
    set(lapackblas_libraries
        /nfs/home/lmoffatt/.nix-profile/lib/libcblas.so
        /nfs/home/lmoffatt/.nix-profile/lib/liblapack.so )      # 
endif()

find_package(GSL)
if(GSL_FOUND)
    set(gsl_library ${GSL_LIBRARY})
    set(gsl_dir     ${GSL_INCLUDE_DIR})
else()
    message(WARNING "GSL not found – falling back to Nix store path")
    set(gsl_library /nfs/home/lmoffatt/.nix-profile/lib/libgsl.so)
    set(gsl_dir     /nix/store/snnl5y7hbk7b49cvfi8mhvgzrknjwz1p-gsl-2.7.1-dev/include/)   # 
endif()

# -------------------------------------------------------
# 5) Public include directory ### NEW
# -------------------------------------------------------
set(MACRODR_PUBLIC_INC ${PROJECT_SOURCE_DIR}/include)

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
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/legacy>
    $<INSTALL_INTERFACE:include>
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

 add_executable(macrodr_cli
     ${MACRODR_CLI_SRC}
     ${MACRODR_CLI_HEADERS})

target_include_directories(macrodr_cli PRIVATE
    ${PROJECT_SOURCE_DIR}/include/macrodr/cli
    ${PROJECT_SOURCE_DIR}/include/macrodr/dsl
    ${PROJECT_SOURCE_DIR}/include/macrodr/interface
    ${PROJECT_SOURCE_DIR}/include/macrodr/cmd)

target_link_libraries(macrodr_cli PRIVATE macrodr_core)

# -------------------------------------------------------
# 9) Data-file groups for Qt Creator (unchanged)
# -------------------------------------------------------
set(DATA_DIRS scripts experiments models slurm clusters docs multi_task spec)
set(PROJECT_DATA_FILES)
foreach(dir IN LISTS DATA_DIRS)
    file(GLOB_RECURSE _tmp CONFIGURE_DEPENDS
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
# -------------------------------------------------------
# 11) Export config so downstream projects can find_package(macrodr)
# -------------------------------------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/macrodrConfigVersion.cmake
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/macrodrConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/macrodrConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/macrodr)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/macrodrConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/macrodrConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/macrodr)

#  ───  install the core library  ───────────────────────────
install(TARGETS macrodr_core
            EXPORT  macrodrTargets          # <- make sure this is present
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    
# (optional) also export the CLI if you want users to link it:
# install(TARGETS macrodr_cli
#         EXPORT  macrodrTargets
#         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(EXPORT macrodrTargets
        NAMESPACE macrodr::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/macrodr)
