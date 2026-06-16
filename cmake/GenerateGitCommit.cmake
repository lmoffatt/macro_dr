# Run at BUILD time (cmake -P) to regenerate the git-commit stamp TU. Queries git,
# then configure_file()s the template into OUTPUT_FILE. configure_file only rewrites
# the output when the rendered content differs (write-if-different), so the generated
# .cpp keeps a stable mtime and the compiler re-runs ONLY when the hash changes.
#
# Inputs (passed via -D on the cmake -P invocation):
#   GIT_SRC_DIR    repo root (where .git lives)
#   TEMPLATE_FILE  cmake/git_commit.cpp.in
#   OUTPUT_FILE    <build>/git_commit.generated.cpp

# OVERRIDE_HASH (optional): an explicit stamp passed by the caller. Used verbatim,
# skipping git entirely — build_cluster.sh sets it to the login-node hash because
# the compile runs on a compute node where git is usually absent (which otherwise
# makes the query below fail silently and fall back to "unknown").
if(NOT "${OVERRIDE_HASH}" STREQUAL "")
    set(MACRODR_GIT_COMMIT_HASH "${OVERRIDE_HASH}")
else()
    set(MACRODR_GIT_COMMIT_HASH "")
    if(EXISTS "${GIT_SRC_DIR}/.git")
        execute_process(
            COMMAND git rev-parse --short HEAD
            WORKING_DIRECTORY "${GIT_SRC_DIR}"
            OUTPUT_VARIABLE MACRODR_GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        # "-dirty" when there are uncommitted tracked changes.
        execute_process(
            COMMAND git status --porcelain --untracked-files=no
            WORKING_DIRECTORY "${GIT_SRC_DIR}"
            OUTPUT_VARIABLE MACRODR_GIT_DIRTY
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(NOT MACRODR_GIT_DIRTY STREQUAL "")
            set(MACRODR_GIT_COMMIT_HASH "${MACRODR_GIT_COMMIT_HASH}-dirty")
        endif()
    endif()
    if(MACRODR_GIT_COMMIT_HASH STREQUAL "")
        set(MACRODR_GIT_COMMIT_HASH "unknown")
    endif()
endif()

configure_file("${TEMPLATE_FILE}" "${OUTPUT_FILE}" @ONLY)
