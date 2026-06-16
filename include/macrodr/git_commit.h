#pragma once

namespace macrodr {

// Git commit hash the binary was built from: the short hash, plus a "-dirty"
// suffix when the build tree had uncommitted tracked changes ("unknown" outside
// a git checkout). It is compiled into ONE build-time-GENERATED translation unit
// (git_commit.generated.cpp, regenerated write-if-different on every build), so a
// new commit recompiles only that single file (then a fast relink) — never the
// whole core library. Used as CSV-provenance row 1 and by the `--commit` CLI flag.
const char* git_commit_hash();

}  // namespace macrodr
