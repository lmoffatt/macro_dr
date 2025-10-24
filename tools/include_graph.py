#!/usr/bin/env python3
"""
Generate a C/C++ include graph (directed: file -> included file) from
compile_commands.json using the compiler's -E -H output. Can print cycles and
emit a Graphviz .dot you can preview inside VS Code with the Graphviz extension.

Typical usage:
  - Project graph (headers only):
      python3 tools/include_graph.py --dot build/include-graph.dot
  - Focus on one file's outgoing includes (transitively):
      python3 tools/include_graph.py --focus legacy/variables.h \
                                     --dot build/variables.dot
  - Show detected header cycles:
      python3 tools/include_graph.py --cycles

Requires a valid compile_commands.json at repo root (the repo already has a
symlink compile_commands.json -> build/gcc-debug/compile_commands.json).
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from collections import defaultdict, deque
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


HeaderExts = {".h", ".hpp", ".hh", ".hxx", ".inc"}


def is_header(path: str) -> bool:
    return os.path.splitext(path)[1].lower() in HeaderExts


def read_compile_commands(path: str) -> List[dict]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def real_is_child(p: str, root: str) -> bool:
    rp = os.path.realpath(p)
    rr = os.path.realpath(root)
    return rp.startswith(rr + os.sep) or rp == rr


def gather_edges(
    compdb: List[dict],
    cwd_root: str,
    *,
    local_only: bool = True,
    headers_only: bool = True,
    verbose: bool = False,
) -> Tuple[Set[Tuple[str, str]], Set[str]]:
    """Return (edges, nodes) where edges are (including -> included).

    Only includes edges where both endpoints are inside the repo when
    local_only=True, and optionally filters to header->header edges.
    Nodes contains all file paths (relative to repo root) observed in edges.
    """
    edges: Set[Tuple[str, str]] = set()
    nodes: Set[str] = set()

    for entry in compdb:
        # Prefer 'command' but fall back to 'arguments'.
        if entry.get("command"):
            args = shlex.split(entry["command"])  # preserves quoting
        else:
            args = entry.get("arguments", [])
        if not args:
            continue
        compiler, *toks = args
        # Rebuild command for preprocessing: strip -c/-o to avoid object writes.
        new_args = [compiler]
        i = 0
        while i < len(toks):
            t = toks[i]
            # Remove output flag pairs and standalone -c
            if t == "-o" and i + 1 < len(toks):
                i += 2
                continue
            if t == "-c":
                i += 1
                continue
            if t.startswith("-o"):
                i += 1
                continue
            new_args.append(t)
            i += 1

        src = entry["file"]
        new_args.extend(["-E", "-H", src])
        try:
            proc = subprocess.run(
                new_args,
                cwd=entry["directory"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        except FileNotFoundError:
            # Compiler not present? Skip.
            if verbose:
                print(f"warn: compiler not found for {src}")
            continue

        # GCC/Clang print one line per include with dot depth and absolute path
        stack: List[str] = [os.path.realpath(src)]
        for line in proc.stderr.splitlines():
            s = line.lstrip(".")
            depth = len(line) - len(s)
            s = s.strip()
            if depth == 0 or not s:
                continue
            # Skip lines that aren't file paths
            if not (s.startswith("/") or (":" in s and "/" in s)):
                continue
            inc_to = os.path.realpath(s.replace("\\", "/"))
            # Adjust stack to depth
            while len(stack) > depth:
                stack.pop()
            if not stack:
                stack = [os.path.realpath(src)]
            inc_by = stack[-1]
            stack.append(inc_to)

            if local_only and not (real_is_child(inc_by, cwd_root) and real_is_child(inc_to, cwd_root)):
                continue
            a = os.path.relpath(inc_by, cwd_root)
            b = os.path.relpath(inc_to, cwd_root)
            if headers_only and not (is_header(a) and is_header(b)):
                continue
            edges.add((a, b))
            nodes.add(a)
            nodes.add(b)

    return edges, nodes


def build_adj(edges: Iterable[Tuple[str, str]]) -> Dict[str, Set[str]]:
    adj: Dict[str, Set[str]] = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
    return adj


def find_cycles(adj: Dict[str, Set[str]]) -> List[List[str]]:
    visited: Set[str] = set()
    in_stack: Set[str] = set()
    stack: List[str] = []
    raw_cycles: List[List[str]] = []

    def dfs(u: str) -> None:
        visited.add(u)
        in_stack.add(u)
        stack.append(u)
        for v in adj.get(u, ()): 
            if v not in visited:
                dfs(v)
            elif v in in_stack:
                # record cycle u->...->v
                if v in stack:
                    i = stack.index(v)
                    raw_cycles.append(stack[i:] + [v])
        stack.pop()
        in_stack.remove(u)

    for n in list(adj.keys()):
        if n not in visited:
            dfs(n)

    # deduplicate by canonical rotation (and reverse)
    uniq: List[List[str]] = []
    seen: Set[Tuple[str, ...]] = set()
    for cyc in raw_cycles:
        cyc = cyc[:-1]
        if not cyc:
            continue
        rots = [tuple(cyc[i:] + cyc[:i]) for i in range(len(cyc))]
        rc = min(rots)
        rcyc = list(reversed(cyc))
        rots_r = [tuple(rcyc[i:] + rcyc[:i]) for i in range(len(rcyc))]
        key = min(rc, min(rots_r))
        if key not in seen:
            seen.add(key)
            uniq.append(list(key))
    return uniq


def reachable(
    adj: Dict[str, Set[str]],
    start: str,
    direction: str = "out",
) -> Set[str]:
    """Nodes reachable from start following edges in given direction.
    direction in {"out", "in", "both"}.
    """
    if start not in adj and direction == "out":
        # still include the start node
        return {start}
    fwd = adj
    rev: Dict[str, Set[str]] = defaultdict(set)
    for a, nbrs in adj.items():
        for b in nbrs:
            rev[b].add(a)
    seen: Set[str] = {start}
    dq: deque[str] = deque([start])
    while dq:
        u = dq.popleft()
        nexts: Iterable[str] = []
        if direction in ("out", "both"):
            nexts = list(fwd.get(u, ()))
        if direction in ("in", "both"):
            nexts = list(set(nexts) | rev.get(u, set()))
        for v in nexts:
            if v not in seen:
                seen.add(v)
                dq.append(v)
    return seen


def write_dot(
    edges: Iterable[Tuple[str, str]],
    nodes: Iterable[str],
    *,
    root: str,
    focus: Optional[str] = None,
    outfile: str,
) -> None:
    # Styling by top-level folder
    def style_for(n: str) -> Tuple[str, str]:
        if n.startswith("include/"):
            return ("#1f7a1f", "box")  # green
        if n.startswith("legacy/"):
            return ("#d9901a", "box")  # orange
        if n.startswith("src/"):
            return ("#1a73e8", "ellipse")  # blue
        if n.startswith("third_party/"):
            return ("#777777", "box")  # grey
        return ("#333333", "box")

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("digraph G {\n")
        f.write("  rankdir=LR;\n")
        f.write("  node [fontsize=10];\n")
        # Nodes
        for n in sorted(set(nodes)):
            color, shape = style_for(n)
            label = n
            if "/" in n:
                label = n  # keep rel path label
            url = n  # VS Code Graphviz preview will display label; URL often ignored
            f.write(
                f'  "{n}" [label="{label}", color="{color}", shape="{shape}", tooltip="{n}"];\n'
            )
        # Edges
        for a, b in sorted(set(edges)):
            f.write(f'  "{a}" -> "{b}";\n')
        # Emphasize focus node visually
        if focus is not None:
            f.write(f'  "{focus}" [style="filled", fillcolor="#ffe082"];\n')
        f.write("}\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Generate C/C++ include graph from compile_commands.json")
    ap.add_argument(
        "--compile-commands",
        default="compile_commands.json",
        help="Path to compile_commands.json (default: compile_commands.json)",
    )
    ap.add_argument("--dot", help="Write Graphviz .dot to this path (optional)")
    ap.add_argument(
        "--focus",
        help="Focus graph on this file (relative path). Shows transitive outgoing includes.",
    )
    ap.add_argument(
        "--direction",
        choices=["out", "in", "both"],
        default="out",
        help="Traversal direction for --focus (default: out)",
    )
    ap.add_argument(
        "--all-files",
        action="store_true",
        help="Include edges for all files, not just header->header",
    )
    ap.add_argument(
        "--include-external",
        action="store_true",
        help="Include edges to system/external headers (default filters to repo files)",
    )
    ap.add_argument(
        "--cycles",
        action="store_true",
        help="Print detected header include cycles to stdout",
    )
    ap.add_argument(
        "--verbose", "-v", action="store_true", help="Verbose logging while scanning"
    )

    args = ap.parse_args(argv)
    root = os.getcwd()
    compdb = read_compile_commands(args.compile_commands)
    edges, nodes = gather_edges(
        compdb,
        root,
        local_only=not args.include_external,
        headers_only=not args.all_files,
        verbose=args.verbose,
    )
    if args.focus:
        focus_rel = os.path.relpath(os.path.realpath(args.focus), root)
        adj = build_adj(edges)
        sub_nodes = reachable(adj, focus_rel, direction=args.direction)
        sub_edges = set((a, b) for a, b in edges if a in sub_nodes and b in sub_nodes)
        nodes = sub_nodes
        edges = sub_edges
        focus = focus_rel
    else:
        focus = None

    if args.cycles:
        cycles = find_cycles(build_adj(edges))
        if not cycles:
            print("No header include cycles found.")
        else:
            print(f"Found {len(cycles)} header include cycle(s).\n")
            for i, cyc in enumerate(cycles, 1):
                print(f"Cycle {i}:")
                for a, b in zip(cyc, cyc[1:] + [cyc[0]]):
                    print(f"  {a} -> {b}")
                print()

    if args.dot:
        write_dot(edges, nodes, root=root, focus=focus, outfile=args.dot)
        if args.verbose:
            print(f"Wrote {args.dot}")

    # Exit code indicates success; not printing edges by default
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

