#!/usr/bin/env python3
"""check.sh item 12: printed sentences whose opening words were swallowed by a % comment.

The failure this catches has happened FIVE times in this manuscript and it is invisible to the
compiler: an author appends a note to the end of a prose line, or wraps a comment onto the line a
sentence starts on, and the words before the % vanish from the PDF. What is left renders as a
sentence with no subject, so the paper reads as if it had a typo where it actually lost a clause.

The detector strips comments the way TeX does (a whole-line comment prints nothing AND does not
break a paragraph, so it is deleted rather than blanked; a trailing comment takes the rest of its
line), then reports two signatures in what remains: a paragraph that begins with a lowercase word,
and a sentence that begins with a lowercase word after a full stop.

Usage: swallowed_sentences.py <dir-of-tex-files>
Prints one line per hit and a final "N swallowed sentence(s)" summary; exit status is always 0, the
caller decides severity.
"""
import re, sys, glob, os

ABBREV = {'et', 'al', 'vs', 'cf', 'eg', 'ie'}

def strip_comments(text: str) -> str:
    kept = []
    for line in text.split('\n'):
        if line.lstrip().startswith('%'):
            continue
        kept.append(re.sub(r'(?<!\\)%.*$', '', line))
    return '\n'.join(kept)

def swallowed_whole(path: str):
    """Prose commented out entire, which leaves no lowercase signature.

    Found on 2026-08-26, when four sentences of the Discussion, the P2X2 verdict among them, turned
    out to have been invisible since the day they were written. The manuscript writes one paragraph
    per line, so a multi-line comment inserted before a paragraph appends the whole rest of that
    paragraph to its last line. The signature is a comment line far longer than the file's own
    comment style that still carries prose markup.
    """
    hits = []
    comment_lengths = [len(l) for l in open(path).read().split('\n') if l.lstrip().startswith('%')]
    if not comment_lengths:
        return hits
    typical = sorted(comment_lengths)[len(comment_lengths) // 2]
    limit = max(3 * typical, 250)
    for n, line in enumerate(open(path).read().split('\n'), 1):
        if not line.lstrip().startswith('%') or len(line) <= limit:
            continue
        if not re.search(r'\\(cite[a-z]*|ref|texttt|emph|citep|citet)\{', line):
            continue
        hits.append((n, len(line), line[:70]))
    return hits


def scan(path: str):
    text = strip_comments(open(path).read())
    hits = []
    for m in re.finditer(r'\n[ \t]*\n[ \t]*([a-z][a-z]+)\b', text):
        hits.append((m.group(1), text[m.start():m.start() + 90].strip().replace('\n', ' ')))
    for m in re.finditer(r'(?<![A-Z])\.\s+([a-z][a-z]{2,})\b', text):
        if m.group(1) in ABBREV:
            continue
        hits.append((m.group(1), text[max(0, m.start() - 50):m.start() + 60].strip().replace('\n', ' ')))
    return hits

def main() -> int:
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    total = 0
    lines = []
    for f in sorted(glob.glob(os.path.join(root, '*.tex'))):
        for word, ctx in scan(f):
            lines.append("%s: sentence starts '%s' ... %s" % (os.path.basename(f), word, ctx))
            total += 1
        for n, length, head in swallowed_whole(f):
            lines.append("%s:%d: %d-char comment line carries prose markup ... %s"
                         % (os.path.basename(f), n, length, head))
            total += 1
    print('%d swallowed sentence(s)' % total)
    for line in lines:
        print(line)
    return 0

if __name__ == '__main__':
    sys.exit(main())
