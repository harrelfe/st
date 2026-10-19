#!/bin/zsh
# check-dois-append-only.sh -- git pre-commit check for the master DOI file.
# Rejects any commit in which dois.tsv is deleted, or in which any byte of the
# previously committed content was modified, removed, or reordered.  Only appending is allowed.
#
# Install (once, from the project root; skip if you already have a pre-commit hook):
#   printf '#!/bin/zsh\nexec "$(git rev-parse --show-toplevel)/check-dois-append-only.sh"\n' > .git/hooks/pre-commit
#   chmod +x .git/hooks/pre-commit check-dois-append-only.sh
# (git commit --no-verify bypasses any hook; this guards against accidents, not intent.)

f=dois.tsv
cd "$(git rev-parse --show-toplevel)" || exit 1

# Nothing to check unless the file exists in HEAD (first commit of it is always fine).
git cat-file -e "HEAD:$f" 2>/dev/null || exit 0

# Only check commits that touch the file.
git diff --cached --name-only | grep -qx "$f" || exit 0

if ! git cat-file -e ":$f" 2>/dev/null; then
  echo "REFUSED: $f is being deleted.  The DOI master file must never be removed." >&2
  exit 1
fi

old=$(mktemp "${TMPDIR:-/tmp}/dois-old.XXXXXX") ; new=$(mktemp "${TMPDIR:-/tmp}/dois-new.XXXXXX")
trap 'rm -f "$old" "$new"' EXIT
git show "HEAD:$f" > "$old"
git show ":$f"     > "$new"
n=$(wc -c < "$old" | tr -d ' ')       # macOS wc pads its output with spaces

if ! head -c "$n" "$new" | cmp -s - "$old"; then
  echo "REFUSED: the committed content of $f was changed.  Only appending new lines is allowed." >&2
  echo "  Compare:  git diff --cached -- $f" >&2
  echo "  Undo:     git restore --staged --worktree -- $f" >&2
  exit 1
fi
exit 0
