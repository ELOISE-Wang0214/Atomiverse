#!/usr/bin/env bash
# Public CLI contract regression (python-on-demand-install-04).
# No full native compile — help / mutual exclusion / --hf / version mismatch.
set -euo pipefail

PYTHON="${1:-python3}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
SRC="$REPO/python/src"
export PYTHONPATH="$SRC${PYTHONPATH:+:$PYTHONPATH}"

echo "== install_cli regression =="
echo "  python: $PYTHON"
echo "  src:    $SRC"

# help
OUT="$("$PYTHON" -m atomiverse install --help 2>&1)" || {
    echo "FAIL: install --help failed"; echo "$OUT"; exit 1
}
echo "$OUT" | grep -q -- '--config' || { echo "FAIL: help missing --config"; exit 1; }
echo "$OUT" | grep -q -- '--src' || { echo "FAIL: help missing --src"; exit 1; }
echo "$OUT" | grep -q -- '--url' || { echo "FAIL: help missing --url"; exit 1; }
echo "$OUT" | grep -qi module || { echo "FAIL: help missing MODULE"; exit 1; }
echo "ok: help lists --config/--src/--url/MODULE"

# mutual exclusion (scheme B: positional module names)
if "$PYTHON" -m atomiverse install --config "$HERE/fixtures/minimal_ok.toml" md >/dev/null 2>&1; then
    echo "FAIL: expected mutual exclusion failure"
    exit 1
fi
echo "ok: --config + md rejected"

# --hf
if "$PYTHON" -m atomiverse install --hf >/dev/null 2>&1; then
    echo "FAIL: expected --hf rejection"
    exit 1
fi
echo "ok: --hf rejected"

# version mismatch fixture (toml only — also need src with wrong version)
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/python"
cat > "$TMP/python/pyproject.toml" <<'EOF'
[project]
name = "atomiverse"
version = "0.0.0-mismatch"
EOF
if "$PYTHON" -m atomiverse install md --src "$TMP" --dry-run >/dev/null 2>&1; then
    echo "FAIL: expected version mismatch rejection"
    exit 1
fi
echo "ok: version mismatch rejected"

echo "REGRESSION PASSED: python/install_cli"
