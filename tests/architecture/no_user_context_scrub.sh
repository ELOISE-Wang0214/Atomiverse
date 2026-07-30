#!/usr/bin/env bash
# Grep gate: user-facing trees never hand-drive the engine bus.
#
# Iron law (CLAUDE.md frame-only-token): the only public I/O token is Frame.
# Out-of-library consumers drive via presets (Evaluate / MD / Minimize / MC /
# single_point / method.run(frame) / HarmonicAnalysis::run / …).
#
# Banned in regressions / examples / app / docs/tutorial:
#   - constructing cpu::Context / LaunchContext / BuildContext
#   - #include of component/context.h
#   - hand Component::launch / brace-init build({ctx,f})
#   - protected Driver loop primitives (setup_pass / step_once / teardown_pass)
#   - constructing cpu::Runtime as a second token
#
# Allowed:
#   - Component author *overrides*: void launch(const cpu::LaunchContext&)
#   - CV build(const Frame&) — Frame-only collective-variable protocol
#   - preset .build() / .run(frame) / method.build(frame) public APIs
#   - unit tests under tests/; docs/developer/
#
# Usage: no_user_context_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: no_user_context_scrub.sh <repo-root>}"
SELF="no_user_context_scrub.sh"
status=0

# Construction / named variables of the engine bus (not override signatures).
CONSTRUCT_PATS=(
    'cpu::Context[[:space:]]+[A-Za-z_{]'
    'cuda::Context[[:space:]]+[A-Za-z_{]'
    'torch::Context[[:space:]]+[A-Za-z_{]'
    'cpu::LaunchContext[[:space:]]+[A-Za-z_{]'
    'cpu::BuildContext[[:space:]]+[A-Za-z_{]'
    'cpu::DestroyContext[[:space:]]+[A-Za-z_{]'
    'cuda::LaunchContext[[:space:]]+[A-Za-z_{]'
    'cuda::BuildContext[[:space:]]+[A-Za-z_{]'
    'ATV::cpu::Context[[:space:]]+[A-Za-z_{]'
    'ATV::cuda::Context[[:space:]]+[A-Za-z_{]'
    'cpu::Runtime[[:space:]]+[A-Za-z_{]'
    'ATV::cpu::Runtime[[:space:]]+[A-Za-z_{]'
)

# Hand-drive of protected / lifecycle entry points.
HAND_DRIVE_PATS=(
    'setup_pass[[:space:]]*\('
    'teardown_pass[[:space:]]*\('
    '\.step_once[[:space:]]*\('
    'step_once[[:space:]]*\('
    '\.launch[[:space:]]*\('
    '\.(build|destroy)[[:space:]]*\([[:space:]]*\{'
)

INCLUDE_PAT='atomiverse/cpu/component/context\.h|atomiverse/cuda/component/context\.cuh'

USER_DIRS=(
    "regressions"
    "examples"
    "app"
    "docs/tutorial"
)

for dir in "${USER_DIRS[@]}"; do
    target="${ROOT}/${dir}"
    [ -d "${target}" ] || continue
    for pat in "${CONSTRUCT_PATS[@]}"; do
        matches="$(grep -rnE --include='*.cpp' --include='*.h' --include='*.cuh' \
            --include='*.cu' --include='*.hpp' --include='*.md' \
            "${pat}" "${target}" 2>/dev/null \
            | grep -v "${SELF}" \
            | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*)' \
            | grep -vE 'Wrong:|banned|Never `|// wrong|forbidden' \
            || true)"
        if [ -n "${matches}" ]; then
            echo "FAIL: user-facing ${dir}/ constructs engine bus /${pat}/:"
            echo "${matches}" | head -40
            echo "  -> Frame + preset run/single_point only (iron law frame-only-token)."
            status=1
        fi
    done
    for pat in "${HAND_DRIVE_PATS[@]}"; do
        matches="$(grep -rnE --include='*.cpp' --include='*.h' --include='*.cuh' \
            --include='*.hpp' --include='*.md' \
            "${pat}" "${target}" 2>/dev/null \
            | grep -v "${SELF}" \
            | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*)' \
            | grep -vE 'Wrong:|banned|Never `|// wrong|forbidden' \
            | grep -vE 'void[[:space:]]+launch[[:space:]]*\(|override' \
            || true)"
        if [ -n "${matches}" ]; then
            echo "FAIL: user-facing ${dir}/ hand-drives lifecycle /${pat}/:"
            echo "${matches}" | head -40
            echo "  -> use Evaluate/MD/Minimize/Hessian/HarmonicAnalysis::run(frame)."
            status=1
        fi
    done
    matches="$(grep -rnE --include='*.cpp' --include='*.h' --include='*.cuh' \
        --include='*.hpp' --include='*.md' \
        "${INCLUDE_PAT}" "${target}" 2>/dev/null \
        | grep -v "${SELF}" \
        | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*)' \
        || true)"
    if [ -n "${matches}" ]; then
        echo "FAIL: user-facing ${dir}/ includes engine context header:"
        echo "${matches}" | head -40
        echo "  -> include component.h / evaluate.h / md.h; never context.h directly."
        status=1
    fi
done

if [ "${status}" -eq 0 ]; then
    echo "PASS: user trees free of Context construction and hand launch/runtime."
fi
exit "${status}"
