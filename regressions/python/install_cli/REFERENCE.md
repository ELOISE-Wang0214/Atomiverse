# install_cli (docs / primary CLI contract)

Literature: N/A.

Expected exit semantics (no full native compile required):

| Case | Exit | Notes |
|------|------|--------|
| `install --help` | 0 | Mentions `--config`, `--md`, `--global`, `--src`, `--url` |
| `--config` + `--md` together | ≠0 | Mutual exclusion |
| `--hf` | ≠0 | SCF/HF bindings not available |
| Source pyproject version ≠ installed | ≠0 | Before pip/compile |

Canonical product docs: `python/docs/getting_started.md`.
