# Contributing

## Branching model (GitFlow)

This project follows a GitFlow-style branching model (2026-09-15 onward, reverting to the
convention used through `2024.03.0` after a period of committing directly to `main`).

- **`main`** — always reflects the latest released state. Nothing is committed here directly;
  it only receives merges from `release/*` and `hotfix/*` branches, each merge tagged with the
  version number (e.g. `2026.9.0`, no `v` prefix).
- **`develop`** — the integration branch for ongoing development. Day-to-day work targets this
  branch, not `main`.
- **`feature/*`** — branch off `develop` for a unit of work, merge back into `develop` once
  reviewed.
- **`release/X.Y.Z`** — branched from `develop` when preparing a release. Final adjustments
  (version bump, changelog, etc.) happen here, then it merges into both `main` (tagged) and
  `develop`.
- **`hotfix/*`** — branched from `main` for urgent fixes to an already-released version. Merges
  into both `main` (tagged) and `develop`.

Versioning follows `YYYY.M.PATCH` (calendar versioning): the first release in a given
year/month is `PATCH=0` (e.g. `2026.9.0`), subsequent releases in the same month increment
`PATCH` (e.g. `2026.9.1`). Keep `proteindf_bridge/_version.py` and the `version` field in
`setup.cfg` in sync.

See `docs/rust-port-handoff.md` for the Rust port (`rust/`) subproject's specific PR/review
workflow, which follows this same branching model.
