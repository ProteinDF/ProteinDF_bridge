# Documentation TODO

- [x] Enable GitHub Pages for this repo (Settings → Pages → Source: GitHub Actions)
      so the `docs` workflow can deploy. (Verified 2026-09-08: already enabled —
      `build_type: "workflow"`, live at http://proteindf.github.io/ProteinDF_bridge/,
      last successful `docs` workflow run 2026-08-25. This item was stale.)
- [ ] Translate the API reference into Japanese
      (`docs/locale/ja/LC_MESSAGES/api/*.po`) — currently untranslated and
      falls back to English.
- [ ] Fix the docstring formatting warning in `proteindf_bridge/gro.py`
      (`SimpleGro` docstring: block quote ends without a blank line).
- [ ] Install a local TeX distribution to build PDFs
      (`brew install --cask mactex-no-gui`, or `texlive-lang-japanese` on
      Linux for the Japanese PDF) and verify `make latexpdf` /
      `make latexpdf-ja`.
- [ ] Review the auto-generated API reference (`docs/api/*.rst`) for modules
      that should be excluded or reorganized (e.g. `dbmanager.py`, `mail.py`
      if not part of the public API).
