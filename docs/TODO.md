# Documentation TODO

- [x] Enable GitHub Pages for this repo (Settings → Pages → Source: GitHub Actions)
      so the `docs` workflow can deploy. (Verified 2026-09-08: already enabled —
      `build_type: "workflow"`, live at http://proteindf.github.io/ProteinDF_bridge/,
      last successful `docs` workflow run 2026-08-25. This item was stale.)
- [x] Translate the API reference into Japanese
      (`docs/locale/ja/LC_MESSAGES/api/*.po`). Regenerated the catalogs
      (`make -C docs intl`) against the current source first, then
      translated all 107 previously-empty entries across the 32 `api/*.po`
      files (module titles, Napoleon `Parameters`/`Returns`/`Return type`
      fields, docstring prose; `:sphinx_autodoc_typehints_type:` markup and
      literal ASCII-diagram lines were copied verbatim rather than
      translated, since they're not natural language). Verified with
      `sphinx-build -b html -D language=ja docs <out>` — builds clean and
      renders the translated strings correctly.
- [x] Fix the docstring formatting warning in `proteindf_bridge/gro.py`
      (`SimpleGro` docstring: block quote ends without a blank line).
      Turned the inconsistently-indented sample into a proper reST
      literal block (`sample::`); verified the warning is gone via
      `sphinx-build -b html -D language=en docs /tmp/docsbuild_check`.
      Found and fixed the same "Unexpected indentation" docutils error in
      `proteindf_bridge/modeling.py` (`Modeling.get_ACE`/`get_NME`
      docstrings — an ACE/NME diagram with the same indentation issue),
      same fix (turned into a literal block).
- [ ] Install a local TeX distribution to build PDFs
      (`brew install --cask mactex-no-gui`, or `texlive-lang-japanese` on
      Linux for the Japanese PDF) and verify `make latexpdf` /
      `make latexpdf-ja`.
- [ ] Review the auto-generated API reference (`docs/api/*.rst`) for modules
      that should be excluded or reorganized (e.g. `dbmanager.py`, `mail.py`
      if not part of the public API).
