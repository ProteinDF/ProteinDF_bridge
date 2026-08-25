# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment, for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile html html-ja html-all latexpdf latexpdf-ja gettext intl clean

# English HTML (default language)
html:
	@$(SPHINXBUILD) -b html -D language=en "$(SOURCEDIR)" "$(BUILDDIR)/html/en" $(SPHINXOPTS) $(O)

# Japanese HTML
html-ja:
	@$(SPHINXBUILD) -b html -D language=ja "$(SOURCEDIR)" "$(BUILDDIR)/html/ja" $(SPHINXOPTS) $(O)

# Both languages
html-all: html html-ja

# English PDF (requires a TeX installation, e.g. `brew install --cask mactex-no-gui`)
latexpdf:
	@$(SPHINXBUILD) -b latex -D language=en "$(SOURCEDIR)" "$(BUILDDIR)/latex/en" $(SPHINXOPTS) $(O)
	$(MAKE) -C "$(BUILDDIR)/latex/en" all-pdf

# Japanese PDF (requires texlive-lang-japanese / uplatex)
latexpdf-ja:
	@$(SPHINXBUILD) -b latex -D language=ja "$(SOURCEDIR)" "$(BUILDDIR)/latex/ja" $(SPHINXOPTS) $(O)
	$(MAKE) -C "$(BUILDDIR)/latex/ja" all-pdf

# Extract translatable strings into a gettext catalog
gettext:
	@$(SPHINXBUILD) -b gettext "$(SOURCEDIR)" "$(BUILDDIR)/gettext" $(SPHINXOPTS) $(O)

# Update the Japanese .po catalogs from the extracted .pot files
intl: gettext
	sphinx-intl update -p "$(BUILDDIR)/gettext" -l ja -d locale

clean:
	rm -rf "$(BUILDDIR)"

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
