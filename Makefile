# Snp.slicer package — devtools workflow
# Run from package root: make <target>

.PHONY: document check test build install vignettes docs upgrade clean

# Regenerate NAMESPACE and man/ from roxygen2
document:
	R -e "devtools::document()"

# Run R CMD check (same as devtools::check())
check:
	R -e "devtools::check()"

# Run testthat tests
test:
	R -e "devtools::test()"

# Build source tarball
build:
	R -e "devtools::build()"

# Install package into the default library
install:
	R -e "devtools::install()"

# Build vignettes
vignettes:
	R -e "devtools::build_vignettes()"

# Build pkgdown site (output in docs/). Favicons from logo if missing.
docs:
	R -e "pkgdown::build_favicons(); pkgdown::build_site()"

# Upgrade all installed R packages (run from a session with writable library)
upgrade:
	R -e "update.packages(ask = FALSE, checkBuilt = TRUE)"

# Full prep: document, then check (recommended before release)
prep: document check

# Clean check artifacts
clean:
	rm -rf ..Rcheck
	rm -f *.tar.gz
