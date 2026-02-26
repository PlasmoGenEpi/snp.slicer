# Snp.slicer package — devtools workflow
# Run from package root: make <target>

.PHONY: document check test build install clean

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

# Full prep: document, then check (recommended before release)
prep: document check

# Clean check artifacts
clean:
	rm -rf ..Rcheck
	rm -f *.tar.gz
