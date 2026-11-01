#!/bin/bash
# pixi does not run conda post-link scripts by default. Bioconductor data
# packages (GenomeInfoDbData, GO.db, org.Hs.eg.db, HDO.db) install their
# contents in such a script (download + md5 check + R CMD INSTALL), so run
# those here. Already installed packages are skipped.
set -e
for f in "$CONDA_PREFIX"/bin/.bioconductor-*-post-link.sh; do
  [ -e "$f" ] || continue
  pkg=$(sed -n 's/^installBiocDataPackage.sh "\(.*\)-[0-9.]*"$/\1/p' "$f")
  if ls "$CONDA_PREFIX/lib/R/library" | grep -qix "$pkg"; then
    echo "$pkg already installed"
  else
    PREFIX="$CONDA_PREFIX" bash "$f"
  fi
done
