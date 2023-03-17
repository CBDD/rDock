#!/bin/bash

if [ $PY3K -eq 1 ]; then
    2to3 -w -n build/test/RBT_HOME/check_test.py
fi

cd build/
chmod a+x p4utils.pl
make linux-g++-64
make test

cd ..
cp lib/libRbt.* "${PREFIX}/lib/"
mkdir -p "${PREFIX}/lib/perl5/site_perl/" "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/bin"
cp lib/*.pl lib/*.pm "${PREFIX}/lib/perl5/site_perl/"
mv data/ "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/"
mkdir -p "${PREFIX}/bin/"
# Create wrappers for binaries that need RBT_ROOT to be in the environment
for f in rbcalcgrid rbcavity rbdock rblist rbmoegrid; do
    mv "bin/$f" "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/bin/"
    sed -e "s|CHANGEME|${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}|" "$RECIPE_DIR/wrapper.sh" > "${PREFIX}/bin/$f"
    chmod +x "${PREFIX}/bin/$f"
done
# Remove unused to_unix
rm bin/to_unix
cp bin/* "${PREFIX}/bin/"