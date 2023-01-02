#!/bin/bash

cd build/
make linux-g++-64
make test

cd ..
cp lib/libRbt.so.rDock.0 "${PREFIX}/lib/libRbt.so.rDock"
PERL_INSTALLSITELIB=$(perl -e 'use Config; print "$Config{installsitelib}"')
mkdir -p "${PERL_INSTALLSITELIB}" "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/bin"
cp lib/*.pl lib/*.pm "${PERL_INSTALLSITELIB}"
mv data/ "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/"
# Create wrappers for binaries that need RBT_ROOT to be in the environment
for f in rbcalcgrid rbcavity rbdock rblist rbmoegrid; do
    mv "bin/$f" "${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}/bin/"
    sed -e "s|CHANGEME|${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}|" "$RECIPE_DIR/wrapper.sh" > "${PREFIX}/bin/$f"
    chmod +x "${PREFIX}/bin/$f"
done
# Remove unused to_unix
rm bin/to_unix
cp bin/* "${PREFIX}/bin/"