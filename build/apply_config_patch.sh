export PREFIX=$PWD
export C_INCLUDE_PATH=${PREFIX}/include
export CPP_INCLUDE_PATH=${PREFIX}/include
export CPLUS_INCLUDE_PATH=${PREFIX}/include
export CXX_INCLUDE_PATH=${PREFIX}/include
export LIBRARY_PATH=${PREFIX}/lib

patch -p1 < bioconda.patch

# patch from https://aur.archlinux.org/cgit/aur.git/tree/rdock-2021.1.patch?h=rdock
patch -p0 < gcc_compat.patch

# patch from my changes
patch -p1 < sruiz_personal.patch
# and notest, cppunit tests usually fail and are not useful
patch -p0 < notest.patch

# rcanovas's changes regarding FORBIDDEN error in residues
patch -p0 < gcc_rcanovas.patch
