import subprocess
from pytest import fixture
from distutils import dir_util
import os


@fixture
def datadir(tmpdir, request):
    """
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def _write_gridfile(datadir, gridfile):
    with open(gridfile, "w") as f:
        f.write(
            f"""RBT_PARAMETER_FILE_V1.00
TITLE TEST_STEFAN

RECEPTOR_FILE {datadir}/protein_prepared.mol2
#RECEPTOR_FLEX 3.0

SECTION LIGAND
	TRANS_MODE FREE
	ROT_MODE FREE
        DIHEDRAL_MODE FREE
END_SECTION


##################################################################
### CAVITY DEFINITION: TWO SPHERE METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {datadir}/ref_lig.sdf
    RADIUS 8
    SMALL_SPHERE 1.5
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION



#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION


#################################
## PHARMACOPHORIC RESTRAINTS
#################################
SECTION PHARMA
    SCORING_FUNCTION RbtPharmaSF
    WEIGHT 1.0
    CONSTRAINTS_FILE {datadir}/mandatory_restraints.txt
    OPTIONAL_FILE {datadir}/optional_restraints.txt
    NOPT 0
END_SECTION"""
        )


def test_rdock(datadir):
    datadir = str(datadir)
    gridfile = os.path.join(datadir, "rdock-grid.prm")

    _write_gridfile(datadir, gridfile)

    cmd = ["rbcavity", "-was", "-d", "-r", gridfile]

    returncode = subprocess.call(cmd)
    if returncode != 0:
        raise RuntimeError("The grid was not created. Check the parameters passed")

    for ff in ["rdock-grid.as", "rdock-grid_cav1.grd"]:
        if not os.path.exists(os.path.join(datadir, ff)):
            raise RuntimeError(f"Could not find result file {ff}")
