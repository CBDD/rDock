#!/bin/csh -f
# Creates vdW grids for all receptor prm files listed on command line
#
set gridstep = 0.3
set border = 1.0

foreach grid ('vdw1' 'vdw5') 
  foreach file ($argv)
      $RBT_ROOT/bin/rbcalcgrid -r${file} -pcalcgrid_${grid}.prm -o_${grid}.grd -g${gridstep} -b${border}
  end
end
