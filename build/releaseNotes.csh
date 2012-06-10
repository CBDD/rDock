#!/bin/csh -f
# Creates automated release notes for rDock
# Requires connection to Perforce depot
#
# 4 files are created:
# CHANGES - all changelists since branch created
# JOBS - all fixed jobs since branch created
# FILES - all changed files since branch created
# VERSIONS - last changelist submitted in each main subdirectory
#
# NOTE: Release notes assume that client is sync'd to head revision
# except for VERSIONS which lists the changelists for the current
# client revision (#have rather than #head)
#
mkdir -p ../docs
set changes = ../docs/CHANGES
set jobs = ../docs/JOBS
set files = ../docs/FILES
set versions = ../docs/VERSIONS
set basever = ../...@`./p4utils.pl BASEVER`
set currentver = ../...\#head
#
# First the changes
if (-e ${changes}1) /bin/rm -f ${changes}1
if (-e ${changes}2) /bin/rm -f ${changes}2
if (-e ${changes}) /bin/rm -f ${changes}
p4 changes -l -i ${basever} > ${changes}1
p4 changes -l -i ${currentver} > ${changes}2
diff ${changes}1 ${changes}2 > ${changes}
/bin/rm -f ${changes}1
/bin/rm -f ${changes}2
#
# Now the jobs
if (-e ${jobs}1) /bin/rm -f ${jobs}1
if (-e ${jobs}2) /bin/rm -f ${jobs}2
if (-e ${jobs}) /bin/rm -f ${jobs}
p4 jobs -l -i ${basever} > ${jobs}1
p4 jobs -l -i ${currentver} > ${jobs}2
diff ${jobs}1 ${jobs}2 > ${jobs}
/bin/rm -f ${jobs}1
/bin/rm -f ${jobs}2
#
# Now the files
if (-e ${files}1) /bin/rm -f ${files}1
if (-e ${files}2) /bin/rm -f ${files}2
if (-e ${files}) /bin/rm -f ${files}
p4 files ${basever} > ${files}1
p4 files ${currentver} > ${files}2
diff ${files}1 ${files}2 > ${files}
/bin/rm -f ${files}1
/bin/rm -f ${files}2
#
# Now the versions
if (-e ${versions}) /bin/rm -f ${versions}
#p4 where ../... > ${versions}
#p4 changes -m 1 -i ${currentver} >> ${versions}
foreach dir (. bin build data docs include src)
  set havever = ../${dir}/...\#have
  echo "" >> ${versions}
  echo "LATEST CHANGELIST FOR ${havever}" >> ${versions}
  p4 where ../${dir}/... >> ${versions}
  p4 changes -m 1 -i ${havever} >> ${versions}
end
