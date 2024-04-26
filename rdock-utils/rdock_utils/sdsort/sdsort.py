import logging

from .parser import SDSortConfig

logger = logging.getLogger("sdsort")


class SDSort:
    def __init__(self, config: SDSortConfig) -> None:
        self.descending_sort = config.descending_sort
        self.numeric_sort = config.numeric_sort
        self.fast_mode = config.fast_mode
        self.name_field = config.name_field
        self.data_field = config.data_field
        self.files = config.files

    pass



import sys
import argparse

# Mocking SDRecord class for demonstration
class SDRecord:
    def __init__(self):
        self.DATA = {}
    def readRec(self, DATA=False, LINES=False):
        pass
    def addData(self, **kwargs):
        pass
    def copy(self, DATA=False, LINES=False):
        pass
    def writeRec(self):
        pass

# Global variables
SDSORTKEY = None
SDSORTASCEND = True  # 1 = ascending, 0 = descending
SDSORTTEXT = True  # 1 = text sort, 0 = numeric sort
FASTFORMAT = False
FASTKEY = "_TITLE1"

def printHelpAndExit():
    print("\nSorts SD records by given data field\n")
    print("Usage:\tsdsort [-n] [-r] [-f <DataField>] [sdFiles]\n")
    print("\t-n\t\tnumeric sort (default is text sort)\n")
    print("\t-r\t\tdescending sort (default is ascending sort)\n")
    print("\t-f <DataField>\tspecifies sort field\n")
    print("\t-s\t\tfast mode. Sorts the records for each named compound independently (must be consecutive)\n")
    print("\t-id <NameField>\tspecifies compound name field (default = 1st title line)\n\n")
    print("Note:\t_REC (record #) is provided as a pseudo-data field\n")
    print("\n\tIf SD file list not given, reads from standard input\n")
    print("\tOutput is to standard output\n")
    print("\tFast mode can be safely used for partial sorting of huge SD files of raw docking hits\n")
    print("\twithout running into memory problems.\n\n")
    sys.exit()

def parse_args():
    parser = argparse.ArgumentParser(description='Sorts SD records by given data field')
    parser.add_argument('-r', action='store_true', help='descending sort (default is ascending sort)')
    parser.add_argument('-n', action='store_true', help='numeric sort (default is text sort)')
    parser.add_argument('-f', metavar='DataField', help='specifies sort field')
    parser.add_argument('-s', action='store_true', help='fast mode. Sorts the records for each named compound independently (must be consecutive)')
    parser.add_argument('-id', metavar='NameField', help='specifies compound name field (default = 1st title line)')
    parser.add_argument('sdFiles', nargs='*', help='SD files to sort')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    return args

args = parse_args()

SDSORTASCEND = not args.r
SDSORTTEXT = not args.n
FASTFORMAT = args.s
if args.f:
    SDSORTKEY = args.f
if args.id:
    FASTKEY = args.id

# Read records
sdRec = SDRecord()
records = []
nRec = 0
lastid = ""

def sortSD(a, b):
    if SDSORTTEXT:
        if SDSORTASCEND:
            return (a['DATA'][SDSORTKEY] > b['DATA'][SDSORTKEY]) - (a['DATA'][SDSORTKEY] < b['DATA'][SDSORTKEY])
        else:
            return (b['DATA'][SDSORTKEY] > a['DATA'][SDSORTKEY]) - (b['DATA'][SDSORTKEY] < a['DATA'][SDSORTKEY])
    else:
        if SDSORTASCEND:
            return a['DATA'][SDSORTKEY] - b['DATA'][SDSORTKEY]
        else:
            return b['DATA'][SDSORTKEY] - a['DATA'][SDSORTKEY]

while sdRec.readRec(DATA=True, LINES=True):
    nRec += 1
    sdRec.addData(_REC=nRec)  # add record# as temp data field
    if FASTFORMAT:
        id_ = sdRec.DATA[FASTKEY]
        if lastid and lastid != id_:
            for rec in sorted(records, key=lambda x: sortSD(x, x)):
                rec.writeRec()
            records = []  # clear the list
        lastid = id_
    records.append(sdRec.copy(DATA=True, LINES=True))

# Write sorted records
for rec in sorted(records, key=lambda x: sortSD(x, x)):
    rec.writeRec()
