import sys

ref=[item2.split() for item2 in [item.rstrip() for item in open(sys.argv[1]).readlines()]]
test=[item2.split() for item2 in [item.rstrip() for item in open(sys.argv[2]).readlines()]]

refcoords=ref[4:48]
testcoords=test[4:48]

refscore=float(ref[123][0])
testscore=float(test[123][0])

checkcoords=[]
for i in range(0,len(refcoords)):
	for j in [0,1,2]:
		checkcoords.append(abs(float(refcoords[i][j])-float(testcoords[i][j])))
error=0
for check in checkcoords:
	if check > 0.01:
		error=1
if abs(refscore - testscore) > 0.01:
	error = 1

if error == 1:
	print "The test failed, please check the compilation is OK and no errors were raised."
else:
	print "The test succeeded! The results agree with the reference ones.\nHave fun using rDock!!"
