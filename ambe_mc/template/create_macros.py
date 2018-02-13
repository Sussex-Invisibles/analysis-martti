#!/usr/bin/python
import sys
if (len(sys.argv) != 4):
    print "ERROR: Number of arguments must be 3! directory, material, number of jobs"
    exit(1)
for num in range (1, int(sys.argv[3])+1):
    fin = open(sys.argv[1]+"/AmBe_"+sys.argv[2]+".mac","r")
    fout = open(sys.argv[1]+"/macros/AmBe_"+sys.argv[2]+"_"+str(num)+".mac","w")
    for line in fin:
        s = line.replace("{0}", str(num))
        fout.write(s)
    fout.close()
    fin.close()
