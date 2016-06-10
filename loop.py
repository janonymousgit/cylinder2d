#!/usr/bin/python3
# -*- coding: utf-8 -*-
from subprocess import call
import datetime
from getpass import getuser
import os
from multiprocessing import cpu_count
import frogress

# create output directory
d = datetime.datetime.now()
timestamp = str(d.year) + str(d.month) + str(d.day) + str(d.hour)  + str(d.minute) + str(d.second)
uid = getuser()
directory = timestamp + "_" + uid
os.mkdir(directory)

# execution parameters
angles = range(360)
angles = [ang for ang  in angles if ang == 54]
# np = str(cpu_count())
np = "3"
print("starting time: " + str(d.time())[0:8])
print("number of cores: " + np)
print("output directory: " + directory)

for ang in frogress.bar(angles):
# for ang in angles:
    #print("running simulation " + str(ang+1) + "/" + str(len(angles)))
    outputfile = directory + "/" + str(ang) + ".txt"
    command = ["mpirun", "-np", np, "cylinder2d", str(ang)]
    with open(outputfile, "w") as f:
        call(command, stdout=f)
