#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

"""

import re
import numpy as np
import matplotlib.pyplot as plt
import os

def dragLiftPlot(drag, lift, filename, filetype):
    fig = plt.figure()
    plt.plot(list(range(len(drag))), drag)
    plt.plot(list(range(len(lift))), lift)
    plt.legend(['drag', 'lift'], loc='center right')
    filename = filename + "." + filetype
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


# ENTER DIRECTORY FROM loop.py HERE
directory = ""

# Otherwise the newest one is chosen
if directory == "":
    all_subdirs = [d for d in os.listdir('.') if os.path.isdir(d) if not d.startswith('.')]
    directory = max(all_subdirs, key=os.path.getmtime)

directory += "/"
print(directory)

# directory for plot output
imgPath = directory + "plots/"
if not os.path.exists(imgPath):
    os.mkdir(imgPath)


reg = "drag=(-*\d*.\d*e*-*\d*); lift=(-*\d*.\d*e*-*\d*)"
endDrag = []
endLift = []

for file in os.listdir(directory):
    if not os.path.isdir(directory + file):
        with open(directory + file, "r") as f:
            contents = f.read()
            hits = re.findall(reg, contents)
            endDrag.append(hits[-1][0])
            endLift.append(hits[-1][1])
            drag = [float(hit[0]) for hit in hits]
            lift = [float(hit[1]) for hit in hits]
            filename = imgPath + file[0:-4]
            dragLiftPlot(drag,lift, filename, "png")


dragLiftPlot(endDrag, endLift, imgPath + "all", "png")
