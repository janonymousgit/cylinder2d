#!/usr/bin/python3
# -*- coding: utf-8 -*-
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



directory = "2016611153745_jan/"
# file = "2016610221350_jan/0.txt"
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
            filename = directory + "plots/" + file[0:-4]
            dragLiftPlot(drag,lift, filename, "png")

dragLiftPlot(endDrag, endLift, directory + "plots/" + "all", "png")
