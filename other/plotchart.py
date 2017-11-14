# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import csv


def loadcsv(inputcsv):
    rows = []
    with open(inputcsv, newline='') as csvfile:
        sreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        headerrow = next(sreader)
        for row in sreader:
            if row:
                rows.append(row)
    return headerrow, rows

def coldata(inputcsv, theset):
    headerrow, rows = loadcsv(inputcsv)
    return [r[theset] for r in rows]

import matplotlib.pyplot as plt
plt.figure(1)

inputcsv = 'Pagano_3lay_cyl_bend_THIN_T10_conv_Stress_errors_default.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "x", markersize = 12, label = "DCE")

inputcsv = 'Pagano_3lay_cyl_bend_THIN_MST10_conv_Stress_errors_extraptrend.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "o", markersize = 12, label = "TBE")

inputcsv = 'Pagano_3lay_cyl_bend_THIN_MST10_conv_Stress_errors_extrapmean.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "d", markersize = 12, label = "MSOE")

plt.xlabel('Element size')
plt.ylabel('Normalized error')
plt.axis([0.004, 0.1, 0.0001, 1.0])
plt.grid(True)
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.show()