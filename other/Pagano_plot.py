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

def mathtext(t):
    return '\\textit{%s}' % t

import matplotlib
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'
from matplotlib import rc
font = {'family' : 'Times New Roman',
        'size'   : 16}
rc('font', **font)
rc('text', usetex=True)
import matplotlib.pyplot as plt

plt.figure(1)

inputcsv = 'Pagano_3lay_cyl_bend_conv_H8_50_default_Stress.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "x", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "DCE")

inputcsv = 'Pagano_3lay_cyl_bend_conv_MSH8_50_extraptrend_Stress.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "o", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "TBE")

inputcsv = 'Pagano_3lay_cyl_bend_conv_MSH8_50_extrapmean_Stress.CSV'
x = coldata(inputcsv, 0)
y = coldata(inputcsv, 3)
plt.loglog(x, y, color = "black", marker = "d", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "MSOE")

plt.axis([0.004, 0.1, 0.01, 0.7])

# inputcsv = 'Pagano_3lay_cyl_bend_conv_H8_4_default_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "x", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "DCE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MSH8_4_extraptrend_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "o", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "TBE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MSH8_4_extrapmean_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "d", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "MSOE")
#
# plt.axis([0.004, 0.1, 0.01, 0.7])

# inputcsv = 'Pagano_3lay_cyl_bend_conv_T10_4_default_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "x", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "DCE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MST10_4_extraptrend_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "o", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "TBE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MST10_4_extrapmean_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "d", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "MSOE")
#
# plt.axis([0.004, 0.1, 0.005, 0.4])

# inputcsv = 'Pagano_3lay_cyl_bend_conv_T10_50_default_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "x", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "DCE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MST10_50_extraptrend_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "o", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "TBE")
#
# inputcsv = 'Pagano_3lay_cyl_bend_conv_MST10_50_extrapmean_Stress.CSV'
# x = coldata(inputcsv, 0)
# y = coldata(inputcsv, 3)
# plt.loglog(x, y, color = "black", marker = "d", markersize = 12, markerfacecolor='none', markeredgecolor='k', label = "MSOE")
#
# plt.axis([0.004, 0.1, 0.005, 0.2])

plt.xlabel('Element size')
plt.ylabel('Normalized error')

plt.grid(True)
legend = plt.legend(loc='best', shadow=False, fontsize='large')
plt.show()