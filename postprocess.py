"""
This contains the post-processing steps
"""
import numpy as np
import matplotlib.pyplot as plt

calpath = "output_files/mass_fits/bf_cal_ps%d.txt"

if __name__ == "__main__":
    ps = [15, 25, 35]
    cal1 = np.loadtxt(calpath%15)
    cal2 = np.loadtxt(calpath%25)
    cal3 = np.loadtxt(calpath%35)

    for i in range(len(cal1)):
        for j in range(len(cal1[i])):
            print cal1[i,j], cal2[i,j], cal3[i,j]
