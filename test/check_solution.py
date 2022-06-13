#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 16:44:44 2022
"""

import sys
import numpy as np


def calc_solution(xyz, nvar):
    ans = np.zeros(nvar)
    for i in range(0,nvar):
        ans[i] = i + 3.1415;
        for j in range(0,len(xyz),2):
            ans[i] += 0.13 * xyz[j]
        for j in range(1,len(xyz),2):
            ans[i] -= 0.10 * xyz[j]
    return ans

# Get Solution Filename
answer_file_name = "answer.txt"
if len(sys.argv) > 1:
    answer_file_name = sys.argv[1]

# Read Answers
header = np.loadtxt(answer_file_name, dtype=int, max_rows=1)
data = np.loadtxt(answer_file_name, dtype=float, skiprows=1)
ndim = header[0]
npnt = header[1]
nvar = header[2]

# Check Answers
max_error = 0.0
for row in range(0,npnt):
    ans = calc_solution( data[row][0:ndim], nvar)
    max_error = np.max( [max_error, np.max( (data[row][ndim:]-ans) / ans )] )

# Display % Error From Exact Solution
print('Maximum %% Error = %12.3e'%(max_error))
