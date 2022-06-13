#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 16:44:44 2022
"""
import numpy as np


def calc_solution(ndim, xyz, nvar):
    ans = np.zeros(nvar)
    for i in range(0,nvar):
        ans[i] = i + 3.1415;
        for j in range(0,ndim,2):
            ans[i] += 0.13 * xyz[j]
        for j in range(1,ndim,2):
            ans[i] -= 0.10 * xyz[j]
    return ans


