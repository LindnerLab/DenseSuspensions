# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:32:05 2020

This function inserts a row (row) at a given index (idx) in a matrix (matrix)

@author: Lars Kool
"""

def insertRow(matrix,row,idx):
    start = matrix[:idx-1]
    end = matrix[idx:]
    result = pd.concat(start,row,end).reset_index(drop=True)
    return result