# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 11:26:42 2020

@author: Charly
"""

def neighboringElements(matrix):
    matrix = list(matrix)
    result = []
    temp = [matrix[0]]
    last = matrix[0]
    for i in range(1,len(matrix)):
        if  matrix[i] - 1 == last:
            temp.append(matrix[i])
        else:
            result.append(temp)
            temp = [matrix[i]]
        last = matrix[i]
        
    result.append(temp)
    return result