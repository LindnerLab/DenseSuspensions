# -*- coding: utf-8 -*-
"""
Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France
               
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
"""

def neighboring_elements(matrix):
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