#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:13:18 2026

@author: marlonalejandrocalderonbalcazar

calculate Anscombe residuals between background and foreground 2D SFS

"""

import sys
from collections import defaultdict

# arguments
sfs_file = sys.argv[1]
output_file = sys.argv[2]

# def calc_residual(bg_count, fg_count):
    
#     if bg_count <= 0 or fg_count <= 0:
#         return float("nan")
    
#     numerator = (bg_count**(2/3) - bg_count**(-1/3)/9) - (fg_count**(2/3) - fg_count**(-1/3)/9)
#     denominator = bg_count**(1/6)
    
#     residual = 3/2*(numerator/denominator)
    
#     return residual

# poisson 
def calc_residual(bg_count, fg_count):
    
    if bg_count <= 0 or fg_count <= 0:
        return float("nan")
    
    numerator = bg_count - fg_count
    denominator = bg_count**(1/2)
    
    residual = numerator/denominator
    
    return residual


# loop over file with sfs data
sfs = defaultdict(lambda: defaultdict(dict))

with open(sfs_file) as f:
    header = next(f)
    
    for line in f:
        region, gen, i, j, count, density = line.strip().split()
        gen = int(gen)
        i = int(i)
        j = int(j)
        count = int(count)
        density = float(density)
        sfs[gen][region][(i,j)] = density
        
with open(output_file, 'w') as out:
    out.write("generation\tfreq_p1\tfreq_p2\tresidual\n")
    
    for gen in sorted(sfs.keys()):
        bg = sfs[gen].get("bg", 0)
        fg = sfs[gen].get("fg", 0)
        
        all_cells = set(bg.keys()) | set(fg.keys())
        
        for i,j in sorted(all_cells):
            bg_count = bg.get((i,j), 0)
            fg_count = fg.get((i,j), 0)
        
            res = calc_residual(bg_count, fg_count)
            
            out.write(f"{gen}\t{i}\t{j}\t{res}\n")
        
    
