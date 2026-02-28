#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 11:33:09 2026

@author: marlonalejandrocalderonbalcazar

calculate Anscombe residuals between background and foreground 2D SFS in empirical data

"""

import sys
from collections import defaultdict

# arguments
sfs_bg = sys.argv[1]
sfs_fg = sys.argv[2]
output_file = sys.argv[3]

# anscombe
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
sfs_bg_dict = defaultdict(lambda: defaultdict(dict))

with open(sfs_bg) as f:
    header = next(f)
    
    for line in f:
        i, j, count, density = line.strip().split()
        i = int(i)
        j = int(j)
        count = int(count)
        density = float(density)
        sfs_bg_dict[(i,j)] = density
        
        
sfs_fg_dict = defaultdict(lambda: defaultdict(dict))

with open(sfs_fg) as f:
    header = next(f)
    
    for line in f:
        i, j, count, density = line.strip().split()
        i = int(i)
        j = int(j)
        count = int(count)
        density = float(density)
        sfs_fg_dict[(i,j)] = density


all_cells = set(sfs_bg_dict.keys()) | set(sfs_fg_dict.keys())

with open(output_file, "w") as out:
    out.write("freq_p1\tfreq_p2\tresidual\n")

    for i, j in sorted(all_cells):
        bg_count = sfs_bg_dict.get((i, j), 0)
        fg_count = sfs_fg_dict.get((i, j), 0)

        res = calc_residual(bg_count, fg_count)

        out.write(f"{i}\t{j}\t{res}\n")