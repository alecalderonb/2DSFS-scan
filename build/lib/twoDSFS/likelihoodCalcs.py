import gzip
import re
import math
from scipy.stats import poisson
from scipy.stats import multinomial
import numpy as np
from scipy.optimize import minimize
import csv
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import bz2
import seaborn as sns
import matplotlib.colors as mcolors

 
class LikelihoodInference_jointSFS():
    
    def __init__(self, vcf_filename='', popinfo_filename='data/popmap.txt', start_position=None, end_position=None, 
                 pop1='uv', pop2='bv', pop1_size=18, pop2_size=14, variant_type=None, fold=True, regen=True, sim_back = []):
        
        self.vcf_filename = vcf_filename 
        self.popinfo_filename = popinfo_filename
        self.pop1 = pop1
        self.pop2 = pop2
        self.pop1_size = pop1_size
        self.pop2_size = pop2_size
        self.start_position = start_position
        self.end_position = end_position
        self.variant_type = variant_type
        self.fold = fold
        self.bg_sfs = {}
        self.regen = regen
        self.sim_back = sim_back
 
    def make_composite_VCF(self, migRate, gen, noDel=False):
 
        # given a directory with simulated VCF files, combine them into one VCF with all the data to run the analysis
 
        delString = ""
        if noDel:
            delString = "_noDel"   
 
        migDir = "slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams"+delString+"/migRate."+str(migRate)+"/"
        outFile = "slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams"+delString+"/migRate."+str(migRate)+"."+str(gen)+".vcf.gz"

        print_header = True
        out =  gzip.open(outFile, 'wt')

        # loop through each file in migDir
        for item_name in os.listdir(migDir):
            full_path = os.path.join(migDir, item_name)  
            for file_name in os.listdir(full_path):
                if gen in file_name and "gz" in file_name: 
                    file_path = os.path.join(full_path, file_name)
                    
                    fh = gzip.open(file_path, 'rt')
                    for line in fh:
                        if line[0] == "#":
                            if print_header:
                                out.write(line.strip()+"\n")
                        else:
                            out.write(line.strip()+"\n")
                    print_header = False    
                    fh.close()

        out.close()    

    def make_data_dict_vcf(self, vcf_filename, popinfo_filename, pickleWrite = False, outFile = ''):
        """
        parse a vcf file and return a dictionary containing 'calls', 'context', 
        and 'segregating' keys for each SNP. 
        - 
        - calls: dictionary where the keys are the population ids, and the values
                 are tuples with two entries: number of individuals with the reference 
                 allele, and number of individuals with the derived allele. 
        - context: reference allele.
        - segregating: tuple with two entries: reference allele and derived allele.

        arguments:
        - vcf_filename: name of the vcf file containing SNP data.
        - popinfo_filename: file containing population info for each sample.
        - filter: if True, only include SNPs that passed the filter criteria.
        
        add another argument that works as a flag to take only whatever category of function I want (missense or synonymous)

        """
        
        # make population map dictionary
        popmap_file = open(popinfo_filename, "r")
        popmap = {}

        for line in popmap_file:
            columns = line.strip().split()
            if len(columns) >= 2:
                popmap[columns[0]] = columns[1]
        popmap_file.close()

        #open files
        vcf_file = gzip.open(vcf_filename, 'rt')
        
        data_dict = {} # initialize output dictionary
        
        poplist = []
        
        for line in vcf_file:

            if line.startswith('##'): # skip vcf metainfo
                continue
            if line.startswith('#'): # read header
                header_cols = line.split()
                
                # map sample IDs to popmap file
                for sample in header_cols[9:]:
                    if sample in popmap:
                        poplist.append(popmap[sample])
                    else:
                        None
                continue

            cols = line.split("\t")
            snp_id = '-'.join(cols[:2]) # keys for data_dict: CHR_POS
 
            if snp_id in data_dict:
                start = '.01'
                while snp_id in data_dict:
                    snp_id = '-'.join(cols[:2])+start
                    start = str(float(start)+0.01)[1:]

            snp_dict = {} 
            
            # filter SNPs based on 'snp_type' annotation
            info_field = cols[7]
            info_field_parts = info_field.split('|')
            if len(info_field_parts) >= 2:
                annotation = info_field_parts[1]
            else: 
                annotation = 'No annotation'

            if cols[6] != 'PASS' and cols[6] != '.':
                continue
            
            # make alleles uppercase
            ref = cols[3].upper()
            alt = cols[4].upper()
            
            if ref not in ['A', 'C', 'G', 'T'] or alt not in ['A', 'C', 'G', 'T']:
                continue

            snp_dict['segregating'] = (ref, alt)
            snp_dict['context'] = '-' + ref + '-'

            calls_dict = {}
            gtindex = cols[8].split(':').index('GT') # extract the index of the GT field

            # pair each pop from poplist with the corresponding sample genotype data from cols[9:]
            for pop, sample in zip(poplist, cols[9:]):
                if pop is None:
                    continue
                
                gt = sample.split(':')[gtindex] # extract genotype info
                if pop not in calls_dict:
                    calls_dict[pop] = (0, 0)
                
                # count ref and alt alleles
                refcalls, altcalls = calls_dict[pop]
                refcalls += gt[::2].count('0') # gt[::2] slices the genotype and skips the '/' or '|' character
                altcalls += gt[::2].count('1')
                calls_dict[pop] = (refcalls, altcalls)

            snp_dict['calls'] = calls_dict
            snp_dict['annotation'] = annotation
            data_dict[snp_id] = snp_dict

        vcf_file.close()
        
        if pickleWrite:
            with bz2.open(outFile, "wb") as file:
                pickle.dump(data_dict, file)

        return data_dict


    def calculate_2d_sfs(self, data_dict):
        """
        calculate the two-dimensional sfs 
        for two populations from a given SNP data dictionary.

        parameters:
        - data_dict: dictionary containing SNP information. Each entry includes allele counts for populations.
        - pop1: name of population 1 in the dict
        - pop2: name of population 2 in the dict
        
        returns:
        - sfs_dict: dictionary where keys are tuples (p1_freq, p2_freq) and values are counts of SNPs with those frequencies.
        
        add a 1 to the bins where I have zero counts
        """
        num_genomes_p1 = self.pop1_size*2
        num_genomes_p2 = self.pop2_size*2
        sfs_dict = {}
        
        for i in range(num_genomes_p1 + 1):
            for j in range(num_genomes_p2 + 1):
                sfs_dict[(i,j)] = 0 
        
        total_sites = 0
        
        # loop through all snps in the data_dict
        for snp_id, snp_info  in data_dict.items():
            
            chr_id, pos = snp_id.split('-')
            pos = float(pos)
            
            # filter by variant type if specified
            snp_annotation = snp_info.get('annotation')
            if self.variant_type is not None and snp_annotation != self.variant_type:
                continue
            
            # get allele counts for pop1 and pop2
            pop1_calls = snp_info['calls'].get(self.pop1, (0, 0))  # (ref_calls, alt_calls)
            pop2_calls = snp_info['calls'].get(self.pop2, (0, 0))

            pop1_calls_list = list(pop1_calls)
            pop2_calls_list = list(pop2_calls)
            
            # fold based on the identity of the MAF for each pop
            
            if self.fold:
                if pop1_calls_list[1] + pop2_calls_list[1] > self.pop1_size + self.pop2_size: 
                       oldAlt1 = pop1_calls_list[1]      
                       pop1_calls_list[1] = pop1_calls_list[0]
                       pop1_calls_list[0] = oldAlt1
                       
                       oldAlt2 = pop2_calls_list[1]      
                       pop2_calls_list[1] = pop2_calls_list[0]
                       pop2_calls_list[0] = oldAlt2
            
            alt_count_pop1 = pop1_calls_list[1]
            alt_count_pop2 = pop2_calls_list[1]

            # skip snps where both pops are missing or have no alternate alleles
            if alt_count_pop1 == 0 and alt_count_pop2 == 0:
                continue

            # filter for rare alleles  
            #if alt_count_pop1 + alt_count_pop2 <= 1 or alt_count_pop1 + alt_count_pop2  >=  2*self.pop1_size+2*self.pop2_size -1:
            #    continue

            # increase the corresponding bin in the sfs
            sfs_dict.setdefault((alt_count_pop1, alt_count_pop2), 0) 
            sfs_dict[(alt_count_pop1, alt_count_pop2)] += 1 
            
            # increase total number of sites
            total_sites += 1
            
        return sfs_dict    
    
    def normalize_2d_sfs(self, sfs):
        
        # sum all the sites
        counts = list(self.sfs.values())
        total = sum(counts[1:-1]) # exclude first and last bin 
        # print(total)
        
        # divide each bin value by the total number of sites
        normalized_sfs = {}
        for coords, values in self.sfs.items():
            normalized_sfs[coords] = values / total
        
        return normalized_sfs
    
    def count_snps(self, window_data, variant_type):
        
        snp_count = 0
        for snp_data in window_data.values():
            if variant_type is None:
                snp_count += 1
            elif snp_data.get("annotation") == variant_type:
                snp_count += 1
        return snp_count
    
    def calculate_1d_sfs(self, data_dict, pop, pop_size, start_position, end_position, variant_type):

        num_genomes = pop_size*2
        sfs_dict = {}
        
        for i in range(num_genomes + 1):
            sfs_dict[i] = 0

        total_sites = 0
        
        for snp_id, snp_info in data_dict.items():
            chr_id, pos = snp_id.split('-')
            pos = float(pos)
    
            if start_position is not None and pos < start_position:
                continue
            if end_position is not None and pos > end_position:
                continue
    
            snp_annotation = snp_info.get('annotation')
            if variant_type is not None and snp_annotation != variant_type:
                continue
    
            pop_calls = snp_info['calls'].get(pop, (0, 0))  # (ref_calls, alt_calls)
            alt_count = pop_calls[1]
    
            if alt_count == 0:
                continue
    
            sfs_dict[alt_count] += 1
    
            total_sites += 1
    
        return sfs_dict        
    
    def fold_1d_sfs(self, sfs_dict):

        # get total number of chromosomes (2N)
        num_chromosomes = max(sfs_dict.keys())
    
        folded_sfs_dict = {}
    
        for freq, count in sfs_dict.items():
            # calculate maf
            minor_freq = min(freq, num_chromosomes - freq)
    
            # add the count to the corresponding bin in the folded sfs
            if minor_freq in folded_sfs_dict:
                folded_sfs_dict[minor_freq] += count
            else:
                folded_sfs_dict[minor_freq] = count
    
        return folded_sfs_dict
    
    def normalize_1d_sfs(self, sfs):
        self.sfs = sfs
        
        counts = list(sfs.values())
        total = sum(counts[1:-1]) # exclude first and last bin 
        # print(total)
        normalized_sfs = {}
        
        for freq, values in sfs.items():
            normalized_sfs[freq] = values / total
            
        return normalized_sfs
    
    def calculate_likelihood_1D(self, foreground_sfs, background_sfs): 
        
        bins = sorted(foreground_sfs.keys())
        
        ''' foreground '''
        # get observed counts from foreground
        counts_fg = []
        for k in bins[1:-1]:
            count = foreground_sfs[k]
            counts_fg.append(int(count))
        # print(counts_fg)
        
        # Calculate total foreground counts
        total_fg = sum(counts_fg)
    
        # Skip if total_fg is zero
        if total_fg == 0:
            # print("Warning: Foreground SFS has zero counts. Skipping calculation.")
            return None
        
        # normalize foreground counts
        probabilities_fg_norm = []
        for count in counts_fg:
            p_norm = count/total_fg
            probabilities_fg_norm.append(p_norm)
        # print(probabilities_fg)
        
        ''' background '''
        # get observed counts from background
        counts_bg = []
        for k in bins[1:-1]:
            count_bg = background_sfs[k]
            counts_bg.append(count_bg)
        # print(probabilities_bg)
        
        # # Calculate total background counts
        total_bg = sum(counts_bg)
        
        # Skip if total_bg is zero
        if total_bg == 0:
            # print("Warning: Background SFS has zero counts. Skipping calculation.")
            return None
        
        total_bg = sum(counts_bg)
        probabilities_bg_norm = []
        for count in counts_bg:
            p_norm = count/total_bg
            probabilities_bg_norm.append(p_norm)
        # print(probabilities_bg)
    
        if total_bg and total_fg is not None:
                
            #fg_counts =  multinomial.rvs(n = total_fg, p = probabilities_fg_norm) 

            #total_fg_denom = total_fg+len(fg_counts)

            #probabilities_fg_norm = []
            #for count in fg_counts:
            #    p_norm = (count+1)/total_fg_denom
            #    probabilities_fg_norm.append(p_norm)
        
            #print(probabilities_fg_norm)
   
            log_likelihood_bg = multinomial.logpmf(x=counts_fg, n=total_fg, p=probabilities_bg_norm)
            log_likelihood_fg = multinomial.logpmf(x=counts_fg, n=total_fg, p=probabilities_fg_norm)
        
            clr = 2*(log_likelihood_fg - log_likelihood_bg) 
        
        return clr 
    
    def calculate_likelihood_2D(self, foreground_2d_sfs, background_2d_sfs):
 
        if self.regen:
            background_2d_sfs  = self.regenBackground(foreground_2d_sfs,background_2d_sfs)
       
        bins = sorted(foreground_2d_sfs.keys())

        ''' foreground '''
        
        # get observed counts from foreground
        counts_fg = []
        for k in bins[1:-1]:
            count = foreground_2d_sfs[k]
            counts_fg.append(int(count))
        # print(counts_fg)
        
        # normalize foreground counts
        total_fg = sum(counts_fg)
    
        # Skip if total_fg is zero
        if total_fg == 0:
            # print("Warning: Foreground SFS has zero counts. Skipping calculation.")
            return None
        
        
        probabilities_fg_norm = []
        for count in counts_fg:
            p_norm = count/total_fg
            probabilities_fg_norm.append(p_norm)
        # print(probabilities_fg)
        
        ''' background '''
        # get observed counts from background
        counts_bg = []
        for k in bins[1:-1]:
            count_bg = background_2d_sfs[k] 
            counts_bg.append(count_bg)
        # print(probabilities_bg)
        
        # normalize background counts
        total_bg = sum(counts_bg)
        
        # Skip if total_bg is zero
        if total_bg == 0:
            # print("Warning: Background SFS has zero counts. Skipping calculation.")
            return None
        
        probabilities_bg_norm = []
        for count in counts_bg:
            p_norm = (count+1./(len(counts_bg)))/(total_bg+1)
            #p_norm = count/total_bg
            probabilities_bg_norm.append(p_norm)
        # print(probabilities_bg)
       
        if total_bg and total_fg is not None:

            #fg_counts =  multinomial.rvs(n = total_fg, p = probabilities_fg_norm)
  
            #total_fg_denom = total_fg+len(fg_counts)       

            #probabilities_fg_norm = []
            #for count in fg_counts:
            #    p_norm = (count+1)/total_fg_denom
            #    probabilities_fg_norm.append(p_norm)

            #print(probabilities_fg_norm)
            log_likelihood_bg = multinomial.logpmf(x=counts_fg, n=total_fg, p=probabilities_bg_norm)
            log_likelihood_fg = multinomial.logpmf(x=counts_fg, n=total_fg, p=probabilities_fg_norm)
        
            clr = 2*(log_likelihood_fg - log_likelihood_bg) 
        
        return clr           
    
    def T2D_scan(self, data_dict, numSnps = 100, chrom = 1):
        
        '''
        genome scan to calculate T2D
        
        arguments: 
        - data_dict:
        '''

        # sort snps by chromosome and position
        sorted_snps = []
        for snp_key in data_dict.keys():
            coords = snp_key.split('-')
            chromosome_id = coords[0]
            position = float(coords[1])
            sorted_snps.append((chromosome_id, position, snp_key))
            
        sorted_snps.sort(key=lambda x: (x[0], x[1]))
        
        chromosome_snps = {}
        # make the background SNPs, only inlcude snps on the same chromosome 
        for snp_key, snp_data in data_dict.items():
            #if snp_key.startswith(f"{chrom}-"):
            chromosome_snps[snp_key] = snp_data
         
        # limit the bg to a specific subset of positions for the simulations        
        bg_snps = {}
        if self.sim_back:
            for snp_key in chromosome_snps:
                chr_pos = snp_key.split('-')
                position = float(chr_pos[1])                    
                if position > self.sim_back[0] and position < self.sim_back[1]:
                    bg_snps[snp_key] = chromosome_snps[snp_key] 
        else:
            bg_snps = chromosome_snps     

        background_2d_sfs = self.calculate_2d_sfs(bg_snps)

        T2D_windows = {}
        current_chromosome = None
        current_window_start = 0
        window_data = {}
        pos = 0

        k = 0
        for snp in sorted_snps:
            snp_key = snp[2]
            chr_pos = snp_key.split('-')
            pos = chr_pos[1]
            if not current_chromosome:
                current_chromosome = chr_pos[0]
                current_window_start = chr_pos[1]

            # check if the SNP is on the current chromosome
            if chr_pos[0] != current_chromosome:
                foreground_2d_sfs = self.calculate_2d_sfs(window_data)
                #print(foreground_2d_sfs)
                T2D = self.calculate_likelihood_2D(foreground_2d_sfs, background_2d_sfs)
                snp_count = self.count_snps(window_data, self.variant_type)
                window_range = f"{current_chromosome} {current_window_start} {pos}"
                T2D_windows[window_range] = {
                    "snp_count": snp_count,
                    "T2D": T2D
                }
                k = 0
                current_chromosome = chr_pos[0]
                current_window_start = chr_pos[1]
                window_data = {}
                continue

            # check if SNP is within the current window
            if k < numSnps:
                window_data[snp_key] = data_dict[snp_key]  # add SNP to current window
                k += 1
            else:
                # calculate p-values for the current window
                foreground_2d_sfs = self.calculate_2d_sfs(window_data)
                T2D = self.calculate_likelihood_2D(foreground_2d_sfs, background_2d_sfs)
                snp_count = self.count_snps(window_data, self.variant_type)
                window_range = f"{current_chromosome} {current_window_start} {pos}"
                T2D_windows[window_range] = {
                    "snp_count": snp_count,
                    "T2D": T2D
                }
                    
                k = 0
                current_chromosome = chr_pos[0]
                current_window_start = chr_pos[1]
                window_data = {}

        #calculate for the last window if it has any data
        if window_data:
             foreground_2d_sfs = self.calculate_2d_sfs(window_data)
             T2D = self.calculate_likelihood_2D(foreground_2d_sfs, background_2d_sfs)
             snp_count = self.count_snps(window_data, self.variant_type)
             window_range = f"{current_chromosome} {current_window_start} {pos}"
             T2D_windows[window_range] = {
                 "snp_count": snp_count,
                 "T2D": T2D
             }
            
        return T2D_windows  

    def regenBackground(self, fg, bg):

        univariate_fg = {}
        univariate_fg_keys = {}
        tot_fg = 0
        # first, get univariate of foreground
        for mybin in fg:
            tot = mybin[0] + mybin[1]
            if tot not in univariate_fg:
                univariate_fg[tot] = 0
                univariate_fg_keys[tot] = []
            univariate_fg[tot] += fg[mybin]
            univariate_fg_keys[tot].append(mybin)
            tot_fg += fg[mybin]
 
        if tot_fg == 0:
            return bg

        bg_1d = {}
        tot_bg = 0
        # next, take 2D of bg and make a dict that stores all the bins of the 2D that correspond to each bin in the 1D fg
        for mybin in bg:
            if mybin[0]+mybin[1] not in bg_1d:
                bg_1d[mybin[0]+mybin[1]]  = []
            bg_1d[mybin[0]+mybin[1]].append(bg[mybin])
            tot_bg += bg[mybin] 

        bg_regen = {}
        # for however many sites are in FG, select sites in BG with matched joint univariate allele frequency but in proportion to their likelihood in the BG
        for freq in univariate_fg:
            for mybin in univariate_fg_keys[freq]:
                bg_regen[mybin] = 0.
            tot = sum(bg_1d[freq])
            if tot == 0:
                continue
            
            choices = np.divide(bg_1d[freq], tot)

            for k in range(0, len(choices)):
                bg_regen[univariate_fg_keys[freq][k]] = choices[k]*univariate_fg[freq]

            #for i in range(0,univariate_fg[freq]):
                  #k = np.random.choice(len(choices), p=choices)         
                  #bg_regen[univariate_fg_keys[freq][k]] += 1

        #for mybin in bg_regen:
        #    bg_regen[mybin] *= tot_bg/tot_fg

        # return new bg
        return bg_regen

