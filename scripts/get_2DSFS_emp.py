#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 14:51:35 2026

@author: marlonalejandrocalderonbalcazar
"""

import gzip
import pickle
import bz2

def make_data_dict_vcf(
    vcf_filename,
    popinfo_filename,
    chrom=None,
    start=None,
    end=None,
    variant_type=None,
    pickleWrite=False,
    outFile=''
):
    """
    Parse a VCF and return a SNP dictionary with allele counts per population.

    Optional filters:
    - chrom, start, end : restrict to genomic region
    - variant_type     : filter annotation (e.g. 'missense', 'synonymous')
    """

    # read popmap
    popmap = {}
    with open(popinfo_filename) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 2:
                popmap[cols[0]] = cols[1]

    data_dict = {}
    poplist = []

    with gzip.open(vcf_filename, "rt") as vcf:
        for line in vcf:

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                header = line.strip().split()
                for sample in header[9:]:
                    poplist.append(popmap.get(sample))
                continue

            cols = line.strip().split("\t")

            chr_id = cols[0]
            pos = int(cols[1])

            # region filter
            if chrom is not None and chr_id != chrom:
                continue
            if start is not None and pos < start:
                continue
            if end is not None and pos > end:
                continue

            # filter status
            if cols[6] not in ("PASS", "."):
                continue

            ref = cols[3].upper()
            alt = cols[4].upper()
            if ref not in "ACGT" or alt not in "ACGT":
                continue

            # annotation
            info_parts = cols[7].split("|")
            annotation = info_parts[1] if len(info_parts) >= 2 else None
            if variant_type is not None and annotation != variant_type:
                continue

            snp_id = f"{chr_id}-{pos}"
            snp_dict = {
                "segregating": (ref, alt),
                "context": f"-{ref}-",
                "annotation": annotation,
                "calls": {}
            }

            gt_index = cols[8].split(":").index("GT")

            for pop, sample in zip(poplist, cols[9:]):
                if pop is None:
                    continue

                gt = sample.split(":")[gt_index]
                refc, altc = snp_dict["calls"].get(pop, (0, 0))
                refc += gt[::2].count("0")
                altc += gt[::2].count("1")
                snp_dict["calls"][pop] = (refc, altc)

            data_dict[snp_id] = snp_dict

    if pickleWrite:
        with bz2.open(outFile, "wb") as f:
            pickle.dump(data_dict, f)

    return data_dict


def calculate_2d_sfs(
    data_dict,
    pop1,
    pop2,
    pop1_size,
    pop2_size,
    fold=False
):
    """
    Compute a dadi-compatible 2D SFS.
    """

    n1 = pop1_size * 2
    n2 = pop2_size * 2
    max_total = (n1 + n2) // 2
    
    sfs = {}
    
    # NOTE: folded bins only go up to N, not 2N
    if fold:
        for i in range(pop1_size + 1):
            for j in range(pop2_size + 1):
                if i + j <= min(pop1_size, pop2_size):
                    sfs[(i, j)] = 0
    else:
        # full rectangular unfolded SFS
        for i in range(n1 + 1):
            for j in range(n2 + 1):
                sfs[(i, j)] = 0

    total_sites = 0

    for snp_info in data_dict.values():

        # allele counts
        ref1, alt1 = snp_info['calls'].get(pop1, (0, 0))
        ref2, alt2 = snp_info['calls'].get(pop2, (0, 0))

        # skip monomorphic
        if alt1 == 0 and alt2 == 0:
            continue

        i, j = alt1, alt2

        if fold:
            # global folding (already correct)
            total_alt = i + j
            if total_alt > max_total:
                i = n1 - i
                j = n2 - j
        
            # convert allele counts → individuals
            i = min(i, n1 - i)
            j = min(j, n2 - j)
        
            # enforce triangular shape
            if i + j > min(pop1_size, pop2_size):
                continue
            
            # now guaranteed:
            # 0 ≤ i ≤ pop1_size
            # 0 ≤ j ≤ pop2_size

        sfs[(i, j)] += 1
        total_sites += 1

    return sfs, total_sites

def write_2d_sfs(
    sfs,
    pop1_size,
    pop2_size,
    outfile,
    fold=False
):
    """
    Write 2D SFS to file with:
    freq_p1 = number of individuals carrying the allele in pop1
    freq_p2 = number of individuals carrying the allele in pop2
    """

    total = sum(sfs.values())

    with open(outfile, "w") as f:
        f.write("freq_p1\tfreq_p2\tcount\tdensity\n")

        for (i, j), count in sorted(sfs.items()):

            if fold:
                # convert allele counts → individuals
                freq_p1 = min(i, 2 * pop1_size - i)
                freq_p2 = min(j, 2 * pop2_size - j)
            else:
                freq_p1 = i
                freq_p2 = j

            density = count / total if total > 0 else 0
            f.write(f"{freq_p1}\t{freq_p2}\t{count}\t{density}\n")
            
            
# load snp data from vcf
data = make_data_dict_vcf(
    vcf_filename="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/data_summer2024/ECBAnnotated.vcf.gz",
    popinfo_filename="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/popmap.txt",
    chrom=None,
    start=None,
    end=None,
    variant_type=None,
)

# load snp dictionary from pickle file
with bz2.open("/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/data_summer2024/genome_data.pkl.bz2", 'rb') as bz_file:
    data = pickle.load(bz_file)

# background (autosome 1)
data_chr1 = {
    snp_id: snp_info
    for snp_id, snp_info in data.items()
    if (
        snp_id.split("-")[0] == "NC_087088.1"
    )
}

sfs, total_sites = calculate_2d_sfs(
    data_chr1,
    pop1="uv",
    pop2="bv",
    pop1_size=18,
    pop2_size=14,
    fold=True
)

write_2d_sfs(
    sfs,
    pop1_size=18,
    pop2_size=14,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/chr1.sfs.txt"
)    

# foreground (chromosome Z: 4000000-4500000)
data_chrZ = {
    snp_id: snp_info
    for snp_id, snp_info in data.items()
    if (
        snp_id.split("-")[0] == "NC_087119.1"
        and 4000000 <= int(snp_id.split("-")[1]) <= 4500000
    )
}

sfs, total_sites = calculate_2d_sfs(
    data_chrZ,
    pop1="uv",
    pop2="bv",
    pop1_size=18,
    pop2_size=14,
    fold=True
)

write_2d_sfs(
    sfs,
    pop1_size=18,
    pop2_size=14,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/chrZ_4000000-4500000.sfs.txt"
)


    
    
