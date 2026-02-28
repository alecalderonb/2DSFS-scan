#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 12:23:25 2026

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


def calculate_1d_sfs(
    data_dict,
    pop,
    pop_size,
    fold=False
):
    """
    compute 1D SFS from a SNP dictionary
    
    pop: population ID
    pop_size: number of diploid individuals
    fold: whether to fold the SFS
    """

    n = pop_size * 2  # total chromosomes

    # initialize bins
    sfs = {i: 0 for i in range(n + 1)}

    total_sites = 0

    for snp_info in data_dict.values():

        ref, alt = snp_info["calls"].get(pop, (0, 0))

        # skip missing population
        if ref + alt == 0:
            continue

        # skip monomorphic sites
        if alt == 0 or alt == n:
            continue

        
        sfs[alt] += 1
        total_sites += 1

    if fold:
        num_chr = max(sfs.keys())
        folded_sfs = {}
        
        for freq, count in sfs.items():
            
            if count == 0:
                continue
            
            minor_freq = min(freq, num_chr - freq)
            
            if minor_freq in folded_sfs:
                folded_sfs[minor_freq] += count
            else:
                folded_sfs[minor_freq] = count
        
        return folded_sfs, total_sites

    return sfs, total_sites

def write_1d_sfs(
    sfs,
    pop_size,
    outfile,
    fold=False
):
    """
    write 1D SFS to file
    """

    total = sum(sfs.values())

    with open(outfile, "w") as f:
        f.write("freq\tcount\tdensity\n")

        for i in sorted(sfs.keys()):

            if fold:
                freq = i
            else:
                freq = i

            count = sfs[i]
            density = count / total if total > 0 else 0

            f.write(f"{freq}\t{count}\t{density}\n")
            
## foreground (chromosome Z: 4000000-4500000)
data = make_data_dict_vcf(
    vcf_filename="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.2.vcf.gz",
    popinfo_filename="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/popmap.txt",
    chrom="NC_087119.1",
    start=4000000,
    end=4500000,
    variant_type=None
)

data_chrZ = {
    snp_id: snp_info
    for snp_id, snp_info in data.items()
    if (
        snp_id.split("-")[0] == "NC_087119.1"
        and 4000000 <= int(snp_id.split("-")[1]) <= 4500000
    )
}

# uv
sfs, total_sites = calculate_1d_sfs(
    data_chrZ,
    pop="uv",
    pop_size=18,
    fold=True
)

write_1d_sfs(
    sfs,
    pop_size=18,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/chrZ_4000000-4500000.uv_sfs.txt"
)

# bv
sfs, total_sites = calculate_1d_sfs(
    data_chrZ,
    pop="bv",
    pop_size=14,
    fold=True
)

write_1d_sfs(
    sfs,
    pop_size=14,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/chrZ_4000000-4500000.bv_sfs.txt"
)

## background (autosome 1)
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

# uv
sfs, total_sites = calculate_1d_sfs(
    data,
    pop="uv",
    pop_size=18,
    fold=True
)

write_1d_sfs(
    sfs,
    pop_size=18,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/chr1.uv_sfs.txt"
)

# bv
sfs, total_sites = calculate_1d_sfs(
    data,
    pop="bv",
    pop_size=14,
    fold=True
)

write_1d_sfs(
    sfs,
    pop_size=14,
    outfile="/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/emp_data/data/chr1.bv_sfs.txt"
)


