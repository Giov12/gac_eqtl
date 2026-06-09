#!/bin/env python3

import gzip
import argparse
import os

def get_arguments():
    """function to return the paths to the input files and max distance between a coordinate & a SNP"""

    parser = argparse.ArgumentParser(description="Filter VCF file for SNPs by linkage markers")
    parser.add_argument("-v", "--vcf", required=True, type=str, help="VCF containing SNPs from individuals in the linkage map")
    parser.add_argument("-l", "--linkage_map", required=True, type=str, help="Linkage map generated from LepMap3")
    parser.add_argument("-m", "--marker_coordinates", required=True, type=str, help = "Marker coordates for markers in linkage map")
    parser.add_argument("-s", "--samples_list", required=True, type=str, help = "list of samples to retain from vcf [1 per line]")
    parser.add_argument("-d", "--distance", required=False, type=int, default=100, help = "Max distance a SNP can be to a linkage marker to retain the SNP")

    args = parser.parse_args()
    vcf = args.vcf
    linkage_map = args.linkage_map
    marker_coords = args.marker_coordinates
    max_dist = args.distance
    samples_file = args.samples_list

    # make sure input files exist
    for f in [vcf, linkage_map, marker_coords, samples_file]:
        assert os.path.isfile(f), f"Could not location {f}"

    # make sure max distance is legal
    assert max_dist >= 0, "Max distance cannot be less than 0"

    return vcf, linkage_map, marker_coords, max_dist, samples_file

def parse_samples_file(samples_file: str):
    """assuming a single column-ed file with sample names to retain"""

    fh = gzip.open(samples_file, "rt") if samples_file.endswith(".gz") else open(samples_file, 'r')

    samples = set()

    for line in fh:
        samples.add(line.strip('\n'))

    fh.close()

    return samples

def make_chrom_name_map():
    """hard code the names of the chromosomes to match between linkage map & other files"""

    chrom_map = {"1": "chrI", "2": "chrII", "3": "chrIII", "4": "chrIV", "5": "chrV",
                 "6": "chrVI", "7": "chrVII", "8": "chrVIII", "9": "chrIX", "10": "chrX",
                 "11": "chrXI", "12": "chrXII", "13": "chrXIII", "14": "chrXIV", "15": "chrXV",
                 "16": "chrXVI", "17": "chrXVII", "18": "chrXVIII", "19": "chrXIX", "20": "chrXX",
                 "21": "chrXXI", "X": "chrM", "U": "chrUn", "Y": "chrY"}
    
    return chrom_map


def get_linkage_markers(linkage_map: str):
    """use only the first two lines to obtain the marker IDS and 
        assert verify that they belong in the same chromosome"""

    # open file & get the two lines
    fh = gzip.open(linkage_map, "rt") if linkage_map.endswith(".gz") else open(linkage_map, 'r')
    id_line = fh.readline()
    chrom_line = fh.readline()
    fh.close()

    # match the IDS and chroms
    id_line = id_line.replace("ID,", '').strip('\n')
    ids = id_line.split(',')
    chrom_line = chrom_line[1:].strip('\n')
    chroms = chrom_line.split(',')
    assert len(ids) == len(chroms), "Error, the number of chromosomes does not match the number of IDS in the linkage map"
    chrom_map = make_chrom_name_map()
    ids_map = {}

    for i in range(len(ids)):
        ids_map[ids[i]] = chrom_map[chroms[i]]

    return ids_map

def parse_marker_coordinates(marker_coords: str, ids_map: dict):
    """get the locations of the linkage markers"""

    # obj to return & its structure chrom : [Marker_Location]
    marker_map = {chrom : [] for chrom in ids_map.values()}

     # open file
    fh = gzip.open(marker_coords, "rt") if marker_coords.endswith(".gz") else open(marker_coords, 'r')

    # now to get the markers
    for line in fh:
        if line.startswith("locus,"):
            continue
        fields = line.strip('\n').split(',')
        # verify matching IDs & chroms
        if fields[0] in ids_map and fields[1] == ids_map[fields[0]]:
            marker_map[fields[1]].append(int(fields[2]))

    fh.close()

    return marker_map

def filter_vcf(vcf: str, marker_map: dict, max_dist: int, samples: set):
    """filter for SNPs no further than max_dist from a linkage marker"""

    msg = "WARNING: Allele frequencies will be wrong since they are not being re-calculated"
    print(msg)
    
    # create output name
    replacement = ".vcf.gz" if vcf.endswith(".gz") else ".vcf"
    outname = "Filtered_" + os.path.basename(vcf).replace(replacement, '') + f"_{max_dist}bp_Dist_SNPs.vcf"

    outfh = open(outname, 'w')

    # now for the filtering
    fh = gzip.open(vcf, "rt") if vcf.endswith(".gz") else open(vcf, 'r')

    for line in fh:
        if line[0] == '#':
            # get matching samples by index
            if line.startswith("#CHROM\t"):
                indices_2_keep = []
                fields = line.strip('\n').split('\t')
                sample_columns = []
                for i, sample in enumerate(fields):
                    if sample in samples:
                        # used to orderly traverse through the line
                        indices_2_keep.append(i)
                        sample_columns.append(sample)
                outline = '\t'.join(fields[:9]) + '\t' 
                outline = outline + '\t'.join(sample_columns) + '\n'
                outfh.write(outline)
            else:
                outfh.write(line)
        else:
            fields = line.strip('\n').split('\t')
            chrom = fields[0]
            if chrom in marker_map:
                snp_pos = int(fields[1])
                for mark_loc in marker_map[chrom]:
                    leftbound = mark_loc - max_dist
                    rightbound = mark_loc + max_dist
                    if leftbound <= snp_pos <= rightbound:
                        num_samples = 0
                        # now add only the samples to be retained
                        for i in indices_2_keep:
                            sample_geno = fields[i]
                            if sample_geno.startswith("./.") == False:
                                num_samples += 1
                        # now to add the columns
                        outline = ''
                        for i, column in enumerate(fields):
                            if i < 7 or i == 8:
                                outline += column + '\t'
                            elif i == 7:
                                ns_samples = f"NS={num_samples}"
                                for c in column.split(';')[1:]:
                                    ns_samples += ';' + c
                                outline += ns_samples + '\t'
                            elif i in indices_2_keep:
                                outline += fields[i] + '\t'
                        # replace last '\t'
                        outline = outline[:-1] + '\n'
                        outfh.write(outline)
                        break # snp recorded
    # close IOs
    fh.close()
    outfh.close()
    

def main():
    """entry point to the pipeline"""

    # get the input
    vcf, linkage_map, marker_coords, max_dist, samples_files = get_arguments()

    # get IDS from linkage map
    ids_map = get_linkage_markers(linkage_map)

    # get samples to retain
    samples = parse_samples_file(samples_files)

    # get locations of IDs
    marker_map = parse_marker_coordinates(marker_coords, ids_map)

    # now to create the new VCF
    filter_vcf(vcf, marker_map, max_dist, samples)

if __name__ == "__main__":
    main()