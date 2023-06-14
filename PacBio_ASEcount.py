#!/usr/bin/env python
# Program: This is a python script for allele-specific junction count in PacBio long reads.
# Date: 20230614
# Input: PacBio BAM file and lead sSNP list
# Output: Allele-specific junction count for each sSNP

# set working directory path
os.chdir('/ebs1/users/wangwenjing/ASE_vertify/make_table')
import pysam
import re

# Function: Check the base at targeted position
def get_base_at_position(read_sequence, cigar, start_position, target_position):
    cigar_ops = re.findall(r'(\d+)(\D)', cigar)
    cigar_ops = [(int(length), op) for length, op in cigar_ops]

    read_pointer = 0
    ref_pointer = start_position
    target_position -= 1 # (Transfer 1-Base to 0-Base)

    for length, op in cigar_ops:
        if op in ['M', '=', 'X']:

            if ref_pointer <= target_position < ref_pointer + length:
                read_offset = target_position - ref_pointer
                return read_sequence[read_pointer + read_offset]
            read_pointer += length
            ref_pointer += length

        elif op == 'I':
            read_pointer += length

        elif op in ['D', 'N']:
            ref_pointer += length

        elif op == 'S':
            read_pointer += length

    return 'N'

# Function: Check intron with perfect match
def check_intron(cigar, read_start, intron_start, intron_end):
    pattern = re.compile(r'(\d+)(\D)')
    ops = pattern.findall(cigar)
    current_position = read_start

    intron_start -= 1
    intron_end -= 1

    for length, op in ops:
        length = int(length)

        if op == 'M' or op == '=' or op == 'X':
            current_position += length
        elif op == 'D' or op == 'N':
            if op == 'N' and current_position == intron_start and current_position + length == intron_end:
                return True
            current_position += length
        elif op == 'I' or op == 'S' or op == 'H':
            pass
        else:
            print(f"Unknown CIGAR operation: {op}")
            break
    return False

# Process file
def process_file(txt_file, bam_file, output_file):
    with open(txt_file, 'r') as txt, open(output_file, 'w') as out:
        bam = pysam.AlignmentFile(bam_file, 'rb')

        out.write('snp_id\tref_with_junction\tref_without_junction\talt_with_junction\talt_without_junction\n')

        for line in txt:
            fields = line.strip().split('\t')

            snp_id = fields[0]
            snp_chr = fields[1]
            snp_pos = int(fields[2])
            intron_start = int(fields[5])
            intron_end = int(fields[6])

            _, _, ref_allele, alt_allele = snp_id.split(':')

            ref_with_junction = ref_without_junction = alt_with_junction = alt_without_junction = 0
            
            for read in bam.fetch(snp_chr, snp_pos, snp_pos+1):

                base_at_snp = get_base_at_position(read.query_sequence,read.cigarstring, read.reference_start, snp_pos)
                contains_intron = check_intron(read.cigarstring, read.reference_start, intron_start, intron_end)

                if base_at_snp == ref_allele:
                    if contains_intron:
                        ref_with_junction += 1
                    else:
                        ref_without_junction += 1
                elif base_at_snp == alt_allele:
                    if contains_intron:
                        alt_with_junction += 1
                    else:
                        alt_without_junction += 1

            out.write(f'{snp_id}\t{ref_with_junction}\t{ref_without_junction}\t{alt_with_junction}\t{alt_without_junction}\n')

        bam.close()

# Main
if __name__ == '__main__':
    txt_file = '/ebs1/users/wangwenjing/ASE_vertify/lead_sSNP_list/SNP_list/organized_lead_SNP/filtered_to_exon/merge_info/pDC_merged.txt'
    bam_file = '/ebs1/users/wangwenjing/barcodes_for_celltype/bam_for_ct/original_sorted/pDC_sorted.bam'
    output_file = '/ebs1/users/wangwenjing/ASE_vertify/make_table/table_for_cells/pDC_for_ase.csv'

    process_file(txt_file, bam_file, output_file)
