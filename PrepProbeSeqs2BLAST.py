#!/usr/bin/python3

"""
Script Name: PrepProbeSeqs2BLAST
Version: 1.0.0
Author: Thomas Bergmann
Description:
    This script takes as input the complete sequences from Illumina SNP chip data (inlcuding the SNP call in brackets)
    and converts it into a multifasta that is compatible with BLAST. The bracket will be replaced with the IUPAC code
    for the respective ambiguous SNP call. For example: GGTTT[A/C] --> GGTTTM
    The sequence will be reduced to a length of 50 bases with the SNP call being always at position 50. 
    Notice: BLAST will ignore the IUPAC code at position 50.

Usage:
    python PrepProbeSeqs2BLAST.py --input <input_file> --output <output_file>
    
Dependencies:
    - Python 3

"""

import sys
import os
import argparse

# Set global variables
count_back = 49 # Adjust this parameter to designated sequence length

# IUPAC mapping for forward variants
IUPAC_fw = {
    "R": "[A/G]",
    "Y": "[C/T]",
    "S": "[G/C]",
    "W": "[A/T]",
    "K": "[G/T]",
    "M": "[A/C]"
}

# IUPAC mapping for reverse variants
IUPAC_rev = {
    "R": "[G/A]",
    "Y": "[T/C]",
    "S": "[C/G]",
    "W": "[T/A]",
    "K": "[T/G]",
    "M": "[C/A]"
}

complement = {"A" : "T",
              "T" : "A",
              "C" : "G",
              "G" : "C",
              "R" : "Y",
              "Y" : "R",
              "S" : "S",
              "W" : "W",
              "K" : "M",
              "M" : "K",
              "N" : "N"}

# Argument Parsing
parser = argparse.ArgumentParser(description ="Convert Illumina SNP probe sequences into BLAST-compatible Multi-FASTA.")
parser.add_argument("-i", "--input", help = "Multi-FASTA with the probe sequences", required=True)
parser.add_argument("-o", "--output", help = "Ready-to-BLAST Multi-Fasta", required=True) 
args = parser.parse_args()

# Input file validation function
def check_file_existence(file_path):
    if not os.path.isfile(file_path):
        sys.exit(f"{file_path} doesn't exist! Please check your input files.")

check_file_existence(args.input)

# Function to loop through the Multi-FASTA and trim it
def trim_fasta(multifasta, count_back):
  # Initialize dictionaries
  seqDic = {}
  directionDic = {}
  counter_fw = 0
  counter_rev = 0
  # Start looping through the Multi-Fasta
  for line in multifasta:
    line = line.rstrip()
    # Remove leading and trailing whitespace
    if line.startswith(">"):
      identifier = line[1:]
      seqDic[identifier] = ""
      directionDic[identifier] = ""
    else:
      # Add the sequence itself
      start_bracket = line.find('[') # Start position of bracket
      end_bracket = line.find(']') + 1 # End position of bracket
      # Check whether you can count backwards
      if start_bracket >= count_back:
        counter_fw += 1
        start_seq = max(0, start_bracket - count_back) # Count 49 backwards from start of bracket
        trimmed_seq = line[start_seq:end_bracket] # Trim the sequence to 49 length + brackets
        seqDic[identifier] += trimmed_seq # Add trimmed seq to dictionary
        directionDic[identifier] += "forward"
      # Else count other direction
      else:
        counter_rev += 1
        end_seq = end_bracket + count_back
        trimmed_seq = line[start_bracket:end_seq] # Trim the sequence to 49 length + brackets
        seqDic[identifier] += trimmed_seq # Add trimmed seq to dictionary
        directionDic[identifier] += "reverse"
        
  return seqDic, directionDic, counter_fw, counter_rev

# Function to replace bracketed variants with IUPAC codes
def replace_brackets_with_iupac(seq, iupac_fw, iupac_rev):
  #Replace using forward variants
  for iupac_code, variant in iupac_fw.items():
    seq = seq.replace(variant, iupac_code)
  # Replace using reverse variants
  for iupac_code, variant in iupac_rev.items():
    seq = seq.replace(variant, iupac_code) 
  return seq

# Open the input file and process sequences
with open(args.input) as fasta_file:
    seqDic, directionDic, counter_fw, counter_rev = trim_fasta(fasta_file, count_back)

# Apply IUPAC replacement for each sequence in seqDic
for identifier, sequence in seqDic.items():
    # Replace the brackets with the iupac code for forward direction
  seqDic[identifier] = replace_brackets_with_iupac(sequence, IUPAC_fw, IUPAC_rev)
  # Replace the brackets with the iupac code for reverse direction
  if directionDic[identifier] == "reverse":
    seqDic[identifier] = seqDic[identifier].translate(str.maketrans(complement))[::-1]

# Write to new Multi-FASTA
with open(args.output, "w") as fasta_file:
  for identifier, sequence in seqDic.items():
    fasta_file.write(f">{identifier}\n{sequence}\n")

print("FASTA format output written to {}".format(args.output))
print("Number of forward sequences: {}".format(counter_fw))
print("Number of reverse sequences: {}".format(counter_rev))



