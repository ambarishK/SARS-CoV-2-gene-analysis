import os
import re
import pylcs
import itertools
def transcribe(sequence):
    return sequence.replace('T', 'U')

def translate_rna(s):
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N",
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S",
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I",

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H",
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R",
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D",
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G",
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V",

                "UAA":"_", "UAC":"Y", "UAG":"_", "UAU":"T",
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S",
                "UGA":"_", "UGC":"C", "UGG":"W", "UGU":"C",
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}

    l = [codon2aa.get(s[n:n+3], 'X') for n in range(0, len(s), 3)]
    return "".join(l)

def find_genes(genome, boundaries, tolerance, reference_genes, gene_names):
    for bounds in boundaries:
        bounds[0] -= tolerance
        bounds[1] += tolerance
        if bounds[0] < 0:
            bounds[0] = 0
        if bounds[1] > len(genome):
            bounds[1] = len(genome)

    genes = list()

    for (interval, ref, name) in zip(boundaries, reference_genes, gene_names):
        ref_prefix = ref[0:64]
        ref_suffix = ref[len(ref)-64:len(ref)]
        seq = genome[interval[0]:interval[1]]
        # Finding potential gene prefixes starting with start codon and suffixes ending with stop codon.
        prefixes = [m.start() for m in re.finditer('(?=ATG.{61})', seq)]
        suffixes = [m.start() for m in re.finditer('(?=.{61}(TGA|TAG|TAA))', seq)]
        # Finding starting index of prefix and suffix with minimal edit distance to reference gene.
        minPrefix = min(prefixes, key=lambda x : pylcs.edit_distance(seq[x:x+64], ref_prefix))
        minSuffix = min(suffixes, key=lambda x : pylcs.edit_distance(seq[x:x+64], ref_suffix))
        genes.append([name, interval[0] + minPrefix, interval[0] + minSuffix + 64, pylcs.edit_distance(ref, seq[minPrefix:minSuffix + 64])])
        genes.append([name, interval[0] + minPrefix, interval[0] + minSuffix + 64, pylcs.edit_distance(translate_rna((transcribe(ref))), translate_rna(transcribe((seq[minPrefix:minSuffix + 64]))))])

    return genes

def get_reference_gene():
    boundaries = list()
    reference_genes = list()
    gene_names = list()
    ref_genome = ''

    with open('./ReferenceGenes/genes.txt', 'r') as ref_info_file, open('./ReferenceGenes/REFERENCE_GENOME.fasta', 'r') as ref_genome_file:
        for line in ref_genome_file:
            if not (line.startswith('>')):
                ref_genome += line.rstrip()

        for line in ref_info_file:
            ref_gene_info = line.split()
            boundaries.append([int(ref_gene_info[1]) - 1, int(ref_gene_info[2])])
            reference_genes.append(ref_genome[int(ref_gene_info[1]) - 1:int(ref_gene_info[2])])
            gene_names.append(ref_gene_info[0])

    return ref_genome, boundaries, reference_genes, gene_names


if __name__ == '__main__':

    ref_genome, boundaries, reference_genes, gene_names = get_reference_gene()

    with open('gisaid_cov2020_sequences.fasta', 'r') as f, open('distances.txt', 'w') as output_file:
        genome = ''
        genome_header = ''
        for line in f:
            if line.startswith('>'):
                if len(genome_header) > 0 and len(genome) > 29000:
                    output_file.write(genome_header)
                    genome = genome.upper()
                    genes = find_genes(genome, boundaries, 128, reference_genes, gene_names)
                    genes.append(['whole_genome', 0, len(genome), pylcs.edit_distance(genome, ref_genome)])
                    output_file.write(str(genes) + '\n')
                genome=''
                genome_header = line
            else:
                genome += line.rstrip()
     
