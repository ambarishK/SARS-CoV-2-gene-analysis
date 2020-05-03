
import Levenshtein



def location(str, substr):
    start = str.index(substr)
    end = start + len(substr) - 1
    return start, end


def find_boundaries(genome, sequence):
    start = 0
    if len(genome.split(sequence)) == 2:
        start, end = location(genome, sequence)


    else:
        temp_min = 10000
        temp_start = 0

        for i in range(len(genome) - len(sequence)):
            if Levenshtein.distance(sequence, genome[i:i + len(sequence)]) < temp_min:
                temp_min = Levenshtein.distance(sequence, genome[i:i + len(sequence)])
                temp_start = i
        start = temp_start
        print(genome[start:start + len(sequence)])

        genome = genome.split(genome[start:start + len(sequence)])[0] + sequence + \
                     genome.split(genome[start:start + len(sequence)])[1]

        print(sequence)
    return (start, start + len(sequence)), genome


ref_genome = ''
with open('./ReferenceGenes/primers_reference_genome.fasta', 'r') as ref_genome_file:
    for line in ref_genome_file:
        if not (line.startswith('>')):
            ref_genome += line.rstrip()




with open('./ReferenceGenes/primers_public.fas', 'r') as ref_primers_file, open('./ReferenceGenes/primers.txt', 'w') as primers_txt, open('./ReferenceGenes/primers_reference_genome.fasta', 'w') as new_genome_file:
    names = list()
    primers = list()
    for line in ref_primers_file:
        if line.startswith('>'):
            names.append(line.rstrip())
        else:
            primers.append(line.rstrip())
    primers = dict(zip(names, primers))
    primers_locations = dict()
    for primer_name in primers:
        print(primer_name)
        primers_locations[primer_name], ref_genome = find_boundaries(ref_genome, primers[primer_name])
    find_boundaries(ref_genome, 'ATATTGCAGCAGTACGCACACA')
    for item in primers_locations:

        primers_txt.write(
            str(item[1:]) + ' ' + str(primers_locations[item][0]) + ' ' + str(primers_locations[item][1]) + '\n')
    new_genome_file.write('>MN908947.3 SARS-COV2 genome changed to compare primers' + '\n')
    new_genome_file.write(ref_genome)
#
# #
# #
# #
# # def find_genes(genome, boundaries, tolerance, reference_genes, gene_names, ref_genome):
# #     for bounds in boundaries:
# #         bounds[0] -= tolerance
# #         bounds[1] += tolerance
# #         if bounds[0] < 0:
# #             bounds[0] = 0
# #         if bounds[1] > len(genome):
# #             bounds[1] = len(genome)
# #
# #     genes = list()
# #
# #     for (interval, ref, name) in zip(boundaries, reference_genes, gene_names):
# #         ref_prefix = ref[0:64]
# #         ref_suffix = ref[len(ref) - 64:len(ref)]
# #         seq = genome[interval[0]:interval[1]]
# #         # Finding potential gene prefixes starting with start codon and suffixes ending with stop codon.
# #         prefixes = [m.start() for m in re.finditer('(?=ATG.{61})', seq)]
# #         suffixes = [m.start() for m in re.finditer('(?=.{61}(TGA|TAG|TAA))', seq)]
# #         # Finding starting index of prefix and suffix with minimal edit distance to reference gene.
# #         if len(prefixes) > 0 and len(suffixes) > 0:
# #             minPrefix = min(prefixes, key=lambda x: Levenshtein.distance(seq[x:x + 64], ref_prefix))
# #             minSuffix = min(suffixes, key=lambda x: Levenshtein.distance(seq[x:x + 64], ref_suffix))
# #             found_gene = seq[minPrefix:minSuffix + 64]
# #             ref_protein = translate_rna(transcribe(ref))
# #             found_protein = translate_rna(transcribe(found_gene))
# #             genes.append(
# #                 [name, interval[0] + minPrefix, interval[0] + minSuffix + 64, Levenshtein.distance(ref, found_gene)])
# #             genes.append([name + '_translation', 0,
# #                           0, Levenshtein.distance(found_protein, ref_protein)])
# #
# #             genes.append([name + '_N', 0, 0, len(re.findall('[^ACTG]', found_gene))])
# #             mutations = []
# #             if len(found_gene) <= len(ref):
# #                 genetic_len = len(found_gene)
# #                 protein_len = len(found_protein)
# #             else:
# #                 genetic_len = len(ref)
# #                 protein_len = len(ref_protein)
# #             for i in range(genetic_len):
# #                 if found_gene[i] != ref[i]:
# #                     mutations.append([i, ref[i], found_gene[i]])
# #
# #             protein_changes = []
# #             for i in range(protein_len):
# #                 if found_protein[i] != ref_protein[i]:
# #                     protein_changes.append([i, ref_protein[i], found_protein[i]])
# #             genes.append([name + '_mutations', 0, 0, mutations])
# #             genes.append([name + '_translation_changes', 0, 0, protein_changes])
# #
# #     edit_distance = Levenshtein.distance(genome, ref_genome)
# #     return genes, edit_distance
# #
# #
# # def get_reference_gene():
# #     boundaries = list()
# #     reference_genes = list()
# #     gene_names = list()
# #     ref_genome = ''
# #
# #     with open('./ReferenceGenes/genes.txt', 'r') as ref_info_file, open('./ReferenceGenes/REFERENCE_GENOME.fasta',
# #                                                                         'r') as ref_genome_file:
# #         for line in ref_genome_file:
# #             if not (line.startswith('>')):
# #                 ref_genome += line.rstrip()
# #
# #         for line in ref_info_file:
# #             ref_gene_info = line.split()
# #             boundaries.append([int(ref_gene_info[1]) - 1, int(ref_gene_info[2])])
# #             reference_genes.append(ref_genome[int(ref_gene_info[1]) - 1:int(ref_gene_info[2])])
# #             gene_names.append(ref_gene_info[0])
# #
# #     return ref_genome, boundaries, reference_genes, gene_names
# #
# #
# # def count_all():
# #     ref_genome, boundaries, reference_genes, gene_names = get_reference_gene()
# #
# #     genomes = []
# #     genomes_headers = []
# #
# #     time_start = time.time()
# #
# #     with open('Cleaned_up_genes.fasta', 'r') as f:
# #         genome = ''
# #         genome_header = ''
# #         for line in tqdm(f):
# #             if line.startswith('>'):
# #                 if len(genome_header) > 0 and len(genome) > 29000:
# #                     genomes_headers.append(genome_header)
# #                     genomes.append(genome.upper())
# #                 genome = ''
# #                 genome_header = line
# #             else:
# #                 genome += line.rstrip()
# #
# #     print("We have in total {} genomes collected from the fasta file".format(len(genomes)))
# #
# #     try:
# #         pool = Pool(cpu_count())
# #         genes_and_distances = pool.map(partial(find_genes,
# #                                                boundaries=boundaries,
# #                                                tolerance=128,
# #                                                reference_genes=reference_genes,
# #                                                gene_names=gene_names, ref_genome=ref_genome), genomes)
# #     finally:
# #         pool.close()
# #         pool.join()
# #
# #     all_genes, edit_distances = zip(*genes_and_distances)
# #
# #     print("In total we have parsed {} genomes".format(len(all_genes)))
# #
# #     with open('distances.txt', 'w') as output_file:
# #         for genome_header, genes, edit_distance, genome in tqdm(
# #                 zip(genomes_headers, all_genes, edit_distances, genomes)):
# #             output_file.write(genome_header)
# #             genes.append(['whole_genome', 0, len(genome), edit_distance])
# #             output_file.write(str(genes) + '\n')
# #
# #     print("Total time {}".format(time.time() - time_start))
# #
#
# # if __name__ == '__main__':
# #     count_all()
