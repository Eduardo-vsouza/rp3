import os
import sys
import multiprocessing
from tqdm import tqdm
import pickle

import numpy as np
import pysam
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure
from ..utils import GTFFilter


def count_alignments_for_read(read, mmExons):
    gene_counts = {}
    for gene_id, exons in mmExons.items():
        for exon in exons:
            if read.reference_name == exon['chrom'] and \
               exon['start'] <= read.reference_start <= exon['end']:
                gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1
                break
        else:
            continue
        break
    return read.query_name, gene_counts

class MMCutoff(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.multimappersGTF = f'{self.countsDir}/multimapping_smorfs.gtf'

        self.mmExons = self.__get_exon_information()
        self.exonCoverage = self.__set_coverage_scaffold()
        self.mmSimulationDir = f'{self.riboSeqDir}/MM_simulation'
        self.intersectedAlignmentsDir = f'{self.mmSimulationDir}/intersected_alignments'
        self.check_dirs([self.mmSimulationDir, self.intersectedAlignmentsDir])
        self.__filter_multimappers_gtf()
        self.mappingsPickle = f'{self.intersectedAlignmentsDir}/mappings.pickle'

    def __filter_multimappers_gtf(self):
        gtf = GTFFilter(args=self.args, gtf=self.rescoredMicroproteinsGTF, fasta=self.microproteinsMM)
        gtf.get_entries()
        gtf.filter_gtf(output=self.multimappersGTF)

    def sort_bam_files(self):
        files = os.listdir(self.riboSeqAlnDir)
        for file in files:
            out = f'{self.riboSeqSortedBamDir}/{file[:-4]}_sorted.bam'
            if file.endswith("out.sam"):
                if not os.path.exists(out) or self.args.overwrite:
                    cmd = (f'samtools sort -@ {self.args.threads} -o {out} '
                           f'{self.riboSeqAlnDir}/{file}')
                    os.system(cmd)

    def intersect_alignments(self):
        print(f"--Intersecting alignment files with multimapping smORFs.")
        files = os.listdir(self.riboSeqSortedBamDir)
        for file in tqdm(files):
            if file.endswith(".bam"):
                intersected = f'{self.intersectedAlignmentsDir}/{file[:-4]}_intersected.bam'
                if not os.path.exists(intersected) or self.args.overwrite:
                    cmd = (f'bedtools intersect -u -f 1 -a {self.riboSeqSortedBamDir}/{file} -b {self.multimappersGTF} >'
                           f' {intersected}')
                    os.system(cmd)

    def index_intersected_bam_files(self):
        print(f"--Indexing intersected bam files.")
        files = os.listdir(self.intersectedAlignmentsDir)
        for file in tqdm(files):
            if file.endswith("bam"):
                if not os.path.exists(f'{self.intersectedAlignmentsDir}/{file}.bai') or self.args.overwrite:
                    cmd = f'samtools index {self.intersectedAlignmentsDir}/{file}'
                    os.system(cmd)

    def get_alignments(self):
        print(f"Obtaining the number of alignments for each read that mapped to a multimapping smORF.")
        intersected_files = os.listdir(self.intersectedAlignmentsDir)
        gtf_bed = pybedtools.BedTool(self.multimappersGTF)
        gene_read_counts = {}
        i = 0
        for bam in intersected_files:
            if bam.endswith(".bam"):
                try:
                    print(f"Processing bam file {i} out of {len(intersected_files)}")
                    bam_file = pysam.AlignmentFile(f"{self.intersectedAlignmentsDir}/{bam}", "rb",
                                                   threads=self.args.threads)
                    if bam not in gene_read_counts:
                        gene_read_counts[bam] = {}
                    for alignment in tqdm(bam_file):
                        # Check if the alignment is mapped
                        if not alignment.is_unmapped:
                            # Get the read name and alignment position
                            read_name = alignment.query_name
                            chrom = alignment.reference_name
                            start = alignment.reference_start
                            end = alignment.reference_end

                            # Create a BedTool object for the alignment
                            alignment_bed = pybedtools.BedTool(f"{chrom}\t{start}\t{end}\t{read_name}", from_string=True)

                            # Intersect the alignment with the GTF file
                            intersected_features = alignment_bed.intersect(gtf_bed, wo=True)

                            # Iterate through intersected features
                            for feature in intersected_features:

                                # Extract gene ID from the GTF attributes
                                gene_id = feature.fields[-2].split(" ")[-1].replace("\"", "").replace(";", "")

                                # read_name = alignment.query_name
                                nh_tag = alignment.get_tag("NH")  # NH holds the number of alignments for that read
                                # mismatches = alignment.get_tag("nM")
                                # if mismatches <= self.args.maxMismatches:
                                    # print("aln", alignment)
                                    # print("nh", nh_tag)

                                if gene_id in gene_read_counts[bam]:
                                    gene_read_counts[bam][gene_id]['nh'].append(nh_tag)
                                    gene_read_counts[bam][gene_id]['coords'].append(f'{start}-{end}')

                                else:
                                    gene_read_counts[bam][gene_id] = {}
                                    gene_read_counts[bam][gene_id]['nh'] = [nh_tag]
                                    gene_read_counts[bam][gene_id]['coords'] = [f'{start}-{end}']

                except:
                    pass
                # break
            # break
        with open(self.mappingsPickle, "wb") as pickle_file:
            pickle.dump(gene_read_counts, pickle_file)

    def order_mappings(self):
        with open(self.mappingsPickle, "rb") as pickle_file:
            self.geneMappings = pickle.load(pickle_file)
        alns = {'gene': [], 'min_mm': [], 'file': []}
        sorted_gene_n = {}
        j = 1
        total_mm = {}
        sorted_gene = []
        coverage = {}
        for bam in self.geneMappings:
            if bam not in alns:
                alns[bam] = {}
            for gene in self.geneMappings[bam]:
                nh_coords_zip = zip(self.geneMappings[bam][gene]['nh'], self.geneMappings[bam][gene]['coords'])

                # Sort based on 'nh' values
                sorted_nh_coords = sorted(nh_coords_zip, key=lambda x: x[0])

                # Unzip the sorted pairs
                sorted_nh, sorted_coords = zip(*sorted_nh_coords)
                print(sorted_nh[:15], gene)
                # Update 'nh' and 'coords' lists in the dictionary
                self.geneMappings[bam][gene]['nh'] = list(sorted_nh)
                self.geneMappings[bam][gene]['coords'] = list(sorted_coords)
                # count =
                # self.geneMappings[bam][gene].sort()
                # print(self.geneMappings[bam][gene][:20])
                if len(self.geneMappings[bam][gene]['coords']) >= self.args.minRawCounts:
                    # if self.geneMappings[bam][gene]['nh'][self.args.minRawCounts-1] >= self.args.minRawCounts:
                        # print("candidate", gene, self.geneMappings[bam][gene]['nh'])
                    # print(self.geneMappings[bam][gene])
                    # print(len(self.geneMappings[bam][gene]))
                    alns['gene'].append(gene)
                    alns['min_mm'].append(self.geneMappings[bam][gene]['nh'][self.args.minRawCounts-1])
                    alns['file'].append(bam)
                    # if self.geneMappings[bam][gene][self.args.minRawCounts-1] > self.args.minRawCounts:
                    if gene not in sorted_gene_n:
                        sorted_gene_n[gene] = j
                        total_mm[gene] = []
                        j += 1
                    sorted_gene.append(sorted_gene_n[gene])
                    # total_mm[gene].append(self.geneMappings[bam][gene][self.args.minRawCounts-1])
                        # print(self.geneMappings[bam][gene])
                        # print(gene, bam)
            # break
        # print(self.geneMappings[bam][])
            # break
        zipped_values = zip(alns['gene'], alns['min_mm'], alns['file'])
        # mm_norm = []
        # genes = []
        #
        # def geo_mean(iterable):
        #     a = np.array(iterable)
        #     return a.prod() ** (1.0 / len(a))
        #
        # for gene in alns['gene']:
        #     total_mm[gene].sort()
        #     # print(total_mm[gene])
        #     if gene not in genes:
        #         mm_norm.append(np.median(total_mm[gene]))
        #         genes.append(gene)
        # Sort the zipped values based on the values of 'min_mm'
        # sorted_zipped_values = sorted(zipped_values, key=lambda x: x[1])

        # Separate the sorted values back into 'gene' and 'min_mm'
        # sorted_gene, sorted_min_mm, sorted_file = zip(*sorted_zipped_values)
        # sorted_alns = sorted(zip(alns['min_mm'], alns['gene'], alns['file']))
        # sorted_alns = sorted(zip(alns['min_mm'], sorted_gene, alns['file']))
        sorted_alns = sorted(zip(alns['min_mm'], sorted_gene, alns['file']))
        # sorted_a = sorted(zip(mm_norm, genes))
        # print(sorted_a)
        # Extract keys and values back into separate lists in the sorted order
        # sorted_min_mm = [item[0] for item in sorted_alns]
        # sorted_gene_names = [item[1] for item in sorted_alns]
        # sorted_file = [item[2] for item in sorted_alns]

        # sorted_gene_i = {}
        # sorted_gene = []
        # data = {'min_mm': mm_norm, 'gene': genes}
        # print(sorted_gene_i)
        # print(sorted_gene)
        # Create the sorted dictionary
        # print(list(sorted_min_mm))
        # sorted_indices = sorted(range(len(genes)), key=lambda k: genes[k])
        # sorted_genes = [genes[i] for i in sorted_indices]
        # sorted_mm_norm = [mm_norm[i] for i in sorted_indices]
        # sorted_alns = {'gene': list(sorted_gene_names), 'min_mm': list(sorted_min_mm), 'file': list(sorted_file)}
        # ax = sns.lineplot(data=sorted_alns, x="gene", y="min_mm", hue='file')
        # ax = sns.lineplot(x="gene", y="min_mm")
        # data = {'Gene': genes, 'mm_norm': mm_norm}
        # df = pd.DataFrame(data)
        #
        # # Sort the DataFrame by 'Gene'
        # df_sorted = df.sort_values(by='mm_norm')
        #
        # # Plot the sorted data using seaborn
        # sns.lineplot(data=df_sorted, x='Gene', y='mm_norm', marker='o')
        # plt.legend([])
        #
        # plt.show()


    def get_coverages(self):
        with open(self.mappingsPickle, "rb") as pickle_file:
            self.geneMappings = pickle.load(pickle_file)
        alns = {'gene': [], 'min_mm': [], 'file': []}
        sorted_gene_n = {}
        j = 1
        total_mm = {}
        sorted_gene = []
        coverage = {}
        for bam in self.geneMappings:
            j += 1
            if bam not in alns:
                alns[bam] = {}
            for gene in self.geneMappings[bam]:
                if len(self.geneMappings[bam][gene]['coords']) >= self.args.minRawCounts:

                    nh_coords_zip = zip(self.geneMappings[bam][gene]['nh'], self.geneMappings[bam][gene]['coords'])

                    # Sort based on 'nh' values
                    sorted_nh_coords = sorted(nh_coords_zip, key=lambda x: x[0])

                    # Unzip the sorted pairs
                    sorted_nh, sorted_coords = zip(*sorted_nh_coords)

                    # Update 'nh' and 'coords' lists in the dictionary
                    self.geneMappings[bam][gene]['nh'] = list(sorted_nh)
                    self.geneMappings[bam][gene]['coords'] = list(sorted_coords)
                    # count =
                    # self.geneMappings[bam][gene].sort()
                    for nh, coords in zip(self.geneMappings[bam][gene]['nh'], self.geneMappings[bam][gene]['coords']):
                        # print(nh, coords)
                        self.__update_coverage(gene=gene, nh=nh, coords=coords)

            if j > 2:
                break
        with open(f'{self.riboSeqDir}/exon_coverage.pickle', "wb") as pickle_file:
            pickle.dump(self.exonCoverage, pickle_file)
        # print(self.exonCoverage)
                    # print(self.geneMappings[bam][gene][:20])
                    # print(self.geneMappings[bam][gene])
                    # print(len(self.geneMappings[bam][gene]))
                    # alns['gene'].append(gene)
                    # alns['min_mm'].append(self.geneMappings[bam][gene][self.args.minRawCounts-1])
                    # alns['file'].append(bam)
                    # # if self.geneMappings[bam][gene][self.args.minRawCounts-1] > self.args.minRawCounts:
                    # if gene not in sorted_gene_n:
                    #     sorted_gene_n[gene] = j
                    #     total_mm[gene] = []
                    #     j += 1
                    # sorted_gene.append(sorted_gene_n[gene])
                    # total_mm[gene].append(self.geneMappings[bam][gene][self.args.minRawCounts-1])
                    # print(self.geneMappings[bam][gene])
                    # print(gene, bam)
        # print(self.geneMappings[bam][])
        #     # break
        # # zipped_values = zip(alns['gene'], alns['min_mm'], alns['file'])
        # mm_norm = []
        # genes = []
        #
        # def geo_mean(iterable):
        #     a = np.array(iterable)
        #     return a.prod() ** (1.0 / len(a))
        #
        # for gene in alns['gene']:
        #     total_mm[gene].sort()
        #     # print(total_mm[gene])
        #     if gene not in genes:
        #         mm_norm.append(np.median(total_mm[gene]))
        #         genes.append(gene)
        # # Sort the zipped values based on the values of 'min_mm'
        # # sorted_zipped_values = sorted(zipped_values, key=lambda x: x[1])
        #
        # # Separate the sorted values back into 'gene' and 'min_mm'
        # # sorted_gene, sorted_min_mm, sorted_file = zip(*sorted_zipped_values)
        # # sorted_alns = sorted(zip(alns['min_mm'], alns['gene'], alns['file']))
        # # sorted_alns = sorted(zip(alns['min_mm'], sorted_gene, alns['file']))
        # sorted_alns = sorted(zip(alns['min_mm'], sorted_gene, alns['file']))
        # # sorted_a = sorted(zip(mm_norm, genes))
        # # print(sorted_a)
        # # Extract keys and values back into separate lists in the sorted order
        # # sorted_min_mm = [item[0] for item in sorted_alns]
        # # sorted_gene_names = [item[1] for item in sorted_alns]
        # # sorted_file = [item[2] for item in sorted_alns]
        #
        # # sorted_gene_i = {}
        # # sorted_gene = []
        # data = {'min_mm': mm_norm, 'gene': genes}
        # # print(sorted_gene_i)
        # # print(sorted_gene)
        # # Create the sorted dictionary
        # # print(list(sorted_min_mm))
        # # sorted_indices = sorted(range(len(genes)), key=lambda k: genes[k])
        # # sorted_genes = [genes[i] for i in sorted_indices]
        # # sorted_mm_norm = [mm_norm[i] for i in sorted_indices]
        # # sorted_alns = {'gene': list(sorted_gene_names), 'min_mm': list(sorted_min_mm), 'file': list(sorted_file)}
        # # ax = sns.lineplot(data=sorted_alns, x="gene", y="min_mm", hue='file')
        # # ax = sns.lineplot(x="gene", y="min_mm")
        # data = {'Gene': genes, 'mm_norm': mm_norm}
        # df = pd.DataFrame(data)
        #
        # # Sort the DataFrame by 'Gene'
        # df_sorted = df.sort_values(by='mm_norm')
        #
        # # Plot the sorted data using seaborn
        # sns.lineplot(data=df_sorted, x='Gene', y='mm_norm', marker='o')
        # plt.legend([])
        #
        # plt.show()

    def __get_exon_information(self):
        exons = {}
        with open(self.multimappersGTF, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split('\t')
                if cols [2] == 'exon':
                    attrs = cols[8].split(";")
                    for a in attrs:
                        if 'transcript_id' in a:
                            rna = a.split(" ")[-1].replace("\"", "")
                            start, end = cols[3], cols[4]
                            if rna not in exons:
                                exons[rna] = {}
                                # exons[rna] = {'coverage': []}
                                exons[rna]['coords'] = [(start, end)]
                            else:
                                exons[rna]['coords'].append((start, end))
        return exons

    def __set_coverage_scaffold(self):
        exon_coverage = {}
        for rna in self.mmExons:
            if rna not in exon_coverage:
                exon_coverage[rna] = {}
            for coords in self.mmExons[rna]['coords']:
                # print(rna)
                # print(self.mmExons[rna])
                # print(self.mmExons[rna]['coords'])
                # print(coords)
                start, end = int(coords[0]), int(coords[1])
                for i in range(start, end):
                    exon_coverage[rna][i] = {'cov': 0, 'nh': 0}
                    # exon_coverage[rna][i]['cov'] = 0
                    # exon_coverage[rna][i]['nh'] = 0
        return exon_coverage

    def __update_coverage(self, gene, nh, coords):
        # print(gene, nh, coords)
        start, end = int(coords.split("-")[0]), int(coords.split("-")[1])
        for i in range(start, end):
            add = False
            for rna in self.mmExons:
                if '_F:' in rna:
                    # print(self.mmExons[rna]['coords'])
                    for coord in self.mmExons[rna]['coords']:
                        starts, ends = int(coord[0]), int(coord[1])
                        # print(starts, ends)
                        if i in range(starts, ends):
                            add = True
                            break
            if add:
                if self.exonCoverage[gene][i]['cov'] == 0:
                    self.exonCoverage[gene][i]['cov'] += 1
                    self.exonCoverage[gene][i]['nh'] = nh
                else:
                    self.exonCoverage[gene][i]['cov'] += 1


    def calculate_min_mm_for_coverage(self):
        with open(f'{self.riboSeqDir}/exon_coverage.pickle', "rb") as pickle_file:
            self.exonCoverage = pickle.load(pickle_file)
        min_mm = {}
        min_cov = 100
        for gene in self.exonCoverage:
            nh = 0
            if gene not in min_mm:
                min_mm[gene] = {'coverage': 0, 'min_mm': 0}
            # print(self.exonCoverage)
            for i in self.exonCoverage[gene]:
                # if len(self.exonCoverage[])
                if min_mm[gene]['coverage'] < min_cov:
                    print(len(self.exonCoverage[gene]))
                    print(100/len(self.exonCoverage[gene]))

                    if self.exonCoverage[gene][i]['cov'] > 0:
                        min_mm[gene]['coverage'] += (100/len(self.exonCoverage[gene]))
                        if self.exonCoverage[gene][i]['nh'] > nh:
                            nh = self.exonCoverage[gene][i]['nh']
            min_mm[gene]['min_mm'] = nh

        data = {'genes': [], 'min_mm': [], 'cov': []}
        for gene in min_mm:
            # if min_mm[gene]['coverage'] >= min_cov:
            data['cov'].append(min_mm[gene]['coverage'])
            data['min_mm'].append(min_mm[gene]['min_mm'])
            data['genes'].append(gene)

        nh_gene_zip = zip(data['min_mm'], data['genes'])

        # Sort based on 'nh' values
        sorted_nh_gene = sorted(nh_gene_zip, key=lambda x: x[0])

        # Unzip the sorted pairs
        sorted_nh, sorted_gene = zip(*sorted_nh_gene)
        ndata = {'genes': sorted_gene, 'min_mm': sorted_nh}
        # ax = sns.lineplot(data=data, x='genes', y='min_mm', marker='o')
        ax = sns.scatterplot(data=data, x='cov', y='min_mm')

        plt.show()


    def extract_exons(self):
        with open(self.args.rescoredMicroproteinsGTF, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'exon':
                        attrs = fields[8].split(";")
                        for a in attrs:
                            if 'transcript_id' in a:
                                smorf = a.split(" ")[-1].replace("\"", "")
                        if smorf in self.multimappers:
                            gene_id = fields[8].split(';')[0].split('"')[1]
                            exon_info = {'chrom': fields[0], 'start': int(fields[3]), 'end': int(fields[4])}
                            self.mmExons.setdefault(gene_id, []).append(exon_info)

    def count_alignments(self):
        alignment_count = {}
        total_alignment_count = {}
        files = os.listdir(self.args.riboSeqSortedBamDir)
        bam_files = [file for file in files if file.endswith(".bam")]
        for file in bam_files:
            with pysam.AlignmentFile(os.path.join(self.args.riboSeqSortedBamDir, file), 'rb') as sam:
                pool = multiprocessing.Pool(processes=128)
                results = pool.starmap(count_alignments_for_read, [(read, self.mmExons) for read in sam.fetch()])
            for read_id, gene_counts in results:
                alignment_count[read_id] = gene_counts
        print(alignment_count)

        # return alignment_count, total_alignment_count
        # for file in bam_files:
        #     with pysam.AlignmentFile(f'{self.riboSeqSortedBamDir}/{file}', 'rb') as sam:
        #         for read in sam.fetch():
        #             total_alignment_count[read.query_name] = total_alignment_count.get(read.query_name,
        #                                                                                0) + 1  # Counting total alignments
        #             for gene_id, exons in self.mmExons.items():
        #                 for exon in exons:
        #                     if read.reference_name == exon['chrom'] and \
        #                             exon['start'] <= read.reference_start <= exon['end']:
        #                         alignment_count[read.query_name] = alignment_count.get(read.query_name, 0) + 1
        #                         break  # Once a read is found to map to an exon, no need to check other exons of the gene
        #                 else:
        #                     continue  # Continue to the next gene
        #                 break  # Once a read is found to map to an exon, no need to check other genes
        #     break
