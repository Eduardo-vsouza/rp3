import os
import sys

import pysam
import pyranges
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from ..pipeline_config import PipelineStructure


class SeqDist(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.orfCoverage = None
        self.genomeCoverage = {}
        self.regions = self.__parse_gtf(self.rescoredMicroproteinsGTF)

        self.smorfCoverages = {}
        self.singleCoverage = {}
        self.covPlots = f'{self.riboSeqCovDir}/coverage_plots'
        self.check_dirs([self.covPlots])

        # self.genes = self.__parse_gtf(f'{self.rescoredMicroproteinsGTF}')
    #
    # @staticmethod
    # def __parse_gtf(gtf_file):
    #     genes = {}
    #     with open(gtf_file) as f:
    #         for line in f:
    #             if line.startswith('#'):
    #                 continue
    #             parts = line.strip().split('\t')
    #             if parts[2] == 'transcript':
    #                 attributes = dict(item.strip().split(' ') for item in parts[8].strip().split(';') if item.strip())
    #                 gene_id = attributes['gene_id'].strip('"')
    #                 chrom = parts[0]
    #                 start = int(parts[3])
    #                 end = int(parts[4])
    #                 genes[gene_id] = {'chrom': chrom, 'start': start, 'end': end}
    #     return genes
    #
    # # Function to associate read coverage data with genes
    # @staticmethod
    # def associate_coverage_with_genes(args):
    #     read_coverage_data, gene_id, gene_info = args
    #     gene_coverage = 0
    #     for (chrom, pos), coverage in read_coverage_data.items():
    #         if gene_info['chrom'] == chrom and gene_info['start'] <= pos <= gene_info['end']:
    #             gene_coverage += coverage
    #     return gene_id, gene_coverage
    #
    # @staticmethod
    # def process_aligned_read(args):
    #     read_info, genome_coverage = args
    #     if not read_info['is_unmapped']:
    #         start_pos = read_info['reference_start'] + 1  # Add 1 to convert 0-based to 1-based position
    #         end_pos = read_info['reference_end'] + 1  # Add 1 to include the end position
    #         for pos in range(start_pos, end_pos):
    #             genome_coverage[(read_info['reference_name'], pos)] = genome_coverage.get(
    #                 (read_info['reference_name'], pos), 0) + 1

    @staticmethod
    def __parse_gtf(gtf_file):
        regions = {}
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip().split('\t')


                if parts[2] == 'exon':
                    chrom = parts[0]
                    # if chrom not in regions:
                    #     regions[chrom] = []
                    attrs = parts[8].split(";")
                    for a in attrs:
                        if 'transcript_id' in a:
                            name = a.split(" ")[2].replace("\"", "")
                            print(name)
                            start = int(parts[3])
                            end = int(parts[4])
                            if name not in regions:
                                regions[name] = {'chrom': chrom, 'positions': [], 'coverage': []}
                            for i in range(start, end):
                                if i not in regions[name]["positions"]:
                                    regions[name]['positions'].append(i)
                    # regions[chrom].append((start, end))
            for smorf in regions:
                # print(smorf)
                if len((regions[smorf]["positions"])) > 450:
                    print(len(regions[smorf]["positions"]))
        return regions

    def __index_bam_files(self):
        files = os.listdir(self.riboSeqSortedBamDir)
        for file in files:
            if file.endswith(".bam"):
                if not os.path.exists(f'{self.riboSeqSortedBamDir}/{file}.bai'):
                    cmd = f'samtools index {self.riboSeqSortedBamDir}/{file}'
                    os.system(cmd)

    def get_genome_coverage(self):
        print("Assessing genome coverage")
        # Initialize dictionary to store read coverage data
        self.genomeCoverage = {}
        self.check_dirs([self.riboSeqCovDir, self.riboSeqSortedBamDir])
        self.__sort_bam_files()
        self.__index_bam_files()
        bam_files = os.listdir(self.riboSeqSortedBamDir)

        # get coverage by position
        for file in bam_files:
            if file.endswith(".bam"):
                for smorf in self.regions:
                    # for region in self.regions[chrom]:
                    splat = smorf.split(":")
                    first_part = splat[0]
                    if '+' in first_part:
                        chrom = first_part.split("+")[1]
                    else:
                        chrom = first_part.split("-")[1]
                    coords = splat[1].split("_")[0].split("-")
                    start, end = coords[0], coords[1]
                    # print(self.regions[chrom])
                    # start, end = region[0], region[1]
                    locus = f'{chrom}:{start}-{end}'
                    cmd = (f'{self.toolPaths["samtools"]} depth {self.riboSeqSortedBamDir}/{file} -r {locus} >> '
                           f'{self.riboSeqCovDir}/{file}_coverage.txt')
                    os.system(cmd)

    # Filter coverage data to include only positions within specified regions
    @staticmethod
    def __filter_coverage_data(positions, coverage, regions):
        filtered_positions = []
        filtered_coverage = []
        for pos, cov in zip(positions, coverage):
            for region_start, region_end in regions:
                if region_start <= pos <= region_end:
                    filtered_positions.append(pos)
                    filtered_coverage.append(cov)
                    break  # Break the loop if the position is within any of the regions
        return filtered_positions, filtered_coverage

    def __sort_bam_files(self):
        print("Sorting bam files")
        files = os.listdir(self.riboSeqAlnDir)
        for file in files:
            if file.endswith(".sam"):
                if not os.path.exists(f'{self.riboSeqSortedBamDir}/{file[:-4]}_sorted.bam'):
                    cmd = (f'{self.toolPaths["samtools"]} sort {self.riboSeqAlnDir}/{file} -@ {self.args.threads} -o '
                           f'{self.riboSeqSortedBamDir}/{file[:-4]}_sorted.bam')
                    os.system(cmd)

    def get_orf_coverage(self):
        print("Assessing ORF coverage.")
        def read_coverage_file(coverage_file):
            coverage = {}
            with open(coverage_file, 'r') as file:
                for line in file:
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    if chrom not in coverage:
                        coverage[chrom] = {}
                    coverage[chrom][int(fields[1])] = int(fields[2])
            return coverage

        files = os.listdir(self.riboSeqCovDir)
        coverage_files = [file for file in files if file.endswith("_coverage.txt")]
        # Iterate through each coverage file
        for coverage_file in coverage_files:
            print(f"Assessing coverage for {coverage_file}.")
            # Read coverage data from the file

            coverage = read_coverage_file(f'{self.riboSeqCovDir}/{coverage_file}')
            for smorf in self.regions:
                chrom = self.regions[smorf]['chrom']
                for i, coord in enumerate(self.regions[smorf]['positions']):
                    if chrom in coverage:
                        if int(coord) in coverage[chrom]:
                            # for pos, cov in zip(coverage[chrom]['position'], coverage[chrom]['coverage']):
                            self.regions[smorf]['coverage'].append(coverage[chrom][coord])
                            if i not in self.smorfCoverages:
                                self.smorfCoverages[i] = []
                            if smorf not in self.singleCoverage:
                                self.singleCoverage[smorf] = {}
                            if i not in self.singleCoverage[smorf]:
                                self.singleCoverage[smorf][i] = []
                            self.singleCoverage[smorf][i].append((coverage[chrom][coord]))

                            # normalizing the coverage by the smorf length, so we dont get our statistics skewed by
                            # smorfs with different sequence length

                            self.smorfCoverages[i].append((coverage[chrom][coord])/len(self.regions[smorf]["positions"]))


                # print(self.regions[smorf])
                # break
            # break
            # print(self.regions)
            # for smorf in self.regions:
            #     if len(self.regions[smorf]['coverage']) > 1:
            #         print(self.regions[smorf])
                # print(self.regions[smorf]['coverage'])

    def plot_coverage(self):
        data = {"positions": [], "coverage": []}
        # for smorf in self.regions:
        #     # print(len(self.regions[smorf]["positions"]))
        #     for i, pos in enumerate(self.regions[smorf]["positions"]):
        #         data["positions"].append(pos)
        #         data["coverage"].append(self.regions[smorf]["coverage"][i])
        # print(data)
        for pos in self.smorfCoverages:
            for cov in self.smorfCoverages[pos]:
                data['positions'].append(pos)
                data['coverage'].append(cov)
        ax = sns.barplot(data=data, x="positions", y="coverage", errwidth=0.1, color='#E57C7C')
        plt.xlabel("Nucleotide position")
        plt.ylabel("Read counts")
        # plt.xticks(rotation=180)
        # print(data["positions"])
        # print(len(data["positions"]))
        plt.xticks(range(0, 451, 50))
        plt.savefig(f'{self.covPlots}/overall_riboSeq_coverage_smorfs.png')
        # plt.show()
            # for i, pos in enumerate(self.regions[smorf]["positions"]):
                # print(len(self.regions[smorf]["positions"]))
                # print(i)
                # print(self.regions[smorf]["coverage"])
                # data[i].append(self.regions[smorf]['coverage'][i])
        # print(data)

    def plot_coverage_single_orfs(self):
        print("Generating Ribo-seq coverage plots for each identified smORF.")
        for smorf in tqdm(self.singleCoverage):
            plt.clf()
            data = {"positions": [], "coverage": []}

            for pos in self.singleCoverage[smorf]:
                for cov in self.singleCoverage[smorf][pos]:
                    data['positions'].append(pos)
                    data['coverage'].append(cov)
            ax = sns.barplot(data=data, x="positions", y="coverage", errwidth=0.2, color='#E57C7C', ci=None)
            plt.xlabel("Nucleotide position")
            plt.ylabel("Read counts")
            seq_len = len(self.singleCoverage[smorf])
            plt.xlim((-1, seq_len))
            plt.xticks(range(0, seq_len, 50))
            plt.savefig(f'{self.covPlots}/{smorf}_riboSeq_coverage.png')
