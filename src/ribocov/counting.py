import os
import sys

import pandas as pd
import numpy as np
from ..pipeline_config import PipelineStructure


class FeatureCounts(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.set_translation_attrs()

        self.referencePlusRescoredGTF = f'{self.countsDir}/rescored_smorfs_plus_reference_annotation.gtf'

    def append_reference_annotation(self):
        cmd = f'cat {self.rescoredMicroproteinsGTF} {self.args.gtf} > {self.referencePlusRescoredGTF}'
        os.system(cmd)

    def run_feature_counts(self):
        folder = self.args.aln
        if self.args.aln is None:
            folder = self.riboSeqAlnDir
        files = os.listdir(folder)
        input_files = ''
        for file in files:
            if self.args.aln is not None:
                if file.endswith(".sam") or file.endswith(".bam"):
                    input_files += f' {folder}/{file}'
            else:
                if file.endswith("Aligned.out.sam"):
                    input_files += f' {folder}/{file}'

        shortcut = self.toolPaths["featureCounts"]
        gtf = self.referencePlusRescoredGTF
        # for gtf in gtf_files:
        # print(input_files)
        print("feature counts")
        gtf_name = gtf.split("/")[-1][:-4]
        cmd = f'{shortcut} -T {self.args.threads} -a {gtf} -o ' \
              f'{self.rawCountsDir}/{gtf_name}_counts.txt{input_files}'
        os.system(cmd)
        print(cmd)
        cmd_mm = f'{shortcut} -T {self.args.threads} -a {gtf} -M -o ' \
                 f'{self.rawCountsDir}/{gtf_name}_multimappers_counts.txt{input_files}'
        os.system(cmd_mm)
        cmd_mm_ambig = f'{shortcut} -T {self.args.threads} -a {gtf} -M -O -o ' \
                       f'{self.rawCountsDir}/{gtf_name}_multimappers_ambiguous_counts.txt{input_files}'
        os.system(cmd_mm_ambig)
        cmd_ambig = f'{shortcut} -T {self.args.threads} -a {gtf} -O -o' \
                    f' {self.rawCountsDir}/{gtf_name}_ambiguous_counts.txt{input_files}'
        os.system(cmd_ambig)


class CountNormalizer:
    def __init__(self, feature_counts, groups=None):
        self.counts = feature_counts
        # if groups is not None:
        #     self.groups = pd.read_csv(groups, sep='\t')
        # else:
        #     self.groups = groups
        self.geneLengths = {}

        self.scalingFactor = {}

        self.RPKMs = {}

    def get_rpkm(self):
        """ counts all reads and divide by 1.000.000 """

        df = pd.read_csv(self.counts, sep='\t', header=1)

        # if self.groups is not None:
        #     groups = self.groups
        #     smorfs = groups["smorf"].tolist()
        #     df = df[df["Geneid"].isin(smorfs)]
        # print(df)
        genes = df["Geneid"].tolist()
        self.genes = genes
        lens = df["Length"].tolist()

        samples = df.columns[6:]
        counts = {}    # counts

        for sample in samples:
            counts[sample] = df[sample].tolist()    # creates a list for each sample containing the counts for each gene
            if sample not in self.RPKMs:
                self.RPKMs[sample] = []

        for sample in counts:    # calculates the scaling factor per 1.000.000
            total = np.sum(counts[sample])
            scaling_factor = total/1000000
            self.scalingFactor[sample] = scaling_factor

        for i, gene in enumerate(genes):
            gene_length = lens[i]/1000
            for sample in samples:
                count = counts[sample][i]

                rpm = count/self.scalingFactor[sample]
                rpkm = rpm/gene_length
                self.RPKMs[sample].append(rpkm)

    def save_table(self, output):
        data = {'Geneid': []}

        for gene in self.genes:
            data['Geneid'].append(gene)
        for sample in self.RPKMs:
            if sample not in data:
                data[sample] = []
            rpkms = self.RPKMs[sample]
            for rpkm in rpkms:
                data[sample].append(rpkm)

        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)




