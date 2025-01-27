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
        run = self.verify_checkpoint(outfile=self.referencePlusRescoredGTF,
                                     step=f"merging of Rp3 and reference GTF files")
        if run:
            if self.args.externalGTF:
                print(f"--Merging external GTF and reference annotations")
                cmd = f'cat {self.args.externalGTF} {self.args.gtf} > {self.referencePlusRescoredGTF}'
                os.system(cmd)
            else:
                if os.path.exists(self.rescoredMicroproteinsFasta):
                    cmd = f'cat {self.rescoredMicroproteinsGTF} {self.args.gtf} > {self.referencePlusRescoredGTF}'
                else:
                    if os.path.exists(self.microproteinsBlast):
                        cmd = f'cat {self.microproteinsBlast} {self.args.gtf} > {self.referencePlusRescoredGTF}'
                    else:
                        # cmd = f'cat {self.args.gtf} {self.mergedResults}/merged_predicted_microproteins.gtf > {self.referencePlusRescoredGTF}'
                        cmd = f'cp {self.uniqueMicroproteinsGTF} {self.referencePlusRescoredGTF}'

            os.system(cmd)

    def run_feature_counts(self):
        folder = self.args.aln
        if self.args.aln is None:
            folder = self.riboSeqAlnDir
        files = os.listdir(folder)
        input_files = ''
        for file in files:
            if self.args.aln is not None:
                if file.endswith("Aligned.out.sam") or file.endswith(".bam"):
                    input_files += f' {folder}/{file}'
            else:
                if file.endswith("Aligned.out.sam"):
                    input_files += f' {folder}/{file}'

        shortcut = self.toolPaths["featureCounts"]
        # print(input_files)
        gtf = self.referencePlusRescoredGTF
        # for gtf in gtf_files:
        # print(input_files)
        # print("feature counts")
        gtf_name = gtf.split("/")[-1][:-4]
        threads = self.args.threads
        if self.args.threads > 64:
            threads = 64
        attribute = 'gene_id'
        feature = 'CDS'
        out_ambig = f'{self.rawCountsDir}/{gtf_name}_ambiguous_counts.txt'
        run = self.verify_checkpoint(outfile=out_ambig, step="read counting with featureCounts")
        if run:
            cmd = f'{shortcut} -T {threads} -t {feature} -a {gtf} -g {attribute} -o ' \
                  f'{self.rawCountsDir}/{gtf_name}_counts.txt{input_files}'
            os.system(cmd)
            # print(cmd)
            cmd_mm = f'{shortcut} -T {threads} -t {feature} -a {gtf} -M -g {attribute} -o ' \
                     f'{self.rawCountsDir}/{gtf_name}_multimappers_counts.txt{input_files}'
            os.system(cmd_mm)
            cmd_mm_ambig = f'{shortcut} -T {threads} -t {feature} -a {gtf} -M -O -g {attribute} -o ' \
                           f'{self.rawCountsDir}/{gtf_name}_multimappers_ambiguous_counts.txt{input_files}'
            os.system(cmd_mm_ambig)
            cmd_ambig = f'{shortcut} -T {threads} -t {feature} -a {gtf} -O -g {attribute} -o' \
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

    def get_rpkm(self, min_raw_counts):
        """ counts all reads and divide by 1.000.000 """

        df = pd.read_csv(self.counts, sep='\t', header=1)
        # df = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
        threshold = min_raw_counts

        # Iterate through each column (sample) and replace values below the threshold with 0
        for column in df.columns[6:]:  # Assuming the first column contains gene names
            df[column] = df[column].apply(lambda x: 0 if x < threshold else x)
        # df = df[df.iloc[:, 6:].apply(lambda x: (x > 10).any(), axis=1)]

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




