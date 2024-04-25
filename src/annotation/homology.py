import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML

from ..stats import Stater, DunnWithTukey
from ..pipeline_config import PipelineStructure


class Hit:
    def __init__(self, entry, start, end, seq):
        self.entry = entry
        self.start = start
        self.end = end
        self.seq = seq


class Homologs(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        """

        :param fasta: fasta with spliced CDSs
        """
        self.outdir = self.args.outdir
        # self.mappingClassDir = f~
        # self.__check_dirs([self.mappingClassDir, self.mappingHomologyDir])

        self.genome = self.args.genome
        # self.__get_transcriptome_sequences()
        self.__filter_split_nuc()
        self.fasta = self.splitMicroproteins
        self.fastaToBlast = f'{self.mappingHomologyDir}/cds_to_blast.fasta'

        self.blastedXML = f'{self.mappingHomologyDir}/smorfs_blasted_to_genome.xml'

        self.__fix_fasta()

        self.homologs = {}

        # self.customPalette = ["#BF7CD5", "#deb887", "#88C783", "#E75151", "#6495ed"]
        self.customPalette = ["#38322C", '#516067', '#C38D4F', '#A74A43', '#D8CFC0']
        self.mappingGroupsHomologs = f'{self.mappingHomologyDir}/mapping_groups_homologs.txt'

    def __filter_split_nuc(self):
        if not os.path.exists(self.splitMicroproteins):
            fasta = self.select_fasta()
            filtered = SeqIO.parse(fasta, 'fasta')
            entries = []
            for record in filtered:
                entries.append(str(record.description))
            to_write = []
            dbs = os.listdir(self.translationDir)
            for db in dbs:
                files = os.listdir(f'{self.translationDir}/{db}')

                for file in files:
                    if file.endswith("split_nuc"):
                        records = SeqIO.parse(f'{self.translationDir}/{db}/{file}', 'fasta')
                        for record in records:
                            # print(str(record.description))

                            if str(record.description) in entries:
                                to_write.append(f'>{str(record.description)}\n{str(record.seq)}\n')
            with open(self.splitMicroproteins, 'w') as out:
                out.writelines(to_write)

    def __check_dirs(self, folders):
        if type(folders) == list:
            for folder in folders:
                if not os.path.exists(folder):
                    os.mkdir(folder)
        elif type(folders) == str:
            os.mkdir(folders)

    def __fix_fasta(self):
        if not os.path.exists(self.fastaToBlast):
            to_write = []
            records = SeqIO.parse(self.fasta, 'fasta')
            for record in records:
                seq = str(record.seq).upper()
                to_write.append(f'>{str(record.description)}\n{seq}\n')
            with open(self.fastaToBlast, 'w') as outfile:
                outfile.writelines(to_write)

    def __get_transcriptome_sequences(self):
        print(f"Generating transcriptome fasta file with nucleotide sequences.")
        dbs = os.listdir(self.translationDir)
        for db in dbs:
            files = os.listdir(f'{self.translationDir}/{db}')
            for file in files:
                if file.endswith("gtf") and not file.endswith("_ORFs.gtf"):
                    cmd = (f'gffread -w {self.transcriptomeFasta} -g {self.args.genome} '
                           f'{self.translationDir}/{db}/{file}')
                    os.system(cmd)

    def blast(self):
        if not self.args.alignToTranscriptome:
            subject = self.genome
        else:
            self.__get_transcriptome_sequences()
            subject = self.transcriptomeFasta
        print(f"Blasting novel smORFs against the transcriptome to identify paralogs.")
        cmd = (f'blastn -query {self.fastaToBlast} -subject {subject} -evalue 0.01 -outfmt 5 -out'
               f' {self.blastedXML} -task blastn-short -soft_masking false -dust no')
        os.system(cmd)

    def parse(self, evalue=0.001, score=50, pc_id=0, qcov=0):
        for record in NCBIXML.parse(open(self.blastedXML)):
            if record.alignments:  # skip queries with no matches
                query = record.query
                if query not in self.homologs:
                    self.homologs[query] = []
                for align in record.alignments:
                    chrom = align.hit_def
                    for hsp in align.hsps:
                        if (hsp.expect <= evalue and hsp.bits >= score and self.count_mismatches(hsp) <= self.args.maxMismatches):
                            hit = Hit(entry=chrom, seq=hsp.sbjct, start=hsp.sbjct_start, end=hsp.sbjct_end)
                            self.homologs[query].append(hit)

    @staticmethod
    def count_mismatches(hsp):
        mismatches = 0
        for q, s in zip(hsp.query, hsp.sbjct):
            if q != s:
                mismatches += 1
        return mismatches

    def save_clusters(self):
        df = pd.read_csv(self.microproteinsMappingGroupsExclusive, sep='\t')
        orfs = {}
        smorfs, clusters = df["smorf"].tolist(), df["group"].tolist()

        cluster_for_plot = {'group': [], 'homologs': []}

        for smorf, groups in zip(smorfs, clusters):
            # if cluster not in cluster_for_plot:
            #     cluster_for_plot[cluster] = []
            # clusters = groups.split(",")
            # for cluster in clusters:
            if smorf in self.homologs:
                cluster_for_plot['group'].append(groups)
                cluster_for_plot['homologs'].append(len(self.homologs[smorf]))
        df = pd.DataFrame(data=cluster_for_plot)
        df.to_csv(self.mappingGroupsHomologs, sep='\t', index=False)

    def plot(self):
        df = pd.read_csv(self.mappingGroupsHomologs, sep='\t')
        # print(df)
        # df['homologs'] = df['homologs'].apply(lambda x: np.log1p(x))
        data = {'group': df["group"].tolist(), "homologs": df["homologs"].tolist()}
        # sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True)
        sns.set(rc={'font.size': 14})
        print(data)
        # Set the style back to the default
        sns.set_style("white")
        sns.set_palette(palette=self.customPalette)
        # print(df)
        # ax = sns.boxplot(data=df, x='group', y='homologs', hue='group', width=0.8, dodge=False, showfliers=False)
        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
        plt.legend().remove()
        dictio = {}
        order = []
        for group, intensity in zip(data['group'], data['homologs']):
            if group not in order:
                order.append(group)
            if group not in dictio:
                dictio[group] = []
            dictio[group].append(intensity)

        dunn = DunnWithTukey(data_frame=df)
        dunn.dunn_with_fdr_bh(val_col="homologs", group_col="group")
        values = df["homologs"].tolist()
        # for value in values:
        #     if value == 0:
        #         print(value)
        # dunn.posthoc_tukey()
        # print(dunn.result_posthoc)
        dunn.plot_with_pvalues(order=['Default', 'MM', 'Amb', 'MM_Amb', 'No coverage'], ylabel='Homologs',
                               xlabel='Groups')