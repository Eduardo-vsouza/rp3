import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from ..stats import Stater, DunnWithTukey
from ..pipeline_config import PipelineStructure


class Isoforms(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.isoformsDir = f'{self.mappingClassDir}/isoforms'
        self.check_dirs([self.isoformsDir])
        self.filteredIsoformsGTF = f'{self.isoformsDir}/filtered_isoforms.gtf'
        self.overlappedGTF = f'{self.isoformsDir}/overlapped.gtf'
        self.overlaps = {}
        self.customPalette = ["#BF7CD5", "#deb887", "#88C783", "#E75151", "#6495ed"]
        self.customPalette = ["#82AC7C", '#6d597A', '#C38D4F', '#A74A43', '#9DA4D4']

    def filter_gtf_features(self):
        cmd = f"grep 'transcript	' {self.args.refGTF} > {self.filteredIsoformsGTF}"
        os.system(cmd)

    def intersect(self):
        if os.path.exists(self.rescoredMicroproteinsGTF):
            smorfs_gtf = self.rescoredMicroproteinsGTF
        else:
            smorfs_gtf = self.uniqueMicroproteinsGTF
        cmd = f'bedtools intersect -c -a {smorfs_gtf} -b {self.filteredIsoformsGTF} > {self.overlappedGTF}'
        print(cmd)
        os.system(cmd)

    def read_intersected(self):
        with open(self.overlappedGTF, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                line = line.rstrip()
                cols = line.split("\t")
                overlaps = cols[-1]
                gene = cols[8].split(";")[0].replace("gene_id ", "").replace("\"", "")
                self.overlaps[gene] = overlaps

    def plot_cluster_overlaps(self):
        overlaps = {'smorf': [], 'Overlapping features': [], 'group': []}
        if self.args.exclusiveMappingGroups:
            cluster_df = self.microproteinsMappingGroupsExclusive
        else:
            cluster_df = self.microproteinMappingGroupsForPlotsUnion
        df = pd.read_csv(cluster_df, sep='\t')
        rename_dict = {'mm': 'MM', 'mm_amb': 'MM_amb', 'amb': 'Amb', 'default': 'No coverage'}
        df['group'].replace(rename_dict, inplace=True)
        print(df)
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            if self.args.exclusiveMappingGroups:
                overlaps['smorf'].append(smorf)
                overlaps['group'].append(group)
                overlaps['Overlapping features'].append(int(self.overlaps[smorf])-1)
            else:
                splat = group.split(",")
                for g in splat:
                    overlaps['smorf'].append(smorf)
                    overlaps['group'].append(g)
                    overlaps['Overlapping features'].append(int(self.overlaps[smorf])-1)

        dictio = {}
        order = []
        for group, intensity in zip(overlaps['group'], overlaps['Overlapping features']):
            if group not in order:
                order.append(group)
            if group not in dictio:
                dictio[group] = []
            dictio[group].append(intensity)   # hi buddy
        df = pd.DataFrame(data=overlaps)    # why are you delving so deep into my code
        # sns.set(rc={'font.size': 12})    # do you need a friend

        sns.set(rc={'font.size': 14})

        # Set the style back to the default
        sns.set_style("white")
        sns.set_palette(palette=self.customPalette)
        sns.set_palette("coolwarm_r")
        dunn = DunnWithTukey(data_frame=df)
        dunn.dunn_with_fdr_bh(val_col="Overlapping features", group_col="group")
        # values = df["homologs"].tolist()
        # for value in values:
        #     if value == 0:
        #         print(value)
        # dunn.posthoc_tukey()
        # print(dunn.result_posthoc)
        dunn.plot_with_pvalues(order=['Default', 'MM', 'Amb', 'MM_Amb', 'No coverage'], ylabel='Overlapping features',
                               xlabel='Groups')
        # ax = sns.boxplot(data=df, x='group', y='Overlapping features', hue='group', width=0.8, dodge=False, showfliers=False)
        # stater = Stater(ax=ax, df=df, groups_x=dictio, x='', xlabel='group', ylabel='Overlapping features')
        # stater.test(order, loc='outside')
        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
        # plt.legend().remove()
        # plt.legend.remove()
        # plt.show()
