import os
import re

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from ..stats import DunnWithTukey
from ..pipeline_config import PipelineStructure


class Repeater(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.repeatsFile = self.args.repeatsFile
        self.repeatsDir = f'{self.mappingClassDir}/repeats'
        self.check_dirs([self.repeatsDir])
        self.fixedRepeatsFile = f'{self.repeatsDir}/fixed_repeats_file.txt'
        # self.fixedRepeatsFile = self.args.repeatsFile
        os.system(f'cp {self.args.repeatsFile} {self.fixedRepeatsFile}')
        self.smORFsRepeatsFile = f'{self.repeatsDir}/smorfs_in_repeats.txt'
        # self.__check_dirs([self.outdir])
        self.repeats = {}
        self.smORFsRepeats = {}

        self.customPalette = ["#BF7CD5", "#deb887", "#88C783", "#E75151", "#6495ed"]
        self.customPalette = ["#82AC7C", '#6d597A', '#C38D4F', '#A74A43', '#9DA4D4']


    @staticmethod
    def __check_dirs(folders):
        for folder in folders:
            if not os.path.exists(folder):
                os.mkdir(folder)


    def fix_repeats_file(self):
        fixed = []
        # Read the file and split the lines
        with open(self.repeatsFile, 'r') as f:
            lines = f.readlines()

        # Define the regular expression pattern to split the columns


        # Loop through the lines and split the columns using the regular expression pattern
        for i, line in enumerate(lines):
                splat = line.split(" ")
                line = '\t'.join([el for el in splat if el != ''])
                # Remove any trailing whitespace and newline characters
                # line = line.strip()
                # Split the columns using the regular expression pattern
                # columns = re.split(pattern, line)
                # Add the split columns to the list
                # fixed.append('\t'.join(columns))
                fixed.append(line)
        with open(self.fixedRepeatsFile, 'w') as outfile:
            outfile.writelines(fixed)

    def get_repeats(self):
        df = pd.read_csv(self.fixedRepeatsFile, sep='\t')
        r_starts, r_ends = df["query_start"].tolist(), df["query_end"].tolist()
        # print(df)
        # print(df["repeat_start"])
        chroms = df["query_seq"].tolist()
        repeats = df["repeat_class_family"].tolist()
        for i, r_start in enumerate(r_starts):
            r_end = r_ends[i]
            chrom = chroms[i]
            family = repeats[i]
            if chrom not in self.repeats:
                self.repeats[chrom] = {'starts': [], 'ends': [], 'repeats': []}
            self.repeats[chrom]['starts'].append(r_start)
            self.repeats[chrom]['ends'].append(r_end)
            self.repeats[chrom]['repeats'].append(family)

    def check_smorfs_within_repeats(self):
        smorf_df = pd.read_csv(self.microproteinsMappingGroupsExclusive, sep='\t')
        print(smorf_df)
        smorfs, groups = smorf_df["smorf"].tolist(), smorf_df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            start, end, chrom = self.__get_smorf_info(smorf)
            for i, r_start in enumerate(self.repeats[chrom]['starts']):
                r_start = int(r_start)
                r_end = int(self.repeats[chrom]['ends'][i])
                if start in range(r_start, r_end) or end in range(r_start, r_end):
                    if smorf not in self.smORFsRepeats:
                        self.smORFsRepeats[smorf] = []
                    self.smORFsRepeats[smorf].append(self.repeats[chrom]['repeats'][i])
        data = {'smorfs': [], 'repeats': []}
        for smorf in self.smORFsRepeats:
            data['smorfs'].append(smorf)
            data['repeats'].append(self.smORFsRepeats[smorf])
        df = pd.DataFrame(data=data)
        df.to_csv(self.smORFsRepeatsFile, sep='\t', index=False)

    @staticmethod
    def __get_smorf_info(smorf):
        splat = smorf.split("_F:")
        first_half = splat[0].split(":")
        coords = first_half[1].split("-")
        start, end = coords[0], coords[1]
        chrom_part = first_half[0]
        if '+' in chrom_part:
            chrom = chrom_part.split("+")[1]
        else:
            chrom = chrom_part.split("-")[1]
        return int(start), int(end), chrom

    def plot_smorfs_in_repeats_based_on_clusters(self):
        repeats_df = pd.read_csv(self.smORFsRepeatsFile, sep='\t')
        if self.args.exclusiveMappingGroups:
            cluster_df = pd.read_csv(self.microproteinsMappingGroupsExclusive, sep='\t')
        else:
            cluster_df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        print(cluster_df)
        rename_dict = {'mm': 'MM', 'mm_amb': 'MM_amb', 'amb': 'Amb', 'default': 'Default', 'No coverage': 'No coverage'}
        cluster_df['group'].replace(rename_dict, inplace=True)
        data = {'smorf': [], 'repeat': [], 'group': []}
        clusters = {}
        totals = {}


        orfs_cl, groups = cluster_df["smorf"].tolist(), cluster_df["group"].tolist()
        # if self.args.exclusiveMappingGroups:
        for orf, group in zip(orfs_cl, groups):
            clusters[orf] = group
            if group not in totals:
                totals[group] = 0
            totals[group] += 1
        print(totals)
        # else:
        #     for orf, grouplist in zip(orfs_cl, groups):
        #         splat = grouplist.split(",")
        #         for group in splat:

        repeat_prop = {'group': [], 'repeats': []}
        repeats = {}
        smorfs_r, repeat_types = repeats_df["smorfs"].tolist(), repeats_df["repeats"].tolist()
        for smorf, repeat in zip(smorfs_r, repeat_types):
            # for r in repeat:
            if not self.args.exclusiveMappingGroups:
                # for grouplist in clusters[smorf]:
                grouplist= clusters[smorf]
                splat = grouplist.split(",")
                for group in splat:
                    if group not in repeats:
                        repeats[group] = 0
                    repeats[group] += 1
            else:
                if clusters[smorf] not in repeats:
                    repeats[clusters[smorf]] = 0
                    # totals[clusters[smorf]] = 0
                repeats[clusters[smorf]] += 1
                # totals[clusters[smorf]] += 1

        for group in repeats:
            repeat_prop['group'].append(group)
            # repeat_prop['repeats'].append(repeats[group]/totals[group])
            repeat_prop['repeats'].append(repeats[group])

            # repeat_prop['group'].append(clusters[smorf])
            # data['repeat'].append(r)

        sns.set(rc={'font.size': 14})

        # Set the style back to the default
        sns.set_style("white")
        sns.set_palette(palette=self.customPalette)
        sns.set_palette(palette='coolwarm_r')
        # colors = ["#D9BFA1","#304D6D"]

        df = pd.DataFrame(data=repeat_prop)
        # ax = sns.barplot(data=df, x='group', y='repeats', edgecolor='black', order=['Default', 'MM', 'Amb', 'MM_Amb',
        #                                                                             'No coverage'])
        ax = sns.barplot(data=df, x='group', y='repeats', edgecolor='black', order=['Default', 'MM', 'Amb', 'MM_Amb',
                                                                                    'No coverage'])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xticks(rotation=45, ha="center")
        plt.ylabel("smORFs in repeat regions")
        # plt.legend()
        plt.show()
