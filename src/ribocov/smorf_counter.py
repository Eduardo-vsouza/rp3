import os
import sys
import collections

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..pipeline_config import PipelineStructure
from ..stats import DunnWithTukey


class ORFCounter(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.groups = pd.read_csv(self.mappingGroups, sep='\t')

        self.ids = {}

        self.customPalette = ["#BF7CD5", "#deb887", "#88C783", "#E75151", "#6495ed"]

    def __filter_smorfs(self, fasta):
        records = SeqIO.parse(fasta, 'fasta')
        entries = []
        for record in records:
            entries.append(str(record.description))
        return entries


    def count_smorfs(self, filter_smorfs=None):
        if filter_smorfs is not None:
            entries = self.__filter_smorfs(fasta=filter_smorfs)
        smorfs = self.groups["smorf"].tolist()
        groups = self.groups["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            # # print(smorf)
            # # if group not in self.ids:
            # #     self.ids[group] = []
            # add = False
            # if filter_smorfs is not None:
            #     if smorf in entries:
            #         add = True
            #     else:
            #         add = False
            # else:
            #     add = True
            # if add:
                self.ids[smorf] = group
            # if smorf not in self.ids[group]:
            #     self.ids[group].append(smorf)
        print(self.ids)
        groups = collections.Counter(self.ids.values())
        # print(gene_counts)
        print(groups)
        df = pd.DataFrame.from_dict(groups, orient='index', columns=['Count'])
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'Group'}, inplace=True)
        sns.set(rc={'font.size': 14})

        # Set the style back to the default
        sns.set_style("white")
        # Plotting with Seaborn
        # sns.set(style='whitegrid')
        # plt.figure(figsize=(10, 6))
        sns.set_palette(palette=self.customPalette)
        sns.set_palette(palette='coolwarm_r')
        # sns.diverging_palette(220, 20, as_cmap=True)


        print(df)
        ax = sns.barplot(x='Group', y='Count', data=df, width=0.6, dodge=False, edgecolor='black',
                         order=['Default', 'MM', 'Amb', 'MM_Amb', 'No coverage'])
        # ax = sns.barplot(x='Group', y='Count', data=df, width=0.6, dodge=False, edgecolor='black',
        #                  order=['Amb', 'MM_Amb', 'No coverage'])
        # Rotate x-axis labels for better readability
        ax.set_xticklabels(ax.get_xticklabels())
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()


        # Add numbers on top of the bars
        for p in ax.patches:
            ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', xytext=(0, 5), textcoords='offset points')

        # Customize the plot
        ax.set_title('Counts by Group')
        # ax.set_xlabel('Group')
        ax.set_ylabel('Number of smORFs')
        # fig = plt.figure()
        plt.savefig(f'{self.plotsDir}/smorf_counts_per_mapping_group.pdf', dpi=600)
        plt.savefig(f'{self.plotsDir}/smorf_counts_per_mapping_group.png', dpi=600)

        plt.show()

        # ax = sns.barplot(data=self.groups, x='group', y='smorf')
