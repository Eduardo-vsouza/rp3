import os
import sys
import re

from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
import seaborn as sns
import scipy
import matplotlib.pyplot as plt

from ..pipeline_config import PipelineStructure


class ProteinCoverage(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.protCovDir = f'{self.outdir}/protein_coverage'
        self.check_dirs([self.protCovDir])
        self.microproteins = {}

    def get_microprotein_sequences(self):
        records = SeqIO.parse(self.rescoredMicroproteinsFasta, 'fasta')
        for record in records:
            entry = str(record.description)
            seq = str(record.seq)
            self.microproteins[entry] = {'seq': seq, 'peptides_anno': [], 'peptides': [], 'utps': []}

    def get_mass_spec_peptides(self):
        groups = os.listdir(self.rescorePostProcessDir)
        for group in groups:
            df = pd.read_csv(f'{self.rescorePostProcessDir}/{group}/peptides_fixed.txt', sep='\t')
            df = df[df["q-value"] <= self.args.qvalue]
            proteins = df["proteinIds"].tolist()
            peptides = df["peptide"].tolist()
            for protein, peptide in zip(proteins, peptides):
                letters_only = re.findall('[a-zA-Z]', peptide)

                # Join the list of letters to form a single string
                pep = ''.join(letters_only)
                anno = False
                utp = False
                protein_list = protein.split(",")

                if 'ANNO' in protein:
                    anno = True
                if ',' not in protein:
                    utp = True
                for prot in protein_list:
                    if prot in self.microproteins:
                        if '_F:' in prot:
                            if anno:
                                self.microproteins[prot]['peptides_anno'].append(pep)
                            else:
                                if utp:
                                    self.microproteins[prot]['utps'].append(pep)
                                else:
                                    self.microproteins[prot]['peptides'].append(pep)

    def plot(self):

        def add_coverage(col, data):
            subsets = {'peptides': 'MP_nonUTP', 'peptides_anno': 'Annotated_nonUTP', 'utps': 'UTP'}
            for pep in self.microproteins[mp][col]:
                index = seq.find(pep)
                for j in range(index, index + len(pep)):
                    data['pos'].append(j)
                    data['coverage'].append(1)
                    data['subset'].append(subsets[col])
                    pos = j
                    if pos not in positions:
                        positions[pos] = {}
                        subset = subsets[col]
                        if subset not in positions[pos]:
                            positions[pos][subset] = 0
                        positions[pos][subset] += 1
            return data

        for mp in self.microproteins:
            positions = {}

            plt.clf()
            data = {'pos': [], 'coverage': [], 'subset': []}
            seq = self.microproteins[mp]['seq']
            # for i, aa in enumerate(seq):
            data = add_coverage(col='peptides', data=data)
            data = add_coverage(col='peptides_anno', data=data)
            data = add_coverage(col='utps', data=data)
            print(data)
            pep_data = {}
            # for i, pos in enumerate(data["positions"]):
            #     if
            df = pd.DataFrame(data=data)
            # unique positions
            udata = {'pos': [], 'coverage': [], 'subset': []}
            for pos in positions:
                for col in positions[pos]:
                    udata['pos'].append(pos)

                    udata['coverage'].append(positions[pos][col])
                    udata['subset'].append(col)
            udf = pd.DataFrame(data=udata)
            df_wide = udf.pivot(index="pos", columns="subset", values="coverage")
            for i, aa in enumerate(seq):
                if i not in udata['pos']:
                    udata['pos'].append(i)
                    udata['coverage'].append(0)
                    udata['subset'].append('UTP')
            # sns.set_style(style="white")
            # sns.set_palette(palette="deep")
            # _, ax = plt.subplots()
            # sns.lineplot(x="pos", y="coverage", hue="subset", data=udf, ax=ax, legend=False)
            # ax.set_prop_cycle(None)
            # df_wide.plot.area(stacked=False, alpha=0.2, ax=ax)
            x = udf["pos"].tolist()
            y = udf["coverage"].tolist()
            data = udata
            # Iterate over subsets
            sorted_data = sorted(zip(data['pos'], data['coverage'], data['subset']))

            # Unzip sorted data
            sorted_pos, sorted_coverage, sorted_subset = zip(*sorted_data)
            import numpy as np
            # Define colors for the subsets
            subset_colors = {'Annotated_nonUTP': '#E3E194', 'UTP': '#8CADDE'}
            plt.figure(figsize=(20, 5))

            # Iterate over the subsets
            for subset in set(sorted_subset):  # Use set() to get unique subsets
                # Extract the subset data for the current subset
                subset_pos = []
                subset_data = []
                for pos, cov, sub in zip(sorted_pos, sorted_coverage, sorted_subset):
                    if sub == subset:
                        subset_pos.append(pos)
                        subset_data.append(cov)

                # Add an extra point at the end to ensure the line drops vertically to zero
                subset_pos.append(max(sorted_pos) + 1)
                subset_data.append(0)

                # Plot the subset line
                plt.step(subset_pos, subset_data, where='post', label=subset, color=subset_colors.get(subset, 'orange'))

                # Fill the area under the subset line
                plt.fill_between(subset_pos, subset_data, step='post', color=subset_colors.get(subset, 'orange'),
                                 alpha=0.5)
            # print(len(sorted_pos))
            # print(len(seq))
            # Show the plot
            plt.xlabel('AA position')
            plt.ylabel('Mass Spec peptides')
            # plt.title('Filled Line Plot')
            plt.legend()
            plt.grid(True)
            # plt.plot(x, y, marker='.', lw=1)
            # ax = plt.gca()
            # d = scipy.zeros(len(y))
            #
            # ax.fill_between(x, y, where=y >= d, interpolate=True, color='blue')
            # ax.fill_between(x, y, where=y <= d, interpolate=True, color='red')

            # ax = sns.kdeplot(data=data, x="pos", hue="subset", fill=True, common_norm=False, palette="crest",
            # alpha=.5, linewidth=0,)
            # plt.stackplot(x=data["pos"], y=data["coverage"], alpha=0.5)
            # print(width)
            # ax = sns.barplot(data=data,  x="pos", y="coverage", hue="subset", estimator=sum, dodge=False)
            # plt.xticks(sorted_pos, [aa for aa in seq])
            plt.show()
            # plt.savefig(f'{self.protCovDir}/{mp}_peptide_coverage.png')
            break







