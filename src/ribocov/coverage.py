import os
import sys
import collections

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from nheatmap import nhm

from .counting import CountNormalizer
from ..pipeline_config import PipelineStructure


class RiboSeqCoverage(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.args = args
        # self.groups = args.groups

        self.counts = {}

    def check_dir(self):
        dirs = [self.rpkmDir, self.plotsDir]
        for folder in dirs:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def feature_counts_to_rpkm(self):
        counts = os.listdir(self.rawCountsDir)
        for count in counts:
            if count.endswith(".txt"):
                data = CountNormalizer(feature_counts=f'{self.rawCountsDir}/{count}')
                data.get_rpkm()
                data.save_table(f'{self.rpkmDir}/{count.replace(".txt", "_rpkm.txt")}')

    def plot(self, cutoff=1):
        files = os.listdir(self.rpkmDir)
        for file in files:
            counts = {}
            df = pd.read_csv(f'{self.rpkmDir}/{file}', sep='\t')
            groups = df.columns[6:]
            group_rpkms = {group: df[group].tolist() for group in groups}
            genes = df["Geneid"].tolist()
            for i, gene in enumerate(genes):
                if '_F:' in gene:
                    if gene not in counts:
                        counts[gene] = []
                    for group in group_rpkms:
                        rpkm = group_rpkms[group][i]
                        counts[gene].append(rpkm)
            counts_mean = []
            genes = []
            if file not in self.counts:
                self.counts[file] = []
            for gene in counts:
                median = np.mean(counts[gene])
                if np.log1p(median) >= cutoff:
                    self.counts[file].append(median)
                # if median >= 100:
                #     median = 100
                counts_mean.append(np.log1p(median))
                genes.append(gene)
            plt.scatter(genes, counts_mean, s=5)
            plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off

            plt.xlabel("Proteogenomics smORFs")
            plt.ylabel("RPKM")
            plt.title(file)
            plt.savefig(f"{self.plotsDir}/{file}.png", dpi=300)
            plt.clf()
            # plt.show()

    def plot_heatmaps(self, cutoff=1):
        files = os.listdir(self.rpkmDir)
        rpkms = {'Gene': []}
        for file in files:
            df = pd.read_csv(f'{self.rpkmDir}/{file}', sep='\t')
            if file.endswith(".txt"):
                file = file.split("/")[-1]
                if 'multimappers' in file:
                    if 'ambiguous' in file:
                        file = 'MM_Amb'
                    else:
                        file = 'MM'
                else:
                    if 'ambiguous' in file:
                        file = 'Amb'
                    else:
                        file = 'Default'

                counts = {}
                if file not in rpkms:
                    rpkms[file] = []
                genes = df["Geneid"].tolist()
                groups = df.columns[6:]
                #
                # file_rpkms = []
                #
                # for group in groups:
                #     ndf = df[group].tolist()

                group_rpkms = {group: df[group].tolist() for group in groups}


                for i, gene in enumerate(genes):
                    if '_F:' in gene or gene.startswith("CUFF."):
                        if gene not in rpkms['Gene']:
                            rpkms['Gene'].append(gene)
                        if gene not in counts:
                            counts[gene] = []
                        for group in group_rpkms:
                            rpkm = group_rpkms[group][i]
                            counts[gene].append(float("{:.3f}".format(rpkm)))

                for gene in rpkms['Gene']:
                    mean = np.median(counts[gene])
                    # if mean < self.args.rpkm:
                        # mean = "Remove"
                        # rpkms[]
                    # print(mean)
                    # if mean == 0:
                    #     mean = 0.000001
                    # else:
                    # rpkms[file].append(np.log1p(mean))
                    rpkms[file].append(mean)


        df = pd.DataFrame(data=rpkms)
        # cols = df.columns
        # default = df[(df < 1).all(axis=1)]

        # for gene in genes:
        #     i += 1
        #     renamed.append(i)
            # splat = gene.split("chr")[1]
            # reduced = splat.split("_F:")[0]
            # renamed.append(reduced)
        # df = df.drop(columns=['Gene'])
        # df.insert(0, "Gene", renamed)
        sys.setrecursionlimit(10000)

        df = df.set_index(df.columns[0])
        df.to_csv(self.mappingGroupsRPKMs, sep='\t')
        df = df[df.gt(self.args.rpkm).any(axis=1)]
        for col in df.columns:
            if col != "Gene":
                df[col] = np.log1p(df[col])

        # default = df[(df < 1).all(axis=1)]
        # default.to_csv(f'{self.outdir}/heatmap_lower_than_1rpkm.txt',
        #                sep='\t')
        # df = df[df.gt(1).any(axis=1)]
        df.to_csv(self.mappingGroupsRPKMsFiltered, sep='\t')
        cmaps = {'gene cluster': 'inferno'}
        cmap_colors = [(0, '#7094EC'), (1, 'red')]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_colors)
        g = nhm(data=df, cmapCenter=cmap, showyticks=False)
        g.hcluster(optimal_ordering=True)

        df_ordered = g.data.iloc[g.rorder, g.corder]
        df_ordered.to_csv(f'{self.args.outdir}/ordered_heatmap.xls')
        print(df_ordered)
        # print(g.dendrogram)
        fig, plots = g.run()
        plt.savefig(f"{self.countsDir}/heatmap.png")
        plt.savefig(f"{self.countsDir}/heatmap.pdf")

        # plt.show()

        # visuz.gene_exp.hmap(df=df, dim=(3, 6), tickfont=(6, 4), zscore=0)
        # sns.clustermap(df)
        # plt.show()
        # print(df)
        # plt.show()
        # sns.clustermap(correlations, method="complete", cmap='RdBu', annot=True,
        #                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15, 12));

    def plot_bars(self):
        # checks number of orfs per file with rpkm higher than the cutoff
        files = []
        orfs = []
        for file in self.counts:
            name = '_'.join(file.split("_")[1:3])
            splat = file.split("_")
            print(splat)
            name += f'_{file.split("_")[7]}'
            files.append(name.replace("counts", "default"))
            orfs.append(len(self.counts[file]))
        plt.barh(files, orfs, edgecolor='black')
        plt.ylabel("Assembly")
        plt.xlabel("Number of covered smORFs")
        plt.show()


class CoverageClassification(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.df = pd.read_csv(self.mappingGroupsRPKMs, sep='\t')
        self.df = self.df[self.df["Gene"] != "Gene"]

        self.classification = {}

    def classify(self):
        genes = self.df["Gene"].tolist()
        cols = self.df.columns[1:]
        for i, gene in enumerate(genes):
            higher = 0
            classification = "No coverage"
            for col in cols:
                rpkms = self.df[col].tolist()
                # if gene not in self.classification:
                #     self.classification[gene] = {}
                rpkm = float(rpkms[i])
                if rpkm > 1 and col == "Default":
                    higher = rpkm
                    classification = col
                if rpkm > higher and classification != 'Default':
                    higher = rpkm
                    classification = col
            self.classification[gene] = classification
        print(self.classification)
        gene_counts = collections.Counter(self.classification.values())
        print(gene_counts)

    def save(self):
        data = {'smorf': [], 'group': []}
        for gene in self.classification:
            data['smorf'].append(gene)
            data['group'].append(self.classification[gene])
        df = pd.DataFrame(data=data)
        df.to_csv(self.mappingGroups, sep='\t', index=False)

