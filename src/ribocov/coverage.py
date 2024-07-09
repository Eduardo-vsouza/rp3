import os
import sys
import collections

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from Bio import SeqIO
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
        self.classificationRawCounts = {}


    def check_dir(self):
        dirs = [self.rpkmDir, self.plotsDir]
        for folder in dirs:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def feature_counts_to_rpkm(self):
        counts = os.listdir(self.rawCountsDir)
        for count in counts:
            if count.endswith(".txt"):
                """ i must get a list of all smorfs with at least 10 read counts before iterating this"""
                data = CountNormalizer(feature_counts=f'{self.rawCountsDir}/{count}')
                data.get_rpkm(min_raw_counts=self.args.minRawCounts)
                data.save_table(f'{self.rpkmDir}/{count.replace(".txt", "_rpkm.txt")}')

    # def __assess_low_counts(self):
    #     counts = os.listdir(self.rawCountsDir)
    #     smorfs = []
    #     for count in counts:
    #         if count.endswith(".txt"):
    #             df = pd.read_csv(f'{self.rawCountsDir}/{count}', sep='\t', header=1)
    #             # df = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    #
    #             df = df[df.iloc[:, 6:].apply(lambda x: (x > self.args.minRawCounts).any(), axis=1)]


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

    def plot_histograms(self, folder, header, title, cutoff=1, log=True, median=False):
        plt.clf()
        print(f"Plotting the distribution of smORFs {title} across mapping groups.")
        files = os.listdir(folder)
        # plt.clf()

        smorfs_cov = {}
        for file in files:
            if file.endswith("txt"):
                df = pd.read_csv(f'{folder}/{file}', sep='\t', header=header)
                if 'Start' in df.columns:
                    to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']
                    df = df.drop(columns=to_drop)
                df = df[df["Geneid"].str.contains("_F:")]
                genes = df["Geneid"].tolist()
                mapping_group = self.__check_mapping_group(file)
                for col in df.columns:
                    if col != "Geneid":
                        counts = df[col].tolist()
                        for gene, count in zip(genes, counts):
                            if gene not in smorfs_cov:
                                smorfs_cov[gene] = {}
                            if mapping_group not in smorfs_cov[gene]:
                                smorfs_cov[gene][mapping_group] = []
                            smorfs_cov[gene][mapping_group].append(count)

        data = {'mapping_group': [], 'read counts': []}
        all_counts = []
        for smorf in smorfs_cov:
            for mapping_group in smorfs_cov[smorf]:
                if any(count > float(cutoff) for count in smorfs_cov[smorf][mapping_group]):
                    data['mapping_group'].append(mapping_group)

                    all_counts.append(smorfs_cov[smorf][mapping_group])
                    if log:
                        if median:
                            data['read counts'].append(np.log1p(np.median(smorfs_cov[smorf][mapping_group])))
                        else:
                            data['read counts'].append(np.log1p(np.mean(smorfs_cov[smorf][mapping_group])))

                    # else:

                    # for count in smorfs_cov[smorf][mapping_group]:
                    #     if count > 0:
                    #         data['read counts'].append(np.log1p(count))
        # print(len(data['mapping_group']))
        # print(len(data['read counts']))
                        # if median:
                        #     data['read counts'].append(np.median(smorfs_cov[smorf][mapping_group]))
                        # else:
                        #     data['read counts'].append(np.mean(smorfs_cov[smorf][mapping_group]))
        from matplotlib import transforms

        set_bins = 40
        if len(set(map(round, data['read counts']))) < set_bins:
            bins = set_bins - (len(set(map(round, data['read counts']))) - set_bins)
        else:
            bins = set_bins
        # print(smorfs_cov)
        # if log:
        sns.set_palette(palette='coolwarm_r')

        ax = sns.histplot(data=data, x="read counts", hue="mapping_group", bins=30, stat="count",
                          multiple="stack", hue_order=['Default', 'MM', 'Amb', 'MM_Amb'])
        # ax = sns.jointplot(data=data, x="read counts", hue="mapping_group")
        # ax = sns.histplot(data=data, x="read counts", hue="mapping_group", bins=30, kde=True)

        #
        # else:
        #     f, (ax_top, ax_bottom) = plt.subplots(ncols=1, nrows=2, sharex=False, gridspec_kw={'hspace': 0.10})
        #     sns.histplot(data=data, x="read counts", hue="mapping_group", bins=30, kde=True, multiple="stack",
        #                  ax=ax_top, discrete=True)
        #     sns.histplot(data=data, x="read counts", hue="mapping_group", bins=30, kde=True, multiple="stack",
        #                  ax=ax_bottom, discrete=True)
        #
        #     # sns.barplot(data=data, y="read counts", x="mapping_group", ax=ax_top)
        #     # sns.barplot(data=data, y="read counts", x="mapping_group", ax=ax_bottom)
        #     ax_top.set_ylim(bottom=10)  # those limits are fake
        #     # ax1 = plt.axes()
        #     # ax1.set_visible(False)
        #     ax_top.set_xlabel('')
        #     # ax_top.title("")
        #     ax_bottom.set_ylim(0, 10)
        #
        #     # ax_top.set_xlim(0, 10)
        #     # ax_bottom.set_xlim(0, 10)
        #     ax_top.set_xticks([])

        #     # ax_top.xticks(False)
        #     ax = ax_top
        #     d = .015  # how big to make the diagonal lines in axes coordinates
        #     # arguments to pass to plot, just so we don't keep repeating them
        #     kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        #     ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
        #     ax2 = ax_bottom
        #     kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        #     ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

            # remove one of the legend
            # ax_bottom.legend_.remove()
        plt.xlabel(f"log({title})")
        plt.ylabel("smORFs")
        plt.title(self.expCutoffsSuffix)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if self.args.manualPlot:
            plt.show()
        else:
            plt.savefig(f'{self.plotsDir}/{title}_distribution_{self.expCutoffsSuffix}.pdf', dpi=600)
            plt.savefig(f'{self.plotsDir}/{title}_distribution_{self.expCutoffsSuffix}.png', dpi=600)

        # plt.show()

    @staticmethod
    def __check_mapping_group(file):
        if 'multimappers' in file:
            if 'ambiguous' in file:
                group = 'MM_Amb'
            else:
                group = 'MM'
        else:
            if 'ambiguous' in file:
                group = 'Amb'
            else:
                group = 'Default'
        return group

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
                df = df[df["Geneid"].str.contains("_F:")]
                # print(df)
                # print("1")
                genes = df["Geneid"].tolist()
                groups = df.columns[6:]
                groups = df.columns[1:]
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
        df = df[df.ge(self.args.rpkm).any(axis=1)]
        for col in df.columns:
            if col != "Gene":
                df[col] = np.log1p(df[col])
        # print("2")
        # default = df[(df < 1).all(axis=1)]
        # default.to_csv(f'{self.outdir}/heatmap_lower_than_1rpkm.txt',
        #                sep='\t')
        # df = df[df.gt(1).any(axis=1)]
        df.to_csv(self.mappingGroupsRPKMsFiltered, sep='\t')
        cmaps = {'gene cluster': 'inferno'}
        cmap_colors = [(0, '#7094EC'), (1, 'red')]
        cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_colors)
        # print("3")
        g = nhm(data=df, cmapCenter=cmap, showyticks=False, linewidths=0)
        # print(4)
        g.hcluster(optimal_ordering=False)
        # print(5)
        df_ordered = g.data.iloc[g.rorder, g.corder]
        df_ordered.to_csv(f'{self.args.outdir}/ordered_heatmap.xls')
        # print(6)
        # print(g.dendrogram)
        fig, plots = g.run()
        # print(7)
        plt.savefig(f"{self.countsDir}/heatmap_{self.expCutoffsSuffix}.png")
        plt.savefig(f"{self.countsDir}/heatmap_{self.expCutoffsSuffix}.pdf")

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
            # print(splat)
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
        self.classificationRawCounts = {}
        self.mappingGroupsGTFDir = f'{self.countsDir}/mapping_groups_gtfs'

    def flatten_smorfs(self, genes):
        from Bio import SeqIO
        locus_cheker = []
        flattened = []
        for prot in genes:
            if '_F:' in prot:
                splat = prot.split(":")
                first = splat[0]
                if '+' in first:
                    chrom = f'+{first.split("+")[1]}'
                else:
                    chrom = f'-{first.split("-")[1]}'
                coords = splat[1].split("_")[0]
                locus = f'{chrom}:{coords}'
                if locus not in locus_cheker:
                    locus_cheker.append(locus)
                    flattened.append(prot)
        entries = []
        seqs = []
        if self.args.rescored:
            records = SeqIO.parse(self.rescoredMicroproteinsFasta, 'fasta')
        else:
            records = SeqIO.parse(self.uniqueMicroproteins, 'fasta')
        for record in records:
            if str(record.seq) not in seqs:
                seqs.append(str(record.seq))
                entries.append(str(record.description))
        filtered = []
        for i in flattened:
            if i in entries:
                filtered.append(i)
        return filtered

    def classify(self, method='union'):
        """
        method: union or separated
        """
        genes = self.df["Gene"].tolist()
        # print("total", len(genes))
        genes = self.flatten_smorfs(genes=genes)
        # print("flattened", len(genes))
        cols = self.df.columns[1:]
        for i, gene in enumerate(genes):
            if method == 'union':
                self.__union_grouping(i, gene, cols)
            if method == 'separated':
                higher = self.args.rpkm
                classification = "No coverage"
                for col in cols:
                    rpkms = self.df[col].tolist()
                    # if gene not in self.classification:
                    #     self.classification[gene] = {}
                    rpkm = float(rpkms[i])
                    if gene not in self.classification:
                        self.classification[gene] = {}
                    self.classification[gene][col] = rpkm
                    # if rpkm is None:
                    #     rpkm = 0
                    # print(rpkm)
                    if rpkm > 1 and col == "Default":
                        higher = rpkm
                        classification = col
                    if rpkm > higher and classification != 'Default':
                        higher = rpkm
                        classification = col

        if method == 'union':
            self.__save_union()
        if method == 'separated':
            self.__save_union()
        # print(self.classification)
        # gene_counts = collections.Counter(self.classification.values())
        # print(gene_counts)

    def __union_grouping(self, i, gene, cols):
        if gene not in self.classification:
            self.classification[gene] = {}
        for col in cols:
            rpkms = self.df[col].tolist()
            # if gene not in self.classification:
            #     self.classification[gene] = {}
            rpkm = float(rpkms[i])
            self.classification[gene][col] = rpkm

    def classify_raw_counts(self):
        files = os.listdir(self.rawCountsDir)
        for file in files:
            if file.endswith("counts.txt"):
                group = self.__check_counts_file_group(file)
                df = pd.read_csv(f'{self.rawCountsDir}/{file}', sep='\t', header=1)
                genes = df["Geneid"].tolist()
                genes = self.flatten_smorfs(genes)
                for col in df.columns[6:]:
                    counts = df[col].tolist()
                    for gene, count in zip(genes, counts):
                        if gene not in self.classificationRawCounts:
                            self.classificationRawCounts[gene] = {}
                        if group not in self.classificationRawCounts[gene]:
                            self.classificationRawCounts[gene][group] = []
                        self.classificationRawCounts[gene][group].append(count)

    @staticmethod
    def __check_counts_file_group(file):
        if 'ambiguous' in file:
            if 'multimappers' not in file:
                group = 'Amb'
            else:
                group = 'MM_Amb'
        else:
            if 'multimappers' in file:
                group = 'MM'
            else:
                group = 'Default'
        return group

    def __save_union(self):
        data = {'mapping_group': [], 'count': []}
        groups = {'MM': 0, 'MM_Amb': 0, 'Default': 0, 'Amb': 0, 'No coverage': 0}

        microprotein_groups = {'microprotein': []}

        for gene in self.classification:
            covered = False
            microprotein_groups['microprotein'].append(gene)
            for col in self.classification[gene]:

                if col not in microprotein_groups:
                    microprotein_groups[col] = []
                # print(self.classificationRawCounts[gene][col])
                arr = np.array(self.classificationRawCounts[gene][col])
                num = arr.max()
                microprotein_groups[col].append(num)
                if float(self.classification[gene][col]) >= self.args.rpkm:
                    if any(count >= self.args.minRawCounts for count in self.classificationRawCounts[gene][col]):
                        covered = True
                        groups[col] += 1

            if not covered:
                groups['No coverage'] += 1
        for group in groups:
            data['mapping_group'].append(group)
            data['count'].append(groups[group])
        df = pd.DataFrame(data=data)
        # if self.args.grouping_method == 'union':
        df.to_csv(self.mappingGroupsUnion, sep='\t', index=False)
        df.to_csv()
        # else:
        gdf = pd.DataFrame(data=microprotein_groups)
        gdf.to_csv(self.microproteinMappingGroups, sep='\t', index=False)
        #     df.to_csv(self.mappingGroups, sep='\t', index=False)

    def save_mapping_classification_union(self):
        df = pd.read_csv(self.microproteinMappingGroups, sep='\t')
        groups = []
        covered_smorfs = []
        mps = df["microprotein"].tolist()
        counts = {col: df[col].tolist() for col in df.columns[1:]}
        for i, mp in enumerate(mps):
            group = []

            for col in counts:
                if int(counts[col][i]) >= self.args.minRawCounts:
                    group.append(col)
            if group:
                groups.append(','.join(group))
            else:
                groups.append('No coverage')
            covered_smorfs.append(mp)
        data = {'smorf': covered_smorfs, 'group': groups}
        outdf = pd.DataFrame(data=data)
        outdf.to_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t', index=False)

    def save(self):
        data = {'smorf': [], 'group': []}
        for gene in self.classification:
            data['smorf'].append(gene)
            data['group'].append(self.classification[gene])
        df = pd.DataFrame(data=data)
        df.to_csv(self.mappingGroups, sep='\t', index=False)

    def save_multimappers(self):
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        mms = []
        amb = []
        for smorf, group in zip(smorfs, groups):
            group_list = group.split(",")
            if 'MM' in group:
                if 'Amb' not in group_list:
                    mms.append(smorf)

        if os.path.exists(self.rescoredMicroproteinsFasta):
            fasta = self.rescoredMicroproteinsFasta
        else:
            if os.path.exists(self.microproteinsBlast):
                fasta = self.microproteinsBlast
            else:
                fasta = self.uniqueMicroproteins
        records = SeqIO.parse(fasta, 'fasta')
        out_fasta = []
        for record in records:
            entry = str(record.description)
            if entry in mms:
                out_fasta.append(record)
        SeqIO.write(out_fasta, self.microproteinsMM, "fasta")

    def save_covered_gtf(self):
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs = df["smorf"].tolist()
        to_write = []
        with open(f'{self.countsDir}/rescored_smorfs_plus_reference_annotation.gtf', 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if any(smorf in line for smorf in smorfs):
                    to_write.append(line)
        with open(f'{self.countsDir}/smorfs_covered_by_riboseq.gtf', 'w') as outfile:
            outfile.writelines(to_write)

    def save_exclusive_classification(self):
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        classes = {'smorf': [], 'group': []}
        for smorf, group in zip(smorfs, groups):
            grouplist = group.split(',')
            if 'No coverage' in grouplist:
                orf_class = 'No coverage'
            elif 'Default' in grouplist:
                orf_class = 'Default'
            else:
                if 'Amb' in grouplist:
                    orf_class = 'Amb'
                else:
                    if 'MM' in grouplist:
                        orf_class = 'MM'
                    else:
                        if 'MM_Amb' in grouplist:
                            orf_class = 'MM_Amb'
            classes['smorf'].append(smorf)
            classes['group'].append(orf_class)
        edf = pd.DataFrame(data=classes)
        edf.to_csv(self.microproteinsMappingGroupsExclusive, sep='\t', index=False)

    def save_gtfs_for_each_mapping_group(self):
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        smorf_groups = {}
        for smorf, group in zip(smorfs, groups):
            smorf_groups[smorf] = group.replace(" ", "_").split(",")

        gtfs = {}

        with open(self.rescoredMicroproteinsGTF, 'r') as handler:
            lines = handler.readlines()

            for line in lines:
                cols = line.split("\t")
                attrs = cols[8].split(";")
                for a in attrs:
                    if 'transcript_id' in a:
                        mp = a.split(" ")[-1].replace("\"", "")
                        if mp in smorf_groups:
                            for group in smorf_groups[mp]:
                                if group not in gtfs:
                                    gtfs[group] = []
                                gtfs[group].append(line)
        self.check_dirs([self.mappingGroupsGTFDir])
        for gtf in gtfs:
            out = f'{self.mappingGroupsGTFDir}/{gtf}.gtf'
            with open(out, 'w') as outfile:
                outfile.writelines(gtfs[gtf])
