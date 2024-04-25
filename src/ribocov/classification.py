import os
import sys
import collections

import pandas as pd


class CoverageClassification:
    def __init__(self, df):
        self.df = pd.read_csv(df, sep='\t')
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
        # print(self.classification)
        gene_counts = collections.Counter(self.classification.values())
        # print(gene_counts)

    def save(self, output):
        data = {'smorf': [], 'group': []}
        for gene in self.classification:
            data['smorf'].append(gene)
            data['group'].append(self.classification[gene])
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)