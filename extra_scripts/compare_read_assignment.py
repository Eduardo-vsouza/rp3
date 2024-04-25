import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


class ReadAssignment:
    def __init__(self, files, groups):
        self.files = files.split(",")
        self.groups = groups.split(",")


        self.features = ['Assigned', 'Unassigned_MultiMapping', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']
        self.data = {'group': []}
        self.totalReads = {}

        self.plotData = {'Reads': [], 'Read assignment': [], 'group': []}

    def read_files(self):
        for file, group in zip(self.files, self.groups):
            total_reads = {}

            with open(file, 'r') as handler:
                lines = handler.readlines()
                for line in lines:
                    line = line.rstrip()

                    cols = line.split('\t')

                    if line.startswith("Status"):
                        for i, col in enumerate(cols[1:]):
                            if i not in total_reads:
                                total_reads[i] = 0

                    if cols[0] in self.features:
                        if cols[0] not in self.data:
                            self.data[cols[0]] = []
                        for i, col in enumerate(cols[1:]):
                            total_reads[i] += int(col)
            with open(file, 'r') as handler:
                lines = handler.readlines()
                for line in lines:
                    line = line.rstrip()
                    cols = line.split('\t')
                    if cols[0] in self.features:
                        for i, col in enumerate(cols[1:]):
                            print(col, cols[0], group)
                            self.plotData['group'].append(group)
                            self.plotData['Reads'].append(int(col)/total_reads[i])
                            self.plotData['Read assignment'].append(cols[0])

    def organize_plot(self, output):
        print(self.data)
        for col in self.data:
            # if col == 'group':
            #     if 'group' not in self.plotData:
            #         self.plotData['group'] = []
            for row in self.data[col]:
                if col == 'group':
                    self.plotData['group'].append(row)
                else:
                    self.plotData['Read assignment'].append(col)
                    self.plotData['Reads'].append(row)
        # print(self.plotData)
        df = pd.DataFrame(self.plotData)
        df.to_csv(output, sep='\t', index=False)

        # print(self.data)
    def plot(self):
        ax = sns.boxplot(data=self.plotData, x='Read assignment', y='Reads', hue='group')
        plt.show()


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: .py <files> <groups> <outdir>")
    else:
        data = ReadAssignment(files=sys.argv[1], groups=sys.argv[2])
        data.read_files()
        # data.organize_plot(output=sys.argv[3])
        data.plot()