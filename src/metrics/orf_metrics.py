import os
import sys

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

from ..pipeline_config import PipelineStructure


class ORFMetrics(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        # self.groups = self.args.groups

        self.metrics = {'subset': [], 'identifications': [], 'database': []}

    def get_metrics(self):
        folders = os.listdir(self.summarizedResultsDir)
        groups = 0
        for folder in folders:
            if os.path.isdir(f'{self.summarizedResultsDir}/{folder}'):
                if folder != 'merged':
                    groups += 1
        for folder in folders:
            if os.path.isdir(f'{self.summarizedResultsDir}/{folder}'):
                if (folder == 'merged' and groups > 1) or folder != 'merged':
                    db = folder.replace("_target_decoy_database.fasta", "").replace("_transcript_assembly", "")
                    db_dir = f'{self.summarizedResultsDir}/{folder}/'
                    files = os.listdir(db_dir)
                    for file in files:
                        if file.endswith(".fasta"):
                            subset = file.replace(".fasta", "")
                            self.__add_id_number(file=f'{db_dir}/{file}', db=db, subset=subset)

    def __add_id_number(self, file, subset, db):
        checker = []

        if os.path.exists(file):
            records = SeqIO.parse(file, 'fasta')
            ids = 0
            for record in records:
                if str(record.seq) not in checker:
                    ids += 1
                    checker.append(str(record.seq))
            self.metrics['subset'].append(subset)
            self.metrics['identifications'].append(ids)
            self.metrics['database'].append(db)

    def save(self):
        df = pd.DataFrame(data=self.metrics)
        df.to_csv(f'{self.metricsDir}/orf_metrics.xls', sep='\t', index=False)

    def plot(self):

        plt.rcParams.update({'font.size': 8})

        df = pd.read_csv(f'{self.metricsDir}/orf_metrics.xls', sep='\t')
        ax = sns.barplot(data=df, hue='database', x='identifications', y='subset', edgecolor='black')
        for i in ax.containers:
            ax.bar_label(i, )
        # ax.tick_params(axis='x', labelrotation=45)
        plt.tight_layout()
        # plt.legend(fontsize='x-small')
        # lgd = ax.legend(loc='upper left', bbox_to_anchor=(0, 4), mode="expand", borderaxespad=0.3)
        ax.autoscale_view()

        # plt.show()
        plt.savefig(f'{self.metricsDir}/orf_metrics.png', dpi=300)