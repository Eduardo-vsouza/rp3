import os
import sys

from Bio import SeqIO
import numpy as np

from ..utils import group_folder_generator, check_multiple_dirs


class DatabaseMetrics:
    def __init__(self, args):
        self.args = args

        self.outdir = self.args.outdir
        self.dbDir = f'{self.outdir}/databases'
        self.summDir = f'{self.outdir}/summarized_results'
        self.metricsDir = f'{self.outdir}/metrics'
        check_multiple_dirs([self.metricsDir])
        self.databasePaths = self.__get_database_paths()
        self.databases = {}


    def __get_database_paths(self):
        db_paths = []
        dbs = os.listdir(self.dbDir)
        for db in dbs:
            if db.endswith(".fasta"):
                db_paths.append(f'{self.dbDir}/{db}')
        return db_paths

    def get_metrics(self):
        for db in self.databasePaths:
            lenghts_predicted = []
            lengths_anno = []
            if db not in self.databases:
                self.databases[db] = {'annotated_target': 0, 'predicted_target': 0, 'entries_target': 0, 'entries_decoy' : 0,
                                      'entries_total': 0, 'median_length_annotated': 0, 'median_length_predicted': 0,
                                      'mean_length_annotated': 0, 'mean_length_predicted': 0, 'annotated_smORFs': 0}
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                entry = str(record.description)
                seq = str(record.seq)
                self.databases[db]['entries_total'] += 1

                # decoy
                if entry.startswith("rev_") or 'contaminant' in entry:
                    self.databases[db]['entries_decoy'] += 1

                # predicted target
                if '_F:' in entry or 'smORF' in entry:
                    if not entry.startswith("rev_"):
                        self.databases[db]['predicted_target'] += 1
                        self.databases[db]['entries_target'] += 1
                        lenghts_predicted.append(len(seq))

                # annotated target
                if '_ANNO' in entry:
                    if not entry.startswith("rev_") and 'contaminant' not in entry:
                        self.databases[db]['annotated_target'] += 1
                        self.databases[db]['entries_target'] += 1
                        lengths_anno.append(len(seq))
                        if len(seq) <= 150:
                            self.databases[db]['annotated_smORFs'] += 1
            # print(lengths_anno)
            # print(lenghts_predicted)
            self.databases[db]['mean_length_annotated'] = np.mean(lengths_anno)
            self.databases[db]['median_length_annotated'] = np.median(lengths_anno)
            self.databases[db]['mean_length_predicted'] = np.mean(lenghts_predicted)
            self.databases[db]['median_length_predicted'] = np.median(lenghts_predicted)

    def save(self):
        lines = []
        for db in self.databases:
            lines.append(f'Database name: {db.split("/")[-1]}\n')
            lines.append(f'Path: {db}\n')
            lines.append(f'Total entries: {self.databases[db]["entries_total"]}\n')

            #lines.append('Target database\n')
            lines.append(f'Target entries: {self.databases[db]["entries_target"]}\n')
            lines.append(f'Decoy entries: {self.databases[db]["entries_decoy"]}\n')

            lines.append("Annotated\n")
            lines.append(f'\tAnnotated entries: {self.databases[db]["annotated_target"]}\n')
            lines.append(f'\tAnnotated mean sequence length: {self.databases[db]["mean_length_annotated"]}\n')
            lines.append(f'\tAnnotated median sequence length: {self.databases[db]["median_length_annotated"]}\n')
            lines.append(f'\tAnnotated smORFs (<= 150 aa): {self.databases[db]["annotated_smORFs"]}\n')

            lines.append("Predicted\n")
            lines.append(f'\tPredicted entries: {self.databases[db]["predicted_target"]}\n')
            lines.append(f'\tPredicted mean sequence length: {self.databases[db]["mean_length_predicted"]}\n')
            lines.append(f'\tPredicted median sequence length: {self.databases[db]["median_length_predicted"]}\n')
            lines.append("----------------------------------------\n")
        with open(f'{self.metricsDir}/database_metrics.txt', 'w') as outfile:
            outfile.writelines(lines)


