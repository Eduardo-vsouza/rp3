import os
import sys

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure


class DeepLoc(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
        self.deepLocDir = f'{self.outdir}/deeploc'
        self.check_dirs([self.deepLocDir])

    def run(self):
        if self.args.externalFasta:
            fasta = self.args.externalFasta
        else:
            fasta = self.select_fasta()

        print(f"--- Running DeepLoc on {fasta}...")
        cmd = f'{self.toolPaths["deeploc"]} -f {fasta} -o {self.deepLocDir} -m {self.args.deepLocModel}'
        os.system(cmd)

    def __get_sequences(self):
        if self.args.externalFasta:
            fasta = self.args.externalFasta
        else:
            fasta = self.select_fasta()

        dictio = {}

        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            entry = str(record.description)
            seq = str(record.seq)
            dictio[entry] = seq
        return dictio
    
    def plot_results(self):
        files = os.listdir(self.deepLocDir)

        seqs = self.__get_sequences()

        for file in files:
            if 'results' in file and file.endswith(".csv") and 'sequences' not in file:
                df = pd.read_csv(f'{self.deepLocDir}/{file}', sep=',')
                df = df[~df['Protein_ID'].str.contains("_ANNO")]
                proteins = df["Protein_ID"].tolist()
                locs = df["Localizations"].tolist()
                signals = df["Signals"].tolist()
                membrane_types = df["Membrane types"].tolist()
                sequences = [seqs[prot] for prot in proteins]
                df.insert(1, "sequence", sequences)
                df.to_csv(f'{self.deepLocDir}/deeploc_results_with_sequences.csv', sep=',', index=False)
                # for i, protein in enumerate(proteins):
                #     loc = locs[i]
                #     signal = signals[i]
                #     membrane_type = membrane_types[i]

                import matplotlib.pyplot as plt

                # Compute frequency counts for each category (unsplit)
                loc_counts = pd.Series(locs).value_counts()
                signal_counts = pd.Series(signals).value_counts()
                membrane_counts = pd.Series(membrane_types).value_counts()

                # Create subplots for three horizontal bar plots (unsplit)
                plt.figure(figsize=(8, 8))
                bars = plt.barh(loc_counts.index, loc_counts.values, color='skyblue', edgecolor='black')
                plt.title('Localization Distribution', fontsize=8)
                plt.yticks(fontsize=6)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Localization', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_localization.pdf', dpi=600)
                plt.close()

                plt.figure(figsize=(8, 8))
                bars = plt.barh(signal_counts.index, signal_counts.values, color='salmon', edgecolor='black')
                plt.title('Signals Distribution', fontsize=8)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Signals', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_signals.pdf', dpi=600)
                plt.close()

                plt.figure(figsize=(8, 8))
                bars = plt.barh(membrane_counts.index, membrane_counts.values, color='lightgreen', edgecolor='black')
                plt.title('Membrane Types Distribution', fontsize=8)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Membrane Type', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_membrane.pdf', dpi=600)
                plt.close()

                # Additional plots: split multiple entries per cell by '|'
                # For Localization column
                loc_split = pd.Series(
                    [val for entry in df["Localizations"].dropna() for val in str(entry).split('|')]
                )
                loc_split_counts = loc_split.value_counts()
                plt.figure(figsize=(8, 8))
                bars = plt.barh(loc_split_counts.index, loc_split_counts.values, color='skyblue', edgecolor='black')
                plt.title('Localization Split Distribution', fontsize=8)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Localization', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_localization_split.pdf', dpi=600)
                plt.close()

                # For Signals column
                signal_split = pd.Series(
                    [val for entry in df["Signals"].dropna() for val in str(entry).split('|')]
                )
                signal_split_counts = signal_split.value_counts()
                plt.figure(figsize=(8, 8))
                bars = plt.barh(signal_split_counts.index, signal_split_counts.values, color='salmon', edgecolor='black')
                plt.title('Signals Split Distribution', fontsize=8)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Signals', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_signals_split.pdf', dpi=600)
                plt.close()

                # For Membrane Types column
                membrane_split = pd.Series(
                    [val for entry in df["Membrane types"].dropna() for val in str(entry).split('|')]
                )
                membrane_split_counts = membrane_split.value_counts()
                plt.figure(figsize=(8, 8))
                bars = plt.barh(membrane_split_counts.index, membrane_split_counts.values, color='lightgreen', edgecolor='black')
                plt.title('Membrane Types Split Distribution', fontsize=8)
                plt.xlabel('Count', fontsize=8)
                plt.ylabel('Membrane Type', fontsize=8)
                for bar in bars:
                    width = bar.get_width()
                    plt.text(width + 0.1, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=8)
                ax = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.tight_layout()
                plt.savefig(f'{self.deepLocDir}/deepLoc_results_membrane_split.pdf', dpi=600)
                plt.close()
