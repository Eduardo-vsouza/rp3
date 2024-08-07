import os
import sys

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure
from ..quantification import MOFF
from ..utils import ProtSplit


class FlashLFQ(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.print_row(word="quant mode")


        self.control = self.args.controlGroup
        self.metadata = {'FileName': [],
                         'Condition': [],
                         'Biorep': [],
                         'Fraction': [],
                         'Techrep': []}

        self.flashLFQDir = f'{self.outdir}/flashLFQ'
        self.mzMLCatDir = f'{self.flashLFQDir}/cat_mzml_files'
        self.comparisonsInputDir = f'{self.flashLFQDir}/comparisons_input'
        self.flashLFQComparisonsDir = f'{self.flashLFQDir}/comparisons'

        self.check_dirs([self.flashLFQDir, self.mzMLCatDir, self.comparisonsInputDir, self.flashLFQComparisonsDir])

        self.metadataFile = f'{self.mzMLCatDir}/ExperimentalDesign.tsv'

        self.flashLFQInputs = {}
        self.catFlashLFQInput = f'{self.flashLFQDir}/cat_flash_lfq_input.tsv'
        self.groupCompFolders = self.__define_outdirs()

        self.mzmlFileByGroup = {}
        self.__prepare_mzml()

    def iterate_groups(self):
        results, groups = self.args.results, self.args.groups
        for result, group in zip(results, groups):
            self.print_row(word=f"{group}", character="-")
            args = self.args
            args.outdir = result
            quant = MOFF(args=args, outdir=result)
            # quant.outdir = result
            quant.get_fdr_peptides()
            quant.generate_flash_lfq_input()
            self.flashLFQInputs[group] = quant.flashLFQInput
        self.print_row(word="comparisons", character="=")

    def __prepare_mzml(self):
        """
        FlashLFQ requires a metadata file to be present in the folder containing mzML files. As such, there can't be
        more than one folder containing mzML files. This function copies mzML files from their original folder and
        organizes the metadata in a way that's acceptable for FlashLFQ
        """
        j = 0
        print(f"--Preparing mzML files and metadata")
        for folder, group in zip(self.args.mzmlFolders, self.args.groups):

            j += 1
            files = os.listdir(folder)
            self.mzmlFileByGroup[group] = files
            self.__prepare_metadata(files=files, group=group, fraction=1)
            run = self.verify_checkpoint(outfile=f'{self.mzMLCatDir}/{files[0]}', step="")
            if run:
                cp = f'cp {folder}/* {self.mzMLCatDir}/.'
                os.system(cp)
        df = pd.DataFrame(self.metadata)
        df.to_csv(self.metadataFile, sep='\t', index=False)

    def __prepare_metadata(self, files, group, fraction):
        for i, file in enumerate(files):
            self.metadata['FileName'].append(file.split(".")[0])
            self.metadata['Condition'].append(group)
            self.metadata['Biorep'].append(i+1)
            self.metadata['Fraction'].append(fraction)
            self.metadata['Techrep'].append(1)

    def __define_outdirs(self):
        folders = {}
        for group in self.args.groups:
            if group != self.args.controlGroup:
                outdir = f'{self.flashLFQComparisonsDir}/{self.args.controlGroup}_x_{group}'
                self.check_dirs([outdir])
                folders[group] = outdir
        return folders

    def prepare_input(self):
        print(f"--Preparing input files")
        dfs = []
        for result, group in zip(self.args.results, self.args.groups):
            df = pd.read_csv(self.flashLFQInputs[group], sep='\t')
            dfs.append(df)
        cat_df = pd.concat(dfs)
        cat_df.to_csv(self.catFlashLFQInput, sep='\t', index=False)
        for result, group in zip(self.args.results, self.args.groups):
            if group != self.args.controlGroup:
                outfile = f'{self.comparisonsInputDir}/{self.args.controlGroup}_{group}_flashLFQ_input.tsv'
                run = self.verify_checkpoint(outfile=outfile, step=f"input preparing for {group}")
                if run:
                    files_to_include = self.mzmlFileByGroup[group] + self.mzmlFileByGroup[self.args.controlGroup]
                    df = cat_df[cat_df["File Name"].isin(files_to_include)]
                    df.to_csv(f'{self.comparisonsInputDir}/{self.args.controlGroup}_{group}_flashLFQ_input.tsv',
                              sep='\t', index=False)

    def run_flash_lfq(self):
        print(f"--Running FlashLFQ")
        for group in self.groupCompFolders:
            outfile = f'{self.groupCompFolders[group]}/BayesianFoldChangeAnalysis.tsv'
            run = self.verify_checkpoint(outfile=outfile, step=f"FlashLFQ quantification for {group}")
            self.print_row(word=f"{group} x {self.args.controlGroup}", character="-")
            if run:
                cmd = (f'{self.toolPaths["FlashLFQ"]} '
                       f'--idt {self.comparisonsInputDir}/{self.args.controlGroup}_{group}_flashLFQ_input.tsv'
                       f' --ppm 20 --rep {self.mzMLCatDir} --nor --ctr {group} '
                       f'--out {self.groupCompFolders[group]} --thr {self.args.threads} --bay')
                os.system(cmd)
        print(f"--Done. Results at {self.flashLFQDir}.")
        self.print_row(word="Finished")

    def split_microproteins(self):
        """
        Splits the results into three subsets: unannotated microproteins, annotated microproteins, and
        annotated standard-sized proteins.
        """
        print(f"--Splitting FlashLFQ results into protein groups")
        protsplit = ProtSplit(args=self.args)
        protsplit.split_protein_groups()
        prot_groups = protsplit.get_protein_groups(order='group_prot')  # dict: protein, group
        groups = os.listdir(self.flashLFQComparisonsDir)
        for group in groups:  # first get all proteins that passed the thresholds. We need their sequences
            outfile = f'{self.flashLFQComparisonsDir}/{group}/BayesianFoldChangeAnalysis.tsv'
            df = pd.read_csv(outfile, sep='\t')
            df = df[df["False Discovery Rate"] <= self.args.quantFDR]
            df = df[abs(df["Protein Log2 Fold-Change"]) > self.args.foldChangeCutoff]
            for prot_group in prot_groups:
                outfile = f'{self.flashLFQComparisonsDir}/{group}/{prot_group}_foldChangeAnalysis.csv'
                run = self.verify_checkpoint(outfile, step="foldChange analysis protein group splitting")
                if run:
                    gdf = df[df["Protein Group"].isin(prot_groups[prot_group])]
                    gdf.to_csv(outfile, sep='\t', index=False)
