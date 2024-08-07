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
        self.flashLFQHeatmapsDir = f'{self.flashLFQDir}/heatmaps'
        self.check_dirs([self.flashLFQDir, self.mzMLCatDir, self.comparisonsInputDir, self.flashLFQComparisonsDir,
                         self.flashLFQHeatmapsDir])

        self.metadataFile = f'{self.mzMLCatDir}/ExperimentalDesign.tsv'

        self.flashLFQInputs = {}
        self.catFlashLFQInput = f'{self.flashLFQDir}/cat_flash_lfq_input.tsv'
        self.groupCompFolders = self.__define_outdirs()
        self.proteinGroupsNames = ['standard', 'microproteins', 'annotated_microproteins']

        self.mzmlFileByGroup = {}
        self.__prepare_mzml()
        self.fileConditions = self.__associate_file_to_condition()

    def __associate_file_to_condition(self):
        """
        Reads metadata file (coldata) used by FlashLFQ and returns a dictionary {file, condition}
        """
        condict = {}
        df = pd.read_csv(self.metadataFile, sep='\t')
        files, conditions = df["FileName"].tolist(), df["Condition"].tolist()
        for file, condition in zip(files, conditions):
            condict[f'{file}.mzML'] = condition
        return condict

    def iterate_groups(self):
        results, groups = self.args.results, self.args.groups
        for result, group in zip(results, groups):
            self.print_row(word=f"{group}", character="-")
            args = self.args
            args.outdir = result
            quant = MOFF(args=args, outdir=result)
            # quant.get_fdr_peptides()
            # quant.generate_flash_lfq_input()
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
            outfile = self.__get_group_foldchange_output(group)
            df = pd.read_csv(outfile, sep='\t')
            # must include ones with 0 intensities in all replicates of a condition

            for prot_group in prot_groups:
                outfile = f'{self.flashLFQComparisonsDir}/{group}/{prot_group}_foldChangeAnalysis.csv'
                run = self.verify_checkpoint(outfile, step="foldChange analysis protein group splitting")
                if run:
                    intensity_df = self.__get_intensities_file(group)
                    missing_genes = self.__get_missing_genes(intensity_df)
                    df = df[(abs(df["Protein Log2 Fold-Change"]) > self.args.foldChangeCutoff) | (
                        df["Protein Group"].isin(missing_genes))]
                    df = df[
                        (df["False Discovery Rate"] <= self.args.quantFDR) | (df["Protein Group"].isin(missing_genes))]

                    gdf = df[df["Protein Group"].isin(prot_groups[prot_group])]
                    gdf.to_csv(outfile, sep='\t', index=False)

    def __get_missing_genes(self, df):
        """
        Gets the names of the genes with all intensities in a condition == 0. FlashLFQ does not include them in the
        foldChange analysis (apparently). We need to include them in the heatmap (sample-specific).
        """
        # df = pd.read_csv(dataframe, sep='\t')
        missing_genes = []
        genes = df["Protein Groups"].tolist()
        for i, gene in enumerate(genes):
            conditions = {}
            all_zero = False
            for col in df.columns:
                if 'Intensity_' in col:
                    intensities = df[col].tolist()
                    condition = self.fileConditions[f'{col.replace("Intensity_", "")}.mzML']
                    if condition not in conditions:
                        conditions[condition] = []
                    conditions[condition].append(intensities[i])
            for con in conditions:
                if all([intensity == 0 for intensity in conditions[con]]):
                    all_zero = True
            if all_zero:
                missing_genes.append(gene)
        return missing_genes




    def __get_group_foldchange_output(self, group):
        file = f'{self.flashLFQComparisonsDir}/{group}/BayesianFoldChangeAnalysis.tsv'
        return file

    def __get_protein_group_df(self, group, prot_group):
        file = f'{self.flashLFQComparisonsDir}/{group}/{prot_group}_foldChangeAnalysis.csv'
        return file

    def __get_intensities_file(self, group):
        file = f'{self.flashLFQComparisonsDir}/{group}/QuantifiedProteins.tsv'
        return file

    def __get_heatmap_df(self, group, prot_group):
        file = f'{self.flashLFQComparisonsDir}/{group}/{prot_group}_intensities_heatmap.csv'
        return file

    def format_intermediate_heatmap_input(self):
        print(f"--Formatting intensities tables for heatmaps")
        groups = os.listdir(self.flashLFQComparisonsDir)
        for group in groups:
            intensities_df_unformatted = pd.read_csv(self.__get_intensities_file(group), sep='\t')  # need to fix,
            # this depends on the foldchanage and the foldchange depends on this
            intensities_df = self.__format_intensities_file(group)
            # intensities_df = pd.read_csv(, sep='\t')  # get intensities in the counts table format
            missing_genes = self.__get_missing_genes(intensities_df_unformatted)
            for prot_group in self.proteinGroupsNames:

                fg_file = self.__get_protein_group_df(group=group, prot_group=prot_group)   # get foldChange for
                                                                                            # protein group
                fg_df = pd.read_csv(fg_file, sep='\t')  # get proteins passing the fold change and FDR cutoff
                prots = fg_df["Protein Group"].tolist()

                # filter by proteins passing the fold change and FDR cutoff
                df = intensities_df[(intensities_df["Protein Groups"].isin(prots)) | (intensities_df["Protein Groups"].isin(missing_genes))]
                heatmap_out = self.__get_heatmap_df(group, prot_group)
                df.to_csv(heatmap_out, sep='\t', index=False)

    def merge_conditions_heatmaps(self):
        print(f"--Merging heatmap input from different conditions")
        groups = os.listdir(self.flashLFQComparisonsDir)
        prot_group_dfs = {}
        for group in groups:
            for prot_group in self.proteinGroupsNames:
                if prot_group not in prot_group_dfs:
                    prot_group_dfs[prot_group] = []

                df = pd.read_csv(self.__get_heatmap_df(group, prot_group), sep='\t')
                prot_group_dfs[prot_group].append(df)
        for prot_group in prot_group_dfs:
            df = None
            for gdf in prot_group_dfs[prot_group]:
                if df is None:
                    df = gdf
                else:
                    df = pd.merge(df, gdf, on='Protein Groups')

            outfile = f'{self.flashLFQHeatmapsDir}/{prot_group}_heatmap.csv'
            df.to_csv(outfile, sep='\t', index=False)



    def __format_intensities_file(self, group):
        """

        """
        group_intensities_file = self.__get_intensities_file(group=group)  # get intensities and format
        group_intensities_df = pd.read_csv(group_intensities_file, sep='\t')  # for gene counts table format
        cols = [col for col in group_intensities_df.columns if 'Intensity' in col
                or "Protein Group" in col]  # get only prot names and intensities
        df = group_intensities_df[cols]
        df.rename(columns={col: col.replace("Intensity_", "") for col in df.columns})
        return df

