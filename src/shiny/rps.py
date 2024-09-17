import os
import sys
import subprocess

import pandas as pd

from ..pipeline_config import PipelineStructure
from .deploy import RPSDeployer


class RPS(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)


    def visualize(self):
        self.print_row(word="RpS")
        if self.args.compare:
            self.__visualize_multiple()
        else:
            self.__visualize_single()

    def __visualize_single(self):
        """
        Generates a shiny App to visualize results from a single Rp3 outdir
        """
        results_args, groups_args = self.__get_results_data_frames()
        # print("res", results_args)
        cmd = f'Rscript {self.shinyRScript} {results_args} /home/microway/projects/microbiome_smORFs/reza_nanopore_metagenomes/muproteins/240729_T2D-S5_muProteIns/annotation/AIP_search/blasts/mps_families_MSA/msa'
        os.system(cmd)


    def __visualize_multiple(self):
        """
        Generates a shiny App to visualize and compare results from multiple Rp3 outdirs. Requires 'compare' mode
        to be run. In this case,
        """
        print(f"--Preparing data to deploy Shiny App")
        fold_change_dfs, alias = self.__get_fold_change_dfs()
        print(fold_change_dfs)
        pgc_cols = ','.join([f'PGContext_{group}' for group in self.args.groups])
        pgc_paths = ','.join([f'{result}/pg_context/context_figures' for result in self.args.results])
        msa_cols = ','.join([f'MSA_{group}' for group in self.args.groups])
        msa_paths = ','.join([f'{result}/homology/MSA_prot' for result in self.args.results])
        self.print_row(word="ShinyApp Deployment")
        if self.args.deploy:
            deploy = RPSDeployer(outdir=self.outdir,
                                 shiny_r=self.shinyRScript,
                                 paths_to_dfs=fold_change_dfs,
                                 image_cols=pgc_cols,
                                 image_dirs=pgc_paths,
                                 df_alias=alias,
                                 pdf_cols=msa_cols,
                                 pdf_dirs=msa_paths)
            deploy.copy_files()
            deploy.edit_rscript()
        else:
            print(f"--Deploying App")
            cmd = f'Rscript {self.shinyRScript} {fold_change_dfs} {pgc_cols} {pgc_paths} {alias} {msa_cols} {msa_paths}'
            print(cmd)
            os.system(cmd)

    def __get_fold_change_dfs(self):
        # fold change comparison
        print(f"--Gathering Rp3 data")
        compdir = f'{self.args.outdir}/flashLFQ/comparisons'
        comps = os.listdir(compdir)
        df_paths = ''
        df_alias = ''
        for comp in comps:
            self.print_row(word=comp, character="-")
            annotated_mp = pd.read_csv(f'{compdir}/{comp}/annotated_microproteins_foldChangeAnalysis.csv', sep='\t')
            annotated_std = pd.read_csv(f'{compdir}/{comp}/standard_foldChangeAnalysis.csv', sep='\t')
            mp = pd.read_csv(f'{compdir}/{comp}/microproteins_foldChangeAnalysis.csv', sep='\t')

            print(f"--Adding ribo-seq coverage")
            annotated_mp = self.__add_ribocov(annotated_mp, protein_col='Protein Group')
            annotated_std = self.__add_ribocov(annotated_std, protein_col='Protein Group')
            mp = self.__add_ribocov(mp, protein_col='Protein Group')

            print(f"--Adding PGContext")
            annotated_mp = self.__add_pg_context(df=annotated_mp, protein_col='Protein Group')
            annotated_std = self.__add_pg_context(df=annotated_std, protein_col='Protein Group')
            mp = self.__add_pg_context(df=mp, protein_col='Protein Group')

            print(f"--Adding MSA for paralogs")
            # if os.path.exists()
            annotated_mp = self.__add_paralogs_msa(df=annotated_mp, protein_col='Protein Group')
            annotated_std = self.__add_paralogs_msa(df=annotated_std, protein_col='Protein Group')
            mp = self.__add_paralogs_msa(df=mp, protein_col='Protein Group')

            print(f"--Adding protein names")
            annotated_mp = self.__add_protein_names(df=annotated_mp, protein_col='Protein Group')
            annotated_std = self.__add_protein_names(df=annotated_std, protein_col='Protein Group')
            mp = self.__add_protein_names(df=mp, protein_col='Protein Group')

            print(f"--Adding differential expression data")
            annotated_mp = self.__add_differential_expression(df=annotated_mp)
            annotated_std = self.__add_differential_expression(df=annotated_std)
            mp = self.__add_differential_expression(df=mp)

            print(f"--Adding MS1 quantification")
            annotated_mp_shiny = self.__filter_fold_change_df(df=annotated_mp,
                                                              file=f'{compdir}/{comp}/annotated_microproteins_foldChangeAnalysis.csv')
            annotated_std_shiny = self.__filter_fold_change_df(annotated_std,
                                                               file=f'{compdir}/{comp}/standard_foldChangeAnalysis.csv')
            mp_shiny = self.__filter_fold_change_df(mp,
                                                    file=f'{compdir}/{comp}/microproteins_foldChangeAnalysis.csv')

            dfs = [annotated_mp_shiny, annotated_std_shiny, mp_shiny]
            for df in dfs:
                df_paths += f'{df},'
                splat = df.split("/")
                suffix = splat[-1].replace("foldChangeAnalysis", "").replace("forShiny", "")
                suffix = suffix.replace("_pgc", "")[:-5]
                alias = f'{splat[-2]}/{suffix}'
                df_alias += f'{alias},'

        # spec counts comparison
        files = os.listdir(self.args.outdir)
        unique_df = ''
        enriched = ''
        for file in files:
            if file.endswith("unique.csv"):
                unique_df_input = pd.read_csv(f'{self.args.outdir}/{file}', sep='\t')
                unique_df_input = self.__add_ribocov(df=unique_df_input)
                unique_df_input = self.__add_paralogs_msa(df=unique_df_input)
                unique_df_input = self.__add_protein_names(df=unique_df_input)
                unique_df_input = self.__add_differential_expression(df=unique_df_input)
                unique_df = self.__add_pg_context(unique_df_input, save=f'{self.args.outdir}/{file}')
            elif file.endswith("upregulated.csv"):
                enriched_input = pd.read_csv(f'{self.args.outdir}/{file}', sep='\t')
                enriched_input = self.__add_ribocov(df=enriched_input)
                enriched_input = self.__add_paralogs_msa(df=enriched_input)
                enriched_input = self.__add_protein_names(df=enriched_input)
                enriched_input = self.__add_differential_expression(df=enriched_input)
                enriched =  self.__add_pg_context(enriched_input, save=f'{self.args.outdir}/{file}')

        full_comp_input = pd.read_csv(f'{self.args.outdir}/group_comparison.csv', sep='\t')
        full_comp_input = self.__add_ribocov(df=full_comp_input)
        full_comp_input = self.__add_paralogs_msa(df=full_comp_input)
        full_comp_input = self.__add_protein_names(df=full_comp_input)
        full_comp_input = self.__add_differential_expression(df=full_comp_input)
        print(f"--Adding Ribo-seq counts")
        # full_comp_input = self.__add_ribo_seq_counts(df=full_comp_input)
        full_comp = self.__add_pg_context(full_comp_input, save=f'{self.args.outdir}/group_comparison.csv')

        dfs = [unique_df, enriched, full_comp]
        for df in dfs:
            df_paths += f'{df},'
            alias = df.split("/")[-1].replace("_pgc", "").replace("_forShiny.csv", "")
            df_alias += f'{alias},'
        return df_paths[:-1], df_alias[:-1]

    def __filter_fold_change_df(self, df, file):
        """
        Removes unnecessary columns for shiny visualization
        """
        # df = pd.read_csv(file, sep='\t')
        to_drop = ['Standard Deviation of Peptide Log2 Fold-Changes',
                   'Bayes Factor',
                   'Gene',
                   'Organism',
                   'Uncertainty in Protein Log2 Fold-Change',
                   'Null Hypothesis Width',
                   'Control Condition',
                   'Treatment Condition']
        df = df.drop(to_drop, axis=1)
        columns_to_round = ['Protein Log2 Fold-Change',
                            'Protein Intensity in Control Condition',
                            'Protein Intensity in Treatment Condition',
                            ]
        df[columns_to_round] = df[columns_to_round].round(3)

        outfile = f'{file[:-4]}_forShiny.csv'
        df.to_csv(outfile, sep='\t', index=False)
        return outfile

    def __add_pg_context(self, df, protein_col='protein', save=None):
        """
        Add paths to the appropriate PGC figures to display on the Shiny App.
        """
        # df = pd.read_csv(file, sep='\t')
        prots = df[protein_col].tolist()
        fig_paths = {}
        for result, group in zip(self.args.results, self.args.groups):
            if group not in fig_paths:
                fig_paths[group] = []
            fig_dir = f'{result}/pg_context/context_figures'
            figures = os.listdir(fig_dir)
            for prot in prots:
                added = False
                for fig in figures:
                    if prot in fig:
                        fig_paths[group].append(os.path.abspath(f'{fig_dir}/{fig}'))
                        added = True
                        break
                if not added:
                    fig_paths[group].append('')
        for group in fig_paths:
            df.insert(4, f'PGContext_{group}', fig_paths[group])
        if save is not None:
            outfile = f'{save[:-4]}_pgc.csv'
            df.to_csv(outfile, sep='\t', index=False)
            return outfile
        else:
            return df

    def __add_ribocov(self, df, protein_col='protein'):
        # df = pd.read_csv(file, sep='\t')
        proteins = df[protein_col].tolist()
        mapping_groups_for_df = {}
        for results, group in zip(self.args.results, self.args.groups):
            mapping_groups_for_df[group] = []
            mappings_df = f'{results}/counts/microprotein_mapping_groups_plots_union.txt'
            if os.path.exists(mappings_df):
                mapping_groups = self.__get_microprotein_mapping_groups(mappings_df)
            else:
                mapping_groups = {}
            for protein in proteins:
                if protein in mapping_groups:
                    mapping_groups_for_df[group].append(mapping_groups[protein])
                else:
                    mapping_groups_for_df[group].append('No coverage')
        for group in mapping_groups_for_df:
            df.insert(4, f'mapping_groups_{group}', mapping_groups_for_df[group])
        return df

    def __get_microprotein_mapping_groups(self, df):
        mapping_groups = {}
        df = pd.read_csv(df, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            mapping_groups[smorf] = group
        return mapping_groups

    def __add_paralogs_msa(self, df, protein_col='protein'):
        """
        :param df: pandas data frame to be used in the Shiny app. This will add data for paralogs in the genome
        and paths to their MSA files
        """
        paths = {}
        proteins = df[protein_col].tolist()

        for results, group in zip(self.args.results, self.args.groups):

            if group not in paths:
                paths[group] = []
                paths[f'homologs_{group}'] = []
            for prot in proteins:
                added_fasta = False

                added_msa = False
                msa_dir = f'{results}/homology/MSA_prot'
                if os.path.exists(msa_dir):
                    files = os.listdir(msa_dir)
                    for file in files:
                        if prot in file:
                            if file.endswith(".pdf") and not added_msa:
                                paths[group].append(os.path.abspath(f'{msa_dir}/{file}'))
                                added_msa = True
                            if file.endswith(".fasta") and not added_fasta:
                                homologs = self.__count_homologs_in_fasta(msa_fasta=os.path.abspath(f'{msa_dir}/{file}'))
                                added_fasta = True
                                paths[f'homologs_{group}'].append(homologs)
                if not added_msa:
                    paths[group].append('')
                if not added_fasta:
                    paths[f'homologs_{group}'].append('')
        for group in paths:
            if 'homologs' not in group:
                df.insert(6, f"MSA_{group}", paths[group])
            else:
                df.insert(7, group, paths[group])
        return df

    def __count_homologs_in_fasta(self, msa_fasta):
        grep_process = subprocess.run(['grep', '>', msa_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      text=True)

        # Pass the grep output to wc -l
        wc_process = subprocess.run(['wc', '-l'], input=grep_process.stdout, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True)

        # Convert the output to an integer and subtract 1
        homologs = int(wc_process.stdout.strip()) - 1

        return homologs

    def __add_protein_names(self, df, protein_col='protein'):
        proteins = df[protein_col].tolist()
        names = []
        for prot in proteins:
            name = prot
            if '|' in prot:
                splat = prot.split("_")
                for a in splat:
                    if 'GN=' in a:
                        name = a.split("=")[-1]
            names.append(name)
        df.insert(1, "GeneName", names)
        return df

    def __add_differential_expression(self, df):
        de = f'{self.outdir}/de_analysis.csv'
        if os.path.exists(de):

            de_df = pd.read_csv(de, sep=',')
            col_names = ["Row.names", 'log2FoldChange', 'pvalue', 'padj']
            de_df = de_df.loc[:, col_names]

            h_dict = de_df.set_index('Row.names').to_dict(orient='index')
            cols = {col: [] for col in col_names if col != 'Row.names'}
            names = df["GeneName"].tolist()

            for mp in names:
                added = False
                if mp in h_dict:
                    added = True
                    for col in cols:
                        if col != 'Row.names':
                            # if type(h_dict[mp][col]) is int or type(h_dict[mp][col]) is float:
                            if h_dict[mp][col] > 0:
                                cols[col].append(h_dict[mp][col])
                            else:
                                if 'oldChange' in col:
                                    value = 0
                                else:
                                    value = 1
                                cols[col].append(value)
                if not added:
                    for col in cols:
                        if col != 'Row.names':
                            if 'oldChange' in col:
                                value = 0
                            else:
                                value = 1
                            cols[col].append(value)

            for col in cols:
                if col != 'Row.names':
                    df.insert(len(df.columns), col, cols[col])
        return df

    def __add_ribo_seq_counts(self, df, protein_col='protein'):
        proteins = df[protein_col].tolist()
        rpkms = {}
        rpkms_for_df = {}

        for results, group in zip(self.args.results, self.args.groups):
            rpkm_dir = f'{results}/counts/rpkm'
            rpkm_files = os.listdir(rpkm_dir)
            if group not in rpkms:
                rpkms[group] = {}
                rpkms_for_df[group] = {}
            for file in rpkm_files:
                mapping_group = self.__check_mapping_group(file)
                df = pd.read_csv(f'{rpkm_dir}/{file}', sep='\t')
                genes = df["Geneid"].tolist()
                cols = [df[col].tolist() for col in df.columns[1:]]
                for col in df.columns[1:]:
                    if col not in rpkms[group]:
                        rpkms[group][col] = {}
                    if mapping_group not in rpkms[group][col]:
                        rpkms[group][col][mapping_group] = {}
                    for i, gene in enumerate(genes):
                        for colist in cols:
                            rpkm = colist[i]
                            rpkms[group][col][mapping_group][gene] = rpkm
                # break
            for rep in rpkms[group]:
                added = False
                if rep not in rpkms_for_df[group]:
                    rpkms_for_df[group][rep] = {}
                for mg in rpkms[group][rep]:
                    if mg not in rpkms_for_df[group][rep]:
                        rpkms_for_df[group][rep][mg] = []
                    for protein in proteins:
                        added = False
                        if protein in rpkms[group][rep][mg]:
                            rpkm = rpkms[group][rep][mg][protein]
                            # added = True
                            rpkms_for_df[group][rep][mg].append(rpkm)
                            # break
                        # if not added:
                        else:
                            rpkms_for_df[group][rep][mg].append(0)
        for group in rpkms_for_df:
            for rep in rpkms_for_df[group]:
                for mg in rpkms_for_df[group][rep]:
                    print(group, rep, mg, len(rpkms_for_df[group][rep][mg]))
                    name = rep.split("/")[-1].replace(".fastq.gz_trimmedUnmapped.out.mate1Aligned.out.sam", "")
                    name = name.replace("aligned_to_genome_no_contaminant_", "")
                    df.insert(8, f'{name}|{group}|{mg}', rpkms_for_df[group][rep][mg])
        return df

    def __check_mapping_group(self, file):
        mg = ''
        if 'ambiguous' in file:
            if 'multimappers' in file:
                mg = 'MM_Amb'
            else:
                mg = 'Amb'
        if 'multimappers' in file:
            if 'ambiguous' not in file:
                mg = 'MM'
        if 'ambiguous' not in file and 'multimappers' not in file:
            mg = 'Default'
        return mg






    def __get_results_data_frames(self):
            paths = {}
            for outdir in self.args.results:
                paths[outdir] = []
                summarized_result = f'{outdir}/summarized_results/merged/microproteins_summary.csv'
                paths[outdir].append(summarized_result)
            results_args = ''
            groups_args = ''
            for outdir in paths:
                for folder in paths[outdir]:
                    results_args += f'{folder},'
                    groups_args += f'{outdir},'
            return results_args[:-1], groups_args[:-1]
