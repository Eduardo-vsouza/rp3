import os
import sys

import pandas as pd

from ..pipeline_config import PipelineStructure


class RPS(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.shinyRScript = f'{sys.path[0]}/src/shiny/shinyapp.R'

    def visualize(self):
        if self.args.compare:
            self.__visualize_multiple()
        else:
            self.__visualize_single()

    def __visualize_single(self):
        """
        Generates a shiny App to visualize results from a single Rp3 outdir
        """
        results_args, groups_args = self.__get_results_data_frames()
        print("res", results_args)
        cmd = f'Rscript {self.shinyRScript} {results_args} /home/microway/projects/microbiome_smORFs/reza_nanopore_metagenomes/muproteins/240729_T2D-S5_muProteIns/annotation/AIP_search/blasts/mps_families_MSA/msa'
        os.system(cmd)


    def __visualize_multiple(self):
        """
        Generates a shiny App to visualize and compare results from multiple Rp3 outdirs. Requires 'compare' mode
        to be run. In this case,
        """
        fold_change_dfs, alias = self.__get_fold_change_dfs()
        pgc_cols = ','.join([f'PGContext_{group}' for group in self.args.groups])
        print(pgc_cols)
        pgc_paths = ','.join([f'{result}/pg_context/context_figures' for result in self.args.results])
        cmd = f'Rscript {self.shinyRScript} {fold_change_dfs} {pgc_cols} {pgc_paths} {alias}'
        print(cmd)
        os.system(cmd)

    def __get_fold_change_dfs(self):
        # fold change comparison
        compdir = f'{self.args.outdir}/flashLFQ/comparisons'
        comps = os.listdir(compdir)
        df_paths = ''
        df_alias = ''
        for comp in comps:
            annotated_mp = f'{compdir}/{comp}/annotated_microproteins_foldChangeAnalysis.csv'
            annotated_std = f'{compdir}/{comp}/standard_foldChangeAnalysis.csv'
            mp = f'{compdir}/{comp}/microproteins_foldChangeAnalysis.csv'


            annotated_mp_pgc = self.__add_pg_context(file=annotated_mp, protein_col='Protein Group')
            annotated_std_pgc = self.__add_pg_context(file=annotated_std, protein_col='Protein Group')
            mp_pgc = self.__add_pg_context(file=mp, protein_col='Protein Group')

            annotated_mp_shiny = self.__filter_fold_change_df(annotated_mp_pgc)
            annotated_std_shiny = self.__filter_fold_change_df(annotated_std_pgc)
            mp_shiny = self.__filter_fold_change_df(mp_pgc)

            dfs = [annotated_mp_shiny, annotated_std_shiny, mp_shiny]
            for df in dfs:
                df_paths += f'{df},'
                splat = df.split("/")
                suffix = splat[-1].replace("foldChangeAnalysis_pgc_forShiny", "")
                alias = f'{splat[-2]}/{suffix}'
                df_alias += f'{alias},'

        # spec counts comparison
        files = os.listdir(self.args.outdir)
        unique_df = ''
        enriched = ''
        for file in files:
            if file.endswith("unique.csv"):
                unique_df_input = f'{self.args.outdir}/{file}'
                unique_df = self.__add_pg_context(unique_df_input)
            elif file.endswith("upregulated.csv"):
                enriched_input = f'{self.args.outdir}/{file}'
                enriched =  self.__add_pg_context(enriched_input)
        full_comp_input = f'{self.args.outdir}/group_comparison.csv'
        full_comp = self.__add_pg_context(full_comp_input)

        dfs = [unique_df, enriched, full_comp]
        for df in dfs:
            df_paths += f'{df},'
            alias = df.split("/")[-1].replace("_pgc", "")
            df_alias += f'{alias},'
        return df_paths[:-1], df_alias[:-1]

    def __filter_fold_change_df(self, file):
        """
        Removes unnecessary columns for shiny visualization
        """
        df = pd.read_csv(file, sep='\t')
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

    def __add_pg_context(self, file, protein_col='protein'):
        """
        Add paths to the appropriate PGC figures to display on the Shiny App.
        """
        df = pd.read_csv(file, sep='\t')
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
                if not added:
                    fig_paths[group].append('')
        for group in fig_paths:
            df.insert(4, f'PGContext_{group}', fig_paths[group])
        outfile = f'{file[:-4]}_pgc.csv'
        df.to_csv(outfile, sep='\t', index=False)
        return outfile



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
