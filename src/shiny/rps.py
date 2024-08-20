import os
import sys


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
        to be run.
        """
        fold_change_dfs = self.__get_fold_change_dfs()
        cmd = f'Rscript {self.shinyRScript} {fold_change_dfs} /home/microway/projects/microbiome_smORFs/reza_nanopore_metagenomes/muproteins/240729_T2D-S5_muProteIns/annotation/AIP_search/blasts/mps_families_MSA/msa'
        print(cmd)
        os.system(cmd)

    def __get_fold_change_dfs(self):
        compdir = f'{self.args.results[0]}/flashLFQ/comparisons'
        comps = os.listdir(compdir)
        df_paths = ''
        for comp in comps:
            annotated_mp = f'{compdir}/{comp}/annotated_microproteins_foldChangeAnalysis.csv'
            annotated_std = f'{compdir}/{comp}/standard_foldChangeAnalysis.csv'
            mp = f'{compdir}/{comp}/microproteins_foldChangeAnalysis.csv'
            dfs = [annotated_mp, annotated_std, mp]
            for df in dfs:
                df_paths += f'{df},'
        unique_df = f'{self.args.results[0]}/ip_unique.csv'
        enriched = f'{self.args.results[0]}/ip_upregulated.csv'
        full_comp = f'{self.args.results[0]}/group_comparison.csv'
        dfs = [unique_df, enriched, full_comp]
        for df in dfs:
            df_paths += f'{df},'
        return df_paths[:-1]


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
