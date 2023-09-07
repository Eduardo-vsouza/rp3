#!/usr/bin/python3

import os
import sys

import argparse

from src.pipeline import Pipeline
from src.metrics import DatabaseMetrics
# from src.pipeline_config import Configure


class RP3:
    def __init__(self):
        self.args = self.__get_args()

    def __get_args(self):
        self.main_parser = argparse.ArgumentParser(description="Ribosome Profiling and Proteogenomics Pipeline (RP3)",
                                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        self.mode_parser = self.main_parser.add_argument_group("Mode input options")

        self.mode_parser.add_argument("mode", metavar="Mode",
                                      help="Mode to run the pipeline for.\nList of Modes: "
                                           "database, "
                                           "ribocov, "
                                           "search, "
                                           "postms, "
                                           "quant, "
                                           "mods, "
                                           "id, "
                                           "id-quant, " 
                                           "extra, "
                                           "metrics, "
                                           "spectra,"
                                           "anno")

        self.args = self.main_parser.parse_args(sys.argv[1:2])
        self.mode = self.args.mode
        self.parser = argparse.ArgumentParser(description="Run pipeline_config in %s mode" % self.mode,
                                              prog=" ".join(sys.argv[0:2]),
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        self.general_args = self.parser.add_argument_group("General Parameters")
        self.general_args.add_argument("mode", metavar=self.mode)
        self.general_args.add_argument("--outdir", "-o", help="Inform the output directory", default="foopipe")
        self.general_args.add_argument("--threads", "-p", help="Number of threads to be used.", default=1)

        self.modeArguments = self.parser.add_argument_group(f"{self.mode} options")

        if self.mode == 'database':
            self.__set_database_mode()

        elif self.mode == 'translation':
            self.__set_translation_mode()

        elif self.mode == 'search' or self.mode == 'postms':
            self.__set_search_mode()

        elif self.mode == 'quant':
            self.__set_quant_mode()

        elif self.mode == 'search' or self.mode == 'mods':
            self.modeArguments.add_argument("--include_high_fdr", action='store_true')

        elif self.mode == 'id':
            self.__set_database_mode()
            self.__set_search_mode()

        elif self.mode == 'id-quant':
            self.__set_database_mode()
            self.__set_search_mode()
            self.__set_quant_mode()

        elif self.mode == 'extra':
            # extra filtering options
            self.modeArguments.add_argument("--refseq")  # blasts results against refseq or another database

        elif self.mode == 'metrics':
            self.__set_metrics_mode()

        elif self.mode == 'spectra':
            self.__set_spectra_mode()

        elif self.mode == 'anno':
            self.__set_anno_mode()

        elif self.mode == 'rescore':
            self.__set_rescore_mode()

        elif self.mode == 'ribocov':
            self.__set_ribocov_mode()

        # def check_mode(self):
        #     if self.mode == 'first_mode':
        #         self.__set_align_mode()
        #     elif self.mode == 'second_mode':
        #         self.__set_assembly_mode()

        args = self.parser.parse_args()
        return args

    def __set_translation_mode(self):
        self.modeArguments.add_argument("--genome")
        self.modeArguments.add_argument("--gtf_folder")
        self.modeArguments.add_argument("--ssh_folder", default='ribocov')

    def __set_database_mode(self):
        self.modeArguments.add_argument("--proteome", help="Reference proteome")
        self.modeArguments.add_argument("--genome")
        self.modeArguments.add_argument("--gtf_folder")
        self.modeArguments.add_argument("--ssh_folder")
        self.modeArguments.add_argument("--external_database")
        self.modeArguments.add_argument("--skip_translation", action='store_true')
        self.modeArguments.add_argument("--cat", action="store_true", help="Generate concatenated target and decoy "
                                                                           "databases.")

    def __set_search_mode(self):
        self.modeArguments.add_argument("--mzml")
        self.modeArguments.add_argument("--digest_max_length", default=7)
        self.modeArguments.add_argument("--digest_min_length", default=50)
        self.modeArguments.add_argument("--std_proteomics", action='store_true')
        self.modeArguments.add_argument("--quantify", action="store_true")
        self.modeArguments.add_argument("--mod")
        self.modeArguments.add_argument("--create_gtf", action="store_true")
        self.modeArguments.add_argument("--cat", action="store_true", help="Perform the search using a concatenated "
                                                                           "target and decoy database. Requires the "
                                                                           "databases to be generated using the 'cat' "
                                                                           "flag.")
        self.modeArguments.add_argument("--tmt_mod")
        self.modeArguments.add_argument("--fragment_mass_tolerance", default=20)
        self.modeArguments.add_argument("--refseq")  # blasts results against refseq or another database
        self.modeArguments.add_argument("--groups", help="Tab-delimited file associating a database with a raw file. "
                                                         "Should contain two columns: files, groups. Groups should have "
                                                         "the same name as the generated databases. If not specified,"
                                                         " the pipeline will search every raw file using every GTF "
                                                         "file provided.")

    def __set_quant_mode(self):
        self.modeArguments.add_argument("--moff_path", default="/home/eduardo/programs/moFF/moFF-2.0.3/moff_all.py")
        self.modeArguments.add_argument("--flash_lfq_path",
                                        default='dotnet /home/eduardo/programs/Flash-LFQ/CMD.dll')
        self.modeArguments.add_argument("--mzml")
        self.modeArguments.add_argument("--no_anno", action='store_true')
        self.modeArguments.add_argument("--not_blasted", action='store_true')
        self.modeArguments.add_argument("--spec_counts", action='store_true', help="Perform the search using spectral "
                                                                                   "counts only, and not MS1 peak "
                                                                                   "intensity. Does not require "
                                                                                   "--mzml.")

    def __set_metrics_mode(self):
        self.modeArguments.add_argument("--no_db", action="store_true")
        self.modeArguments.add_argument("--no_orfs", action='store_true')
        self.modeArguments.add_argument("--groups", help="If not specified, the group names will appear as the "
                                                         "database names. Provide a comma sep list.")

    def __set_assembly_mode(self):
        self.modeArguments.add_argument("--gtf", help="reference gtf file")
        self.modeArguments.add_argument("--strandness")
        self.modeArguments.add_argument("--step")

    def __set_spectra_mode(self):
        self.modeArguments.add_argument("--mzml")

    def __set_ribocov_mode(self):
        self.modeArguments.add_argument("--fastq", help="Provide the path to the folder containing "
                                                        "fastq files to be aligned to the genome. If the --aln "
                                                        "argument is provided, this is not necessary.")
        self.modeArguments.add_argument("--gtf", help="Reference gtf file containing coordinates for annotated genes. "
                                                      "The novel smORFs sequences from the proteogenomics analysis "
                                                      "will be appended to it.")
        self.modeArguments.add_argument("--genome_index", help="Path to the genome STAR index. If not provided, "
                                                         "it will use the human hg19 index available at /data/")
        self.modeArguments.add_argument("--cont_index", help="STAR index containing the contaminants "
                                                             "(tRNA/rRNA sequences). Reads mapped to these will be "
                                                             "excluded from the analysis.")
        self.modeArguments.add_argument("--aln", help="Folder containing bam or sam files with Ribo-Seq reads aligned to the"
                                                      " genome. In case this is provided, indexes are not required "
                                                      "and the alignment step will be skipped.")
        self.modeArguments.add_argument("--rpkm", help="RPKM cutoff to consider whether a smORF is sufficiently covered"
                                                       " by RPFs or not.", default=1)
        self.modeArguments.add_argument("--multimappings", help="max number of multimappings to be "
                                                                "allowed.", default=99)
        self.modeArguments.add_argument("--adapter", help="Provide the adapter sequence to be removed.",
                                        default="AGATCGGAAGAGCACACGTCT")
        self.modeArguments.add_argument("--plots", action="store_true")
        self.modeArguments.add_argument("--fastx_clipper_path",
                                        default=f"{sys.path[0]}/dependencies/fastx_clipper")
        self.modeArguments.add_argument("--fastx_trimmer_path",
                                        default=f'{sys.path[0]}/dependencies/fastx_trimmer')

    def __set_anno_mode(self):
        self.modeArguments.add_argument("--signalP", action="store_true")
        self.modeArguments.add_argument("--organism", default='eukarya')
        self.modeArguments.add_argument("--conservation", action='store_true')
        self.modeArguments.add_argument("--blast_db", default='/home/microway/blast_databases/homo_pan_mus_rat_equus_ovis_bos_sus_canis_database')
        self.modeArguments.add_argument("--rescored", action='store_true', help="Use this flag if the "
                                                                                "'rescore' mode was used to perform a "
                                                                                "second round of search using the "
                                                                                "results from the first search. Only "
                                                                                "the rescored microproteins will be "
                                                                                "analyzed for conservation in this "
                                                                                "case.")
    def __set_rescore_mode(self):
        self.modeArguments.add_argument("--mzml")
        self.modeArguments.add_argument("--proteome")
        self.modeArguments.add_argument("--msPattern")

    def __define_paths(self):
        pathconfig = Configure(args=self.args)
        self.toolPaths = pathconfig.read_config_file()

    def execute(self):
        pipe = Pipeline(args=self.args)
        # self.__define_paths()
        if self.mode == 'database':
            pipe.genome = self.args.genome
            pipe.proteome = self.args.proteome
            if not self.args.skip_translation:
                pipe.translate_in_silico()

            pipe.generate_databases()
            db = DatabaseMetrics(args=self.args)
            db.get_metrics()
            db.save()

        # elif self.mode == 'ribocov':
        #     pipe.genome = self.args.genome
        #     pipe.translate_in_silico()

        elif self.mode == 'search':
            pipe.search_mass_spec()
            pipe.post_process_searches()
            if self.args.mod is not None:
                pipe.check_mods()
            pipe.calculate_metrics()
            # if self.args.

        elif self.mode == 'postms':
            pipe.post_process_searches()
            # pipe.calculate_metrics()

        elif self.mode == 'quant':
            pipe.quantify()

        elif self.mode == 'id':
            pipe.genome = self.args.genome
            pipe.proteome = self.args.proteome
            pipe.translate_in_silico()
            pipe.generate_databases()
            pipe.search_mass_spec()
            pipe.post_process_searches()
            pipe.calculate_metrics()

        elif self.mode == 'metrics':
            pipe.calculate_metrics()

        elif self.mode == 'extra':
            pipe.extra_filters()

        elif self.mode == 'mods':
            pipe.check_mods()

        elif self.mode == 'spectra':
            pipe.annotate_spectra()

        elif self.mode == 'anno':
            pipe.annotate()

        elif self.mode == 'rescore':
            pipe.rescore()

        elif self.mode == 'ribocov':
            pipe.check_riboseq_coverage()


if __name__ == '__main__':
    data = RP3()
    data.execute()
