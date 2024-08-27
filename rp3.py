#!/usr/bin/python3

import os
import sys

import argparse

from src.pipeline import Pipeline
from src.metrics import DatabaseMetrics
from src.demo import Demo
# from src.pipeline_config import Configure

class StoreMultipleFiles(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)


class StoreMultipleGroups(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)


class RP3:
    def __init__(self):
        self.args = self.__get_args()

    def __get_args(self):
        self.main_parser = argparse.ArgumentParser(description="Ribosome Profiling and Proteogenomics Pipeline (RP3)\n",
                                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        self.mode_parser = self.main_parser.add_argument_group("Mode input options")

        self.mode_parser.add_argument("mode", metavar="Mode",
                                      help="Mode to run the pipeline for.\nList of Modes: "
                                           "database, "
                                           "search, "
                                           "postms, "
                                           "ribocov, "
                                           "id, "
                                           "id-quant, " 
                                           "extra, "
                                           "metrics, "
                                           "spectrum,"
                                           "anno")

        self.args = self.main_parser.parse_args(sys.argv[1:2])
        self.mode = self.args.mode
        self.parser = argparse.ArgumentParser(description="Run pipeline_config in %s mode" % self.mode,
                                              prog=" ".join(sys.argv[0:2]),
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        self.general_args = self.parser.add_argument_group("General Parameters")
        self.general_args.add_argument("mode", metavar=self.mode)
        self.general_args.add_argument("--outdir", "-o", help="Inform the output directory")
        self.general_args.add_argument("--threads", "-p", help="Number of threads to be used.", default=1,
                                       type=int)
        self.general_args.add_argument("--overwrite", action="store_true")
        self.general_args.add_argument("--genomeAssembly", help="available reference assemblies: hg38.")
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

        elif self.mode == 'spectrum':
            self.__set_spectra_mode()

        elif self.mode == 'anno':
            self.__set_anno_mode()

        elif self.mode == 'rescore':
            self.__set_rescore_mode()

        elif self.mode == 'ribocov':
            self.__set_ribocov_mode()

        elif self.mode == 'rna':
            self.__set_rna_mode()

        elif self.mode == 'demo':
            self.__set_demo_mode()

        elif self.mode == 'duo':
            self.__set_duo__mode()

        elif self.mode == 'homology':
            self.__set_homology_mode()

        elif self.mode == 'compare':
            self.__set_compare_mode()

        elif self.mode == 'pgc':
            self.__set_pgc_mode()

        elif self.mode == 'rps':
            self.__set_rps_mode()

        # def check_mode(self):
        #     if self.mode == 'first_mode':
        #         self.__set_align_mode()
        #     elif self.mode == 'second_mode':
        #         self.__set_assembly_mode()

        args = self.parser.parse_args()
        return args

    def __set_demo_mode(self):
        self.modeArguments.add_argument("--skip_rna", action='store_true')
        self.modeArguments.add_argument("--skip_database", action='store_true')
        self.modeArguments.add_argument("--skip_search", action='store_true')
        self.modeArguments.add_argument("--skip_postms", action='store_true')
        self.modeArguments.add_argument("--skip_ribocov", action='store_true')

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
        self.modeArguments.add_argument("--uniprotAnnotation", help="table containing uniprot entry, "
                                                                  "sequence and annotation levels. Only required "
                                                                    "if keeping uncharacterized microproteins "
                                                                    "with an annotation level <= --annotationLevel")
        self.modeArguments.add_argument("--annotationLevel", help="Annotation level cutoff for 'annotated'"
                                                                  " but uncharacterized proteins to be included in the "
                                                                  "analysis. Includes this and lower.", default=4)
        self.modeArguments.add_argument("--maxLength", help="Max ORF length.", default=150, type=int)
        self.modeArguments.add_argument("--highHomologyDB", action="store_true", help="Generate the "
                                                                                      "database including only "
                                                                                      "smORFs with high homology to"
                                                                                      " the rest of the genome.")
        self.modeArguments.add_argument("--minHomologs", default=10)
        self.altDBsArguments = self.parser.add_argument_group("Alternative database methods")
        self.altDBsArguments.add_argument("--repeats", action="store_true")
        self.altDBsArguments.add_argument("--repeatsFile")


    def __set_search_mode(self):
        self.modeArguments.add_argument("--proteinFDR", action='store_true')
        self.modeArguments.add_argument("--mzml")
        self.modeArguments.add_argument("--digest_max_length", default=7)
        self.modeArguments.add_argument("--digest_min_length", default=50)
        self.modeArguments.add_argument("--std_proteomics", action='store_true')
        self.modeArguments.add_argument("--quantify", action="store_true")
        self.modeArguments.add_argument("--quantifyOnly", action="store_true",
                                        help="Do not generate pin files for percolator. This is only relevant when "
                                             "running quantification later on after finishing all the analysis, and "
                                             "generating pin files for proteogenomics searches is not needed.")
        self.modeArguments.add_argument("--mod")
        self.modeArguments.add_argument("--create_gtf", action="store_true")
        self.modeArguments.add_argument("--cat", action="store_true", help="Perform the search using a concatenated "
                                                                           "target and decoy database. Requires the "
                                                                           "databases to be generated using the 'cat' "
                                                                           "flag.")
        self.modeArguments.add_argument("--tmt_mod")
        self.modeArguments.add_argument("--fragment_mass_tolerance", default=600)
        self.modeArguments.add_argument("--refseq")  # blasts results against refseq or another database
        self.modeArguments.add_argument("--groups", help="Tab-delimited file associating a database with a raw file. "
                                                         "Should contain two columns: files, groups. Groups should have "
                                                         "the same name as the generated databases. If not specified,"
                                                         " the pipeline will search every raw file using every GTF "
                                                         "file provided.")
        self.modeArguments.add_argument("--qvalue", help="Set the Percolator FDR cutoff to be applied "
                                                         "during the post-processing of the searches.", default=0.01,
                                        type=float)
        self.modeArguments.add_argument("--hlaPeptidomics", help="Sets parameters in the peptide search that are "
                                                                 "adequate for HLA peptidomics datasets.",
                                        action="store_true")
        self.modeArguments.add_argument("--rescore", action="store_true")
        self.modeArguments.add_argument("--msBooster", action="store_true", help="Run MSBooster on "
                                                                                 "MSFragger pin files to predict "
                                                                                 "retention times before running "
                                                                                 "Percolator.")
        self.modeArguments.add_argument("--recalculateFDR", action="store_true")
        self.modeArguments.add_argument("--groupedFDR", action="store_true", help="Assess the "
                                                                                  "FDR for predicted microproteins "
                                                                                  "and canonical proteins "
                                                                                  "separately.")
        self.modeArguments.add_argument("--postms_mode", help="cat or single.")
        self.modeArguments.add_argument("--smorfUTPs", action="store_true")
        self.modeArguments.add_argument("--proteome", help="required if rescoring")
        self.modeArguments.add_argument("--minReplicates", help="minimum number of replicates a peptide appear in "
                                                                "for it to be considered. It considers all "
                                                                "modifications for a peptide as the same peptide.",
                                        default=1, type=int)
        self.modeArguments.add_argument("--includeLowAnnotation", action="store_true", help="Whether to include "
                                                                                            "uniprot proteins "
                                                                                            "with low annotation "
                                                                                            "level. Requires database "
                                                                                            "to have been generated "
                                                                                            "with --uniprotAnnotation.")
        self.modeArguments.add_argument("--memory", help="RAM available for MSFragger.", default=64)
        self.modeArguments.add_argument("--keepAnnotated", action="store_true")
        self.modeArguments.add_argument("--phosphorylation", action="store_true")
        self.modeArguments.add_argument("--maxORFLength", type=int, default=150)

    def __set_quant_mode(self):
        self.modeArguments.add_argument("--moff_path", default="/home/eduardo/programs/moFF/moFF-2.0.3/moff_all.py")
        self.modeArguments.add_argument("--flash_lfq_path",
                                        default='FlashLFQ')
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
        self.modeArguments.add_argument("--noRibocov", action='store_true')

    def __set_assembly_mode(self):
        self.modeArguments.add_argument("--gtf", help="reference gtf file")
        # self.modeArguments.add_argument("--strandness")
        self.modeArguments.add_argument("--step")

    def __set_spectra_mode(self):

        self.modeArguments.add_argument("--annotateSpectra", help="Generate the annotated fragmentation "
                                                                  "spectra for each identified peptide.",
                                        action="store_true")
        self.modeArguments.add_argument("--predictRT", action="store_true", help="Predict retention times "
                                                                                 "(RT).")
        self.modeArguments.add_argument("--msBooster", action="store_true")
        self.annoArguments = self.parser.add_argument_group("Spectrum annotation parameters. Requires "
                                                            "--annotateSpectra")
        self.annoArguments.add_argument("--mzml")
        self.annoArguments.add_argument("--qvalue", default=0.01, type=float)

        self.rtArguments = self.parser.add_argument_group("RT prediction parameters. Requires --predictRT")
        self.rtArguments.add_argument("--protease", default="trypsin")
        self.rtArguments.add_argument("--minPepLen", default=7, type=int)
        self.rtArguments.add_argument("--maxPepLen", default=50, type=int)
        self.rtArguments.add_argument("--missCleavages", default=1)
        self.rtArguments.add_argument("--maxVarMods", default=1)
        self.rtArguments.add_argument("--precursorCharge", default=(2, 4), help="(min,max)")





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
                                                       " by RPFs or not.", default=1, type=float)
        self.modeArguments.add_argument("--minRawCounts", default=10, type=int)
        self.modeArguments.add_argument("--multimappings", help="max number of multimappings to be "
                                                                "allowed.", default=9999)
        self.modeArguments.add_argument("--adapter", help="Provide the adapter sequence to be removed.",
                                        default="AGATCGGAAGAGCACACGTCT")
        self.modeArguments.add_argument("--plots", action="store_true")
        self.modeArguments.add_argument("--fastx_clipper_path",
                                        default=f"{sys.path[0]}/dependencies/fastx_clipper")
        self.modeArguments.add_argument("--fastx_trimmer_path",
                                        default=f'{sys.path[0]}/dependencies/fastx_trimmer')
        self.modeArguments.add_argument("--grouping_method", help="Choose the way smORFS are "
                                                                   "assigned to their respective mapping groups. \n"
                                                                   "union,\n"
                                                                   "separated", default='union')
        self.modeArguments.add_argument("--skip_trimming", action="store_true")
        self.modeArguments.add_argument("--rescored", action="store_true")
        self.modeArguments.add_argument("--proteinFDR", action="store_true")
        self.modeArguments.add_argument("--qvalue", default=0.01, type=float)
        self.modeArguments.add_argument("--manualPlot", action="store_true")
        # self.modeArguments.add_argument("--simulateMM", action="store_true")
        # self.mmSimulationArguments = self.parser.add_argument_group("MM simulation parameters")
        # self.mmSimulationArguments.add_argument("--maxMismatches", default=2)
        # self.mmSimulationArguments.add_argument("--overwrite", action="store_true")



    def __set_anno_mode(self):
        self.modeArguments.add_argument("--signalP", action="store_true")
        self.modeArguments.add_argument("--organism", default='eukarya')
        self.modeArguments.add_argument("--conservation", action='store_true')
        self.modeArguments.add_argument("--blast_db", default='/home/microway/blast_databases/archita/homo_pan_mus_rat_equus_ovis_bos_sus_canis_macaca_loxodonta_balaenoptera_danio_drosophila_database')
        self.modeArguments.add_argument("--rescored", action='store_true', help="Use this flag if the "
                                                                                "'rescore' mode was used to perform a "
                                                                                "second round of search using the "
                                                                                "results from the first search. Only "
                                                                                "the rescored microproteins will be "
                                                                                "analyzed for conservation in this "
                                                                                "case.")
        self.modeArguments.add_argument("--uniprotTable")
        self.modeArguments.add_argument("--orfClass", action="store_true")
        self.modeArguments.add_argument("--paralogy", action="store_true")
        self.modeArguments.add_argument("--mhc", action="store_true")
        self.modeArguments.add_argument("--repeats", action="store_true")
        self.modeArguments.add_argument("--isoforms", action="store_true")
        self.modeArguments.add_argument("--exclusiveMappingGroups", action="store_true")
        self.__set_mhc_arguments()
        self.__set_paralogy_arguments()
        self.__set_orf_class_arguments()
        self.__set_repeats_arguments()
        self.__set_isoforms_arguments()

    def __set_isoforms_arguments(self):
        self.isoformsArguments = self.parser.add_argument_group("Isoforms parameters.")
        self.isoformsArguments.add_argument("--refGTF")

    def __set_repeats_arguments(self):
        self.repeatsArguments = self.parser.add_argument_group("Repeats parameters.")
        self.repeatsArguments.add_argument("--repeatsFile")


    def __set_paralogy_arguments(self):
        self.paralogyArguments = self.parser.add_argument_group("Paralogy parameters.")
        self.paralogyArguments.add_argument("--genome")
        self.paralogyArguments.add_argument("--alignToTranscriptome", action="store_true")
        self.paralogyArguments.add_argument("--maxMismatches", default=2, type=int)

    def __set_mhc_arguments(self):
        self.mhcArguments = self.parser.add_argument_group("MHC detection parameters.")
        self.mhcArguments.add_argument("--affinity", default=500)
        self.mhcArguments.add_argument("--affinityPercentile", default=0.02)
        self.mhcArguments.add_argument("--filterPipeResults", action="store_true")

    def __set_orf_class_arguments(self):
        self.orfClassArguments = self.parser.add_argument_group("ORF Classification parameters.")
        self.orfClassArguments.add_argument("--gtf", help="reference GTF file. For better accuracy in "
                                                          "annotation, this should be a GTF file from Ensembl. They "
                                                          "contain more terms that help better classifying the smORF.")

    def __set_rescore_mode(self):
        self.modeArguments.add_argument("--proteinFDR", action='store_true')
        self.modeArguments.add_argument("--groupedFDR", action="store_true", help="Assess the "
                                                                                  "FDR for predicted microproteins "
                                                                                  "and canonical proteins "
                                                                                  "separately.")
        self.modeArguments.add_argument("--mzml")
        self.modeArguments.add_argument("--proteome")
        self.modeArguments.add_argument("--msPattern", default='mzML')
        self.modeArguments.add_argument("--msBooster", action="store_true", help="Run MSBooster on "
                                                                                 "MSFragger pin files to predict "
                                                                                 "retention times before running "
                                                                                 "Percolator.")
        self.modeArguments.add_argument("--hlaPeptidomics", help="Sets parameters in the peptide search that are "
                                                                 "adequate for HLA peptidomics datasets.",
                                        action="store_true")
        # self.modeArguments.add_argument("--no_ambiguous", action="store_true", help="If turned on, microproteins "
        #                                                                             "with a peptide matching another "
        #                                                                             "predicted microprotein from "
        #                                                                             "the 3-ft will be discarded.")
        self.modeArguments.add_argument("--smorfUTPs", action="store_true")
        self.modeArguments.add_argument("--qvalue", default=0.01)
        self.modeArguments.add_argument("--postms_mode", help="cat, sep")
        self.modeArguments.add_argument("--keepAnnotated", action="store_true")
        self.modeArguments.add_argument("--phosphorylation", action="store_true")
        self.modeArguments.add_argument("--memory", default=32)
        self.modeArguments.add_argument("--maxORFLength", default=150)
        self.modeArguments.add_argument("--quantify", action="store_true")
        self.modeArguments.add_argument("--quantifyOnly", action="store_true")


    def __set_rna_mode(self):
        self.modeArguments.add_argument("--strandness", help="Inform the strandness of the experiment."
                                                             "--rf: fr-firststrand, \n"
                                                             "--fr: fr-secondstrand")
        self.modeArguments.add_argument("--gtf")
        self.modeArguments.add_argument("--genomeIndex")
        self.modeArguments.add_argument("--contIndex")
        self.modeArguments.add_argument("--pairedPattern", default='R')
        self.modeArguments.add_argument("--libType", help="Single or Paired.")
        self.modeArguments.add_argument("--fastq", help="Folder containing sequencing reads in fastq format.")
        # self.modeArguments.add_argument("genome_index", help="The genome index generated for hisat2 with the --ss "
        #                                                      "and --exon arguments. If not provided, a new index will "
        #                                                      "be generated. This might take a while.")

    def __set_duo__mode(self):
        self.modeArguments.add_argument("--hostPattern")
        self.modeArguments.add_argument("--pathoPattern")

    def __set_homology_mode(self):
        self.modeArguments.add_argument("--eValue", type=float, default=0.001)
        self.modeArguments.add_argument("--bitScore", type=int, default=50)
        self.modeArguments.add_argument("--genome")
        self.modeArguments.add_argument("--blastN", action="store_true")
        self.modeArguments.add_argument("--tBlastN", action="store_true")
        # self.modeArguments.add_argument("--minRawCounts", type=int, default=10)

    def __set_compare_mode(self):
        self.modeArguments.add_argument("--results", help="accepts multiple Rp3 output directories",
                                        nargs='+',
                                        action=StoreMultipleFiles)
        self.modeArguments.add_argument("--groups", help="accepts multiple Rp3 output group names",
                                        nargs='+',
                                        action=StoreMultipleFiles)

        self.modeArguments.add_argument("--rescored", action="store_true")
        self.modeArguments.add_argument("--predictedOnly", action="store_true")
        self.modeArguments.add_argument("--microproteins", action="store_true")
        self.modeArguments.add_argument("--flashLFQ", action="store_true", help="Run FlashLFQ "
                                                                                "to quantify proteins.")
        self.quantArguments = self.parser.add_argument_group("FlashLFQ arguments")
        # self.quantArguments.add_argument("--mzml")
        self.quantArguments.add_argument("--no_anno", action='store_true')
        self.quantArguments.add_argument("--specCounts", action='store_true', help="Perform the search using spectral "
                                                                                   "counts only, and not MS1 peak "
                                                                                   "intensity. Does not require "
                                                                                   "--mzml.")
        self.quantArguments.add_argument("--mzmlFolders", help="accepts multiple mzml folder paths. "
                                                              "Provide them in the same order as --results and "
                                                              "--groups.",
                                        nargs='+',
                                        action=StoreMultipleFiles)
        self.quantArguments.add_argument("--controlGroup")
        self.quantArguments.add_argument("--quantFDR", default=0.05, type=float)
        self.quantArguments.add_argument("--foldChangeCutoff", default=1, type=float)

    def __set_pgc_mode(self):
        self.modeArguments.add_argument("--gtf", help="gtf to intersect the Rp3 smORFs with")
        self.modeArguments.add_argument("--chromSizes", help="file with chromosome sizes")

        self.modeArguments.add_argument("--neighLength", default=500)
        self.modeArguments.add_argument("--pgViz", action="store_true")
        self.modeArguments.add_argument("--noPep", action="store_true")
        self.modeArguments.add_argument("--noPepSeq", action="store_true")

    def __set_rps_mode(self):
        self.modeArguments.add_argument("--results", help="accepts multiple Rp3 output directories",
                                        nargs='+',
                                        action=StoreMultipleFiles)
        self.modeArguments.add_argument("--groups", help="accepts multiple Rp3 output group names",
                                        nargs='+',
                                        action=StoreMultipleFiles)
        self.modeArguments.add_argument("--compare", action="store_true")
        self.modeArguments.add_argument("--deploy", action="store_true")

    def execute(self):
        if not self.mode == 'demo':
            pipe = Pipeline(args=self.args)
            os.system(f'date >> {self.args.outdir}/args.txt')
            for key, value in vars(self.args).items():
                os.system(f"echo '{key}: {value}' >> {self.args.outdir}/args.txt")

        if self.mode == 'database':
            pipe.genome = self.args.genome
            pipe.proteome = self.args.proteome
            if not self.args.skip_translation:
                pipe.translate_in_silico()

            pipe.generate_databases()
            db = DatabaseMetrics(args=self.args)
            db.get_metrics()
            db.save()

        elif self.mode == 'duo':
            pipe.duo_proteogenoimcs()

        elif self.mode == 'search':
            pipe.search_mass_spec()
            pipe.post_process_searches()
            if self.args.mod is not None:
                pipe.check_mods()
            # pipe.calculate_metrics()
            print("Finished peptide search and FDR assessment.")
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

        elif self.mode == 'spectrum':
            pipe.annotate_spectra()

        elif self.mode == 'anno':
            pipe.annotate()

        elif self.mode == 'rescore':
            pipe.rescore()

        elif self.mode == 'ribocov':
            pipe.check_riboseq_coverage()

        elif self.mode == 'rna':
            # self.__set_rna_mode()
            pipe.assemble_transcriptomes()

        elif self.mode == 'demo':
            demo = Demo(args=self.args)
            demo.test()
            demo.print_status()

        elif self.mode == 'homology':
            pipe.find_homologs()

        elif self.mode == 'compare':
            pipe.compare_results()

        elif self.mode == 'pgc':
            pipe.visualize_context()

        elif self.mode == 'rps':
            pipe.create_rps()

if __name__ == '__main__':
    print(" ____       _____\n"
          "|  _ \ _ __|___ /\n"
          "| |_) | '_ \ |_ \\\n"
          "|  _ <| |_) |__) |\n"
          "|_| \_\ .__/____/\n"
          "      |_|  ")
    print("RP3 v1.1.0")
    data = RP3()
    data.execute()

