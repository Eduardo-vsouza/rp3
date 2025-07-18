import os

from src.search.peptide_search import PeptideSearch
from .post_process import PercolatorPostProcessing
from .database import Database
from .translation_ssh import GTFtoFasta
from .utils import Params, ProtSplit
from .extra_filter import ExtraFilter
from .mods import FormylSummarizer
from .search import PeptideReScoring
from .transcriptomics import StringTieAssembly
from .dualpg import DuoMetrics
from .paralogy import HomologyFinder
from .quant import FlashLFQ


class Pipeline:
    def __init__(self, args):
        self.args = args

        # files

        # directories
        self.outdir = args.outdir
        self.__check_dirs()


        # params
        self.threads = args.threads
        self.parameters = Params(args=self.args)


    def __check_dirs(self):
        if self.args.mode == 'rphub':
            import sys
            self.outdir = f'{sys.path[0]}/dump'

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.translationFolder = f'{self.outdir}/translation'

        if not os.path.exists(self.translationFolder):

            os.mkdir(self.translationFolder)

    def translate_in_silico(self):
        # assembly = f'{self.outdir}/transcriptomics/assembly/merged_assembly.gtf'
        # if os.path.exists(assembly):
            # self.args.gtf_folder = assembly
        translation = GTFtoFasta(folder=self.args.gtf_folder, genome=self.args.genome,
                                 local_outdir=self.translationFolder, args=self.args)
        translation.translate()
        self.parameters.add_mode_parameters(translation, args=self.args)
        self.parameters.update_params_file()

    def generate_databases(self):
        self.args.cat = True
        database = Database(outdir=self.outdir, reference_proteome=self.args.proteome,
                            translation_folder=self.translationFolder, external_database=self.args.external_database,
                            args=self.args)

        if self.args.external_database is not None:
            database.prepare_external_database()
        database.select_highly_homologous()
        database.append_reference()
        database.save_target_dbs()
        database.create_decoy_dbs()
        database.split_databases()
        self.parameters.add_mode_parameters(database, self.args)
        self.parameters.update_params_file()
        

    def search_mass_spec(self):
        self.postMSMode = self.args.postms_mode
        self.mzMLFolder = self.args.mzml

        quantify = False
        if self.args.quantify:
            quantify = True
        print("searching")

        # if self.args.engine == 'msfragger':
        search = PeptideSearch(quantify=quantify, args=self.args)

        if self.args.searchMode == 'cat':
            search.iterate_searches_cat()
        else:
            search.search_files()
        

        self.parameters.add_mode_parameters(search, self.args)
        self.parameters.update_params_file()

    def post_process_searches(self):
        if self.args.msBooster:
            from .spectra import Booster

            msb = Booster(args=self.args, rescore=False)
            msb.prepare_pin_files()
            msb.configure_parameters()
            msb.run()
            msb.merge_pin_files()

        self.postMSMode = self.args.postms_mode
        if not self.args.quantifyOnly:
            postms = PercolatorPostProcessing(args=self.args)

            if not self.args.recalculateFDR:  # this arg will re-do the post-processing after the percolator step
                if self.args.cascadeFDRmethod == 'sep':
                    postms.percolate_cascade()
                else:
                    if self.postMSMode == 'sep':
                        postms.percolate_single()
                    else:
                        postms.merge_all_pins()
                        
                        postms.percolate_all_pins()

            postms.fix_multiple_columns()
            postms.remove_annotated()

            """
            # postms.filter_microproteins()
            """
            postms.merge_fasta()
            if self.args.refseq is not None:
                self.extra_filters()

            if not self.args.std_proteomics:
                # postms.merge_gtf_results()
                # postms.create_merged_gtf()
                # postms.create_cds_gtf()
                # postms.get_utps()
                ...

            if self.args.rescore:
                self.args.msPattern = 'mzML'

                self.rescore()
            # counter = ORFCounter(args=self.args)
            # counter.count_smorfs()
            self.args.noRibocov = False

            from .metrics import MicroproteinCombiner

            summ = MicroproteinCombiner(args=self.args)
            summ.gather_microprotein_info()
            summ.save()
            protsplit = ProtSplit(outdir=self.args.outdir)
            protsplit.split_protein_groups()
            protsplit.split_fasta()

            from .metrics import MS2Summ

            ms2summ = MS2Summ(args=self.args)
            ms2summ.gather_peptides()
            ms2summ.gather_spectral_counts()
            ms2summ.gather_protein_sequences()
            ms2summ.save()


    def quantify(self):
        from .quantification import MOFF
        quant = MOFF(args=self.args)
        if self.args.spec_counts:
            quant.count_spectra()
            quant.count_spectra_annotated()
        else:
            quant.get_fdr_peptides()
            quant.generate_flash_lfq_input()
            quant.run_flashlfq()
    #     self.parameters.add_mode_parameters(quant, self.args)
    #     self.parameters.update_params_file()

    def duo_proteogenoimcs(self):
        duo = DuoMetrics(args=self.args)
        print("Calculating metrics for host and pathogen for Duo Proteogenomics")
        duo.separate_results()
        duo.count_orfs()

    def calculate_metrics(self):
        from .metrics import MicroproteinCombiner

        summ = MicroproteinCombiner(args=self.args)
        summ.gather_microprotein_info()
        summ.save()

    def annotate_spectra(self):

        if self.args.annotateSpectra:
            from .spectra import SpectrumAnnotator, Booster

            data = SpectrumAnnotator(args=self.args)
            data.prepare_input_files()
            data.annotate_spectra_parallel()

        if self.args.msBooster:
            from .spectra import SpectrumAnnotator, Booster

            msb = Booster(args=self.args)
            msb.prepare_pin_files()

        if self.args.butterfly:

            from .spectra import Butterfly
            butterfly = Butterfly(args=self.args)
            butterfly.fly()
            # plotter = ButterflyPlotter(['a.pin'], 'run.mzML', 'Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')\# 
            # df_pred = plotter.generate('AAAAAKAK', out_file='butterfly.png', annotate=True)

        # cov = ProteinCoverage(args=self.args)
        # cov.get_microprotein_sequences()
        # cov.get_mass_spec_peptides()
        # cov.plot()

    def extra_filters(self):
        extra = ExtraFilter(args=self.args)
        extra.blast_filter()
        extra.filter_results()

    def check_mods(self):
        mods = FormylSummarizer(args=self.args)
        # mods.get_modified_proteins()
        # mods.summarize_data()

    def annotate(self):
        if self.args.signalP:
            from .annotation import SignalP
            signal = SignalP(args=self.args)
            signal.run()
            signal.organize_files()
            signal.save_metrics()
        if self.args.conservation:
            from .annotation import Conservation
            conserv = Conservation(args=self.args)
            conserv.generate_non_redundant_fasta()
            conserv.blast_microproteins()
            conserv.parse_blast_results()
            conserv.create_evolview_input()
            conserv.generate_data_frame()
            # conserv.classify_conservation_by_mapping_groups()
        if self.args.uniprotTable:
            from .annotation import UniprotAnnotation
            anno = UniprotAnnotation(args=self.args)
            anno.filter_annotations()
            anno.intersect_mass_spec()
            anno.intersect_overall()
            anno.generate_plots()
        if self.args.orfClass:
            from .annotation import ORFClassification, ORFClassVis
            orfclass = ORFClassification(args=self.args)
            orfclass.intersect()
            orfclass.annotate()
            class_vis = ORFClassVis(args=self.args)
            file = class_vis.microproteinMappingGroupsForPlotsUnion
            if os.path.exists(file):
                class_vis.classify_by_groups()
                class_vis.get_annotations()
                class_vis.annotate_groups()
                # class_vis.plot_group_per_anno()
                class_vis.get_annotation_percentages()
                class_vis.plot_stacked_bar()
        if self.args.mhc:
            from .annotation import MHCDetector
            mhc = MHCDetector(args=self.args)
            mhc.run_mhc_flurry()
            mhc.filter_results()
        if self.args.paralogy:
            from .annotation import Homologs
            paralogy = Homologs(args=self.args)
            paralogy.blast()
            paralogy.parse(evalue=0.01)
            paralogy.save_clusters()
            paralogy.plot()
        if self.args.repeats:
            from .annotation import Repeater
            repeats = Repeater(args=self.args)
            repeats.fix_repeats_file()
            repeats.get_repeats()
            repeats.check_smorfs_within_repeats()
            repeats.plot_smorfs_in_repeats_based_on_clusters()
        if self.args.isoforms:
            from .annotation import Isoforms
            isoforms = Isoforms(args=self.args)
            isoforms.filter_gtf_features()
            isoforms.intersect()
            isoforms.read_intersected()
            isoforms.plot_cluster_overlaps()
        if self.args.deepLoc:
            from .annotation import DeepLoc
            deep_loc = DeepLoc(args=self.args)
            deep_loc.run()
            deep_loc.plot_results()


        # comparison = ResultsIntersection(args=self.args)
        # comparison.compare_groups()
        # summary = ResultsSummary(args=self.args)
        # summary.gather_smorf_data()
        # summary.get_signalp()
        # summary.get_conservation()
        # summary.get_riboseq_cov()
        # summary.save()

    def rescore(self):
        rescore = PeptideReScoring(args=self.args)
        if not self.args.quantifyOnly:
            rescore.generate_databases()
        rescore.re_search_peptides()
        if self.args.msBooster and not self.args.quantifyOnly:
            from .spectra import Booster

            msb = Booster(args=self.args)
            msb.prepare_pin_files()
            msb.configure_parameters()
            msb.run()
            msb.merge_pin_files()
        if self.args.postms_mode == 'sep' and not self.args.quantifyOnly:
            rescore.percolate_single()
        else:
            if not self.args.quantifyOnly:
                rescore.re_percolate_all_pins()
        if self.args.groupedFDR:
            if not self.args.quantifyOnly:
                rescore.re_assess_fdr_grouped()
        else:
            if not self.args.quantifyOnly:
                rescore.re_assess_fdr()
        if not self.args.quantifyOnly:
            rescore.merge_results()
            # rescore.filter_gtf()

    def check_riboseq_coverage(self):
        from .ribocov import FeatureCounts, RiboSeqCoverage, CoverageClassification, ORFCounter, RiboSeqAlign, SeqDist

        aln = RiboSeqAlign(args=self.args)
        if self.args.aln is None:
            aln.trim_reads()
            aln.run_trimming_parallel()
            aln.remove_ribosome()
            aln.align_rpfs()
        resc = PeptideReScoring(args=self.args)
        if self.args.rescored:
            rescored = True
        else:
            rescored = False
        resc.filter_gtf(rescored=rescored)

        feature_counts = FeatureCounts(args=self.args)
        feature_counts.append_reference_annotation()
        feature_counts.run_feature_counts()

        cov = RiboSeqCoverage(args=self.args)
        cov.feature_counts_to_rpkm()
        cov.plot_heatmaps(cutoff=self.args.rpkm)
        # cov.plot_histograms(folder=cov.rawCountsDir, header=1, title='Raw counts', cutoff=self.args.minRawCounts, log=True)
        cov.plot_histograms(folder=cov.rpkmDir, header=0, title='RPKM', cutoff=self.args.rpkm, log=True)

        cov_class = CoverageClassification(args=self.args)
        cov_class.classify_raw_counts()
        cov_class.classify(method=self.args.grouping_method)
        cov_class.save()
        cov_class.save_mapping_classification_union()
        cov_class.save_multimappers()
        cov_class.save_covered_gtf()
        cov_class.save_exclusive_classification()
        cov_class.save_gtfs_for_each_mapping_group()

        # if self.args.simulateMM:
            # mmcutoff = MMCutoff(args=self.args)
            # mmcutoff.sort_bam_files()
            # mmcutoff.intersect_alignments()
            # mmcutoff.index_intersected_bam_files()
            # mmcutoff.get_alignments()
            # mmcutoff.order_mappings()
            # mmcutoff.get_coverages()
            # mmcutoff.calculate_min_mm_for_coverage()
        if self.args.plots:
            counter = ORFCounter(args=self.args)
            counter.count_smorfs_union()

            # seqdist = SeqDist(args=self.args)
            # seqdist.get_genome_coverage()
            # seqdist.get_orf_coverage()
            # seqdist.plot_coverage()
            # seqdist.plot_coverage_single_orfs()

    def assemble_transcriptomes(self):
        assembly = StringTieAssembly(args=self.args)
        assembly.clean_reads()
        assembly.remove_contaminants()
        assembly.align_to_genome()
        assembly.assemble_transcriptome()
        assembly.merge_transcriptomes()

    def find_homologs(self):
        homo = HomologyFinder(args=self.args)
        print(f"-Identifying homologs in {self.args.genome} for the identified microproteins.")
        if self.args.tBlastN:
            homo.tblastn_mm_microproteins()
            homo.extract_aligned_sequences(blast='tblastn')
            homo.prepare_alignment_files(blast='tblastn')
            homo.perform_msa(blast='tblastn')
        if self.args.blastN:
            homo.get_mm_nucleotide_sequences()
            homo.blastn_mm_microproteins()
            homo.extract_aligned_sequences(blast='blastn')
            homo.prepare_alignment_files(blast='blastn')
            homo.perform_msa(blast='blastn')

    def compare_results(self):
        if self.args.specCounts:
            from .quant import SpecComparison

            comp = SpecComparison(args=self.args)
            comp.get_spec_counts()
            comp.add_overlaps()
            comp.create_data_frame()
            # comp.separate_up_regulated()
        if self.args.flashLFQ:
            quant = FlashLFQ(args=self.args)
            quant.iterate_groups()
            quant.prepare_input()
            quant.run_flash_lfq()
            quant.split_microproteins()
            quant.format_intermediate_heatmap_input()
            # quant.merge_conditions_heatmaps()
            # quant.generate_heatmaps()
            quant.erupt_volcanoes()


    def visualize_context(self):
        from .pgcontext import PGContext
        if self.args.pgViz:
            from .pgcontext import PGViz
            import tkinter as tk
            root = tk.Tk()
            interface = PGViz(root=root, args=self.args)
            root.mainloop()
        else:
            pgc = PGContext(args=self.args)
            pgc.expand_genes()
            pgc.intersect()
            pgc.gather_microprotein_sequences()
            pgc.gather_overlaps()
            pgc.gather_microproteins_data()
            pgc.define_smorf_limits()
            pgc.gather_ms_peptides()
            pgc.analyze_context()

    def create_rps(self):
        from .shiny import RPS
        rps = RPS(args=self.args)
        rps.visualize()

    def analyze_wgs(self):
        from .wgs import Variant
        variant = Variant(args=self.args)
        variant.pre_process_reads()
        variant.prepare_annotation_files()
        variant.align_reads()
        variant.convert_to_sorted_bam()
        variant.mark_duplicates()
        variant.recalibrate_base_score()


    def enter_rphub(self):
        from .rphub import RpHub
        rphub = RpHub(args=self.args)
        if self.args.results is not None:
            rphub.integrate_results()
        if self.args.summary:
            rphub.generate_summary()
        rphub.fetch_protein_seq()
        rphub.save()
