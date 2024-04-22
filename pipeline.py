import os

from src.search.peptide_search import MSFragger
from .post_process import PercolatorPostProcessing
from .database import Database
from .translation_ssh import GTFtoFasta
from .metrics import DatabaseMetrics, ORFMetrics, MicroproteinCombiner
from .utils import Params
from .extra_filter import ExtraFilter
from .mods import FormylSummarizer
from .spectra import SpectrumAnnotator, ProteinCoverage, RTPred, Booster
# from .annotation import SignalP, Conservation, ResultsSummary, UniprotAnnotation, ORFClassification, ORFClassVis, MHCDetector, Homologs, Repeater, Isoforms
from .annotation import *
from .search import PeptideReScoring
from .ribocov import FeatureCounts, RiboSeqCoverage, CoverageClassification, ORFCounter, RiboSeqAlign, SeqDist, MMCutoff
from .transcriptomics import StringTieAssembly
from .dualpg import DuoMetrics
from .paralogy import HomologyFinder


class Pipeline:
    def __init__(self, args):
        self.args = args

        # files

        # directories
        self.outdir = args.outdir
        self.translationFolder = f'{self.outdir}/translation'
        self.__check_dirs()


        # params
        self.threads = args.threads
        self.parameters = Params(args=self.args)


    def __check_dirs(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.exists(self.translationFolder):
            os.mkdir(self.translationFolder)

    def translate_in_silico(self):
        translation = GTFtoFasta(folder=self.args.gtf_folder, genome=self.args.genome,
                                 local_outdir=self.translationFolder)
        translation.translate()
        self.parameters.add_mode_parameters(translation, args=self.args)
        self.parameters.update_params_file()

    def generate_databases(self):
        database = Database(outdir=self.outdir, reference_proteome=self.args.proteome,
                            translation_folder=self.translationFolder, external_database=self.args.external_database,
                            args=self.args)
        if self.args.external_database is None:
            # database.unzip_assemblies()
            ...
        else:
            database.prepare_external_database()
        database.select_highly_homologous()
        database.append_reference()
        database.save_target_dbs()
        database.create_decoy_dbs()
        self.parameters.add_mode_parameters(database, self.args)
        self.parameters.update_params_file()

    def search_mass_spec(self):
        self.postMSMode = self.args.postms_mode
        self.mzMLFolder = self.args.mzml

        quantify = False
        if self.args.quantify:
            quantify = True
        print("searching")
        search = MSFragger(mzml_folder=self.mzMLFolder, outdir=self.outdir,
                           threads=self.threads, mod=self.args.mod, quantify=quantify, args=self.args)
        if self.postMSMode == 'cat':
            search.iterate_searches_cat()
        elif self.postMSMode == 'sep':
            search.iterate_searches_cat()
        else:
            search.iterate_searches_multi()
        self.parameters.add_mode_parameters(search, self.args)
        self.parameters.update_params_file()

    def post_process_searches(self):
        self.postMSMode = self.args.postms_mode

        postms = PercolatorPostProcessing(args=self.args)

        if not self.args.recalculateFDR:
            if self.postMSMode == 'cat':
                # postms.merge_pin_replicates()
                postms.merge_all_pins()
                # postms.percolate_multi()
                postms.percolate_all_pins()
            elif self.postMSMode == 'sep':
                postms.percolate_single()

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

    # def quantify(self):
        # quant = MOFF(args=self.args)
        # if self.args.spec_counts:
        #     quant.count_spectra()
            # quant.count_spectra_annotated()
    #     else:
    #         quant.get_fdr_peptides()
    #         quant.generate_flash_lfq_input()
    #         quant.run_flashlfq()
    #     self.parameters.add_mode_parameters(quant, self.args)
    #     self.parameters.update_params_file()

    def duo_proteogenoimcs(self):
        duo = DuoMetrics(args=self.args)
        print("Calculating metrics for host and pathogen for Duo Proteogenomics")
        duo.separate_results()
        duo.count_orfs()

    def calculate_metrics(self):
        # if self.args.mode == 'metrics':
        #     if not self.args.no_db:
        #         db = True
        #     else:
        #         db = False
        #     if not self.args.no_orfs:
        #         orfs = True
        #     else:
        #         orfs = False
        # else:
        #     db = True
        #     orfs = True
        # if db:
        #     db = DatabaseMetrics(args=self.args)
        #     db.get_metrics()
        #     db.save()
        # if orfs:
        #     orf = ORFMetrics(args=self.args)
        #     orf.get_metrics()
        #     orf.save()
        #     orf.plot()
        summ = MicroproteinCombiner(args=self.args)
        summ.gather_microprotein_info()
        summ.save()

    def annotate_spectra(self):
        if self.args.annotateSpectra:
            data = SpectrumAnnotator(args=self.args)
            data.prepare_input_files()
            data.annotate_spectra_parallel()
        if self.args.predictRT:
            rt = RTPred(args=self.args)
            rt.prepare_library()
            rt.predict_rts()
            rt.save_library()
            rt.compare_rts()
        if self.args.msBooster:
            msb = Booster(args=self.args)
            msb.prepare_pin_files()

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
        mods.get_modified_proteins()
        mods.summarize_data()

    def annotate(self):
        if self.args.signalP:
            signal = SignalP(args=self.args)
            signal.run()
            signal.organize_files()
            signal.save_metrics()
        if self.args.conservation:
            conserv = Conservation(args=self.args)
            conserv.generate_non_redundant_fasta()
            conserv.blast_microproteins()
            conserv.parse_blast_results()
            conserv.create_evolview_input()
            conserv.classify_conservation_by_mapping_groups()
        if self.args.uniprotTable:
            anno = UniprotAnnotation(args=self.args)
            anno.filter_annotations()
            anno.intersect_mass_spec()
            anno.intersect_overall()
            anno.generate_plots()
        if self.args.orfClass:
            orfclass = ORFClassification(args=self.args)
            orfclass.intersect()
            orfclass.annotate()
            class_vis = ORFClassVis(args=self.args)
            class_vis.classify_by_groups()
            class_vis.get_annotations()
            class_vis.annotate_groups()
            # class_vis.plot_group_per_anno()
            class_vis.get_annotation_percentages()
            class_vis.plot_stacked_bar()
        if self.args.mhc:
            mhc = MHCDetector(args=self.args)
            mhc.run_mhc_flurry()
            mhc.filter_results()
        if self.args.paralogy:
            paralogy = Homologs(args=self.args)
            paralogy.blast()
            paralogy.parse(evalue=0.01)
            paralogy.save_clusters()
            paralogy.plot()
        if self.args.repeats:
            repeats = Repeater(args=self.args)
            repeats.fix_repeats_file()
            repeats.get_repeats()
            repeats.check_smorfs_within_repeats()
            repeats.plot_smorfs_in_repeats_based_on_clusters()
        if self.args.isoforms:
            isoforms = Isoforms(args=self.args)
            isoforms.filter_gtf_features()
            isoforms.intersect()
            isoforms.read_intersected()
            isoforms.plot_cluster_overlaps()



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
        rescore.generate_databases()
        rescore.re_search_peptides()
        if self.args.msBooster:
            msb = Booster(args=self.args)
            msb.prepare_pin_files()
            msb.configure_parameters()
            msb.run()
            msb.merge_pin_files()
        if self.args.postms_mode == 'sep':
            rescore.re_percolate()
        else:
            rescore.re_percolate_all_pins()
        rescore.re_assess_fdr()
        rescore.merge_results()
        rescore.filter_gtf()

    def check_riboseq_coverage(self):
        aln = RiboSeqAlign(args=self.args)
        if self.args.aln is None:
            aln.trim_reads()
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
        cov.plot_histograms(folder=cov.rawCountsDir, header=1, title='Raw counts', cutoff=self.args.minRawCounts, log=True)
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
        #     mmcutoff = MMCutoff(args=self.args)
        #     mmcutoff.sort_bam_files()
        #     mmcutoff.intersect_alignments()
        #     mmcutoff.index_intersected_bam_files()
        #     mmcutoff.get_alignments()
        #     # mmcutoff.order_mappings()
        #     mmcutoff.get_coverages()
        #     mmcutoff.calculate_min_mm_for_coverage()
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
        assembly.check_index()
        assembly.align()
        assembly.assemble_transcriptomes()
        assembly.merge_assemblies()

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

