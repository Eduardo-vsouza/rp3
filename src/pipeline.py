import os
import sys

from .peptide_search import MSFragger
from .post_process import PercolatorPostProcessing
from .database import Database
from .translation_ssh import GTFtoFasta
from .quantification import MOFF
from .metrics import DatabaseMetrics, ORFMetrics
from .utils import Params
from .extra_filter import ExtraFilter
from .mods import FormylSummarizer
from .spectra import SpectrumAnnotator
from .annotation import SignalP, Conservation, ResultsIntersection, ResultsSummary
from .search import PeptideReScoring
from .ribocov import FeatureCounts, RiboSeqCoverage, CoverageClassification, ORFCounter, RiboSeqAlign


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
                                 local_outdir=self.translationFolder, ssh_folder=self.args.ssh_folder)
        translation.upload_files()
        self.parameters.add_mode_parameters(translation, args=self.args)
        self.parameters.update_params_file()

    def generate_databases(self):
        database = Database(outdir=self.outdir, reference_proteome=self.args.proteome,
                            translation_folder=self.translationFolder, external_database=self.args.external_database,
                            args=self.args)
        if self.args.external_database is None:
            database.unzip_assemblies()
        else:
            database.prepare_external_database()
        database.append_reference()
        database.save_target_dbs()
        database.create_decoy_dbs()
        self.parameters.add_mode_parameters(database, self.args)
        self.parameters.update_params_file()

    def search_mass_spec(self):
        self.mzMLFolder = self.args.mzml

        quantify = False
        if self.args.quantify:
            quantify = True
        print("searching")
        search = MSFragger(mzml_folder=self.mzMLFolder, outdir=self.outdir,
                           threads=self.threads, mod=self.args.mod, quantify=quantify, args=self.args)
        # search.iterate_searches()
        if self.args.cat:
            search.iterate_searches_cat()
        else:
            search.iterate_searches_multi()
        self.parameters.add_mode_parameters(search, self.args)
        self.parameters.update_params_file()

    def post_process_searches(self):
        postms = PercolatorPostProcessing(args=self.args)
        postms.merge_pin_replicates()
        if self.args.cat:
            postms.percolate_multi()
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
        # counter = ORFCounter(args=self.args)
        # counter.count_smorfs()

    # def quantify(self):
    #     quant = MOFF(args=self.args)
    #     if self.args.spec_counts:
    #         # quant.count_spectra()
    #         quant.count_spectra_annotated()
    #     else:
    #         quant.get_fdr_peptides()
    #         quant.generate_flash_lfq_input()
    #         quant.run_flashlfq()
    #     self.parameters.add_mode_parameters(quant, self.args)
    #     self.parameters.update_params_file()

    def calculate_metrics(self):
        if self.args.mode == 'metrics':
            if not self.args.no_db:
                db = True
            else:
                db = False
            if not self.args.no_orfs:
                orfs = True
            else:
                orfs = False
        else:
            db = True
            orfs = True
        if db:
            db = DatabaseMetrics(args=self.args)
            db.get_metrics()
            db.save()
        if orfs:
            orf = ORFMetrics(args=self.args)
            orf.get_metrics()
            orf.save()
            orf.plot()

    def annotate_spectra(self):
        data = SpectrumAnnotator(args=self.args)
        data.prepare_input_files()
        data.annotate_spectra()

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
        # comparison = ResultsIntersection(args=self.args)
        # comparison.compare_groups()
        summary = ResultsSummary(args=self.args)
        summary.gather_smorf_data()
        summary.get_signalp()
        summary.get_conservation()
        summary.get_riboseq_cov()
        summary.save()


    def rescore(self):
        rescore = PeptideReScoring(args=self.args)
        rescore.generate_databases()
        rescore.re_search_peptides()
        rescore.re_percolate()
        rescore.re_assess_fdr()
        rescore.merge_results()
        rescore.filter_gtf()

    def check_riboseq_coverage(self):
        aln = RiboSeqAlign(args=self.args)
        if self.args.aln is None:
            aln.trim_reads()
            aln.remove_ribosome()
            aln.align_rpfs()
        feature_counts = FeatureCounts(args=self.args)
        feature_counts.append_reference_annotation()
        feature_counts.run_feature_counts()
        cov = RiboSeqCoverage(args=self.args)
        cov.feature_counts_to_rpkm()
        cov.plot_heatmaps(cutoff=self.args.rpkm)
        cov_class = CoverageClassification(args=self.args)
        cov_class.classify()
        cov_class.save()
        if self.args.plots:
            counter = ORFCounter(args=self.args)
            counter.count_smorfs()

