import os
import sys
import subprocess

from ..utils import check_multiple_dirs



class PipelineStructure:
    def __init__(self, args):
        self.args = args
        self.pipelineDir = f'{sys.path[0]}'
        # self.check_dirs(self.outdir)///
        if self.args.mode == 'demo':
            self.testDir = f'{sys.path[0]}/demo_data'
            self.__set_demo_args()
        self.outdir = self.args.outdir
        self.check_dirs(self.outdir)
        self.translationDir = f'{self.outdir}/translation'

        # db
        self.databaseDir = f'{self.outdir}/databases'
        self.homologyDBDir = f'{self.outdir}/homology_database'
        self.repeatsDBDir = f'{self.outdir}/repeats_database'
        self.repeatsProteogenomicsDBs = f'{self.outdir}/proteogenomics_databases'

        # folders
        self.resultsDir = f'{self.outdir}/results'
        self.searchDir = f'{self.outdir}/peptide_search'
        self.postProcessDir = f'{self.outdir}/post_processing'
        self.quantDir = f'{self.outdir}/quantification'
        self.summarizedResultsDir = f'{self.outdir}/summarized_results'
        self.mergedResults = f'{self.summarizedResultsDir}/merged'
        self.spectrumDir = f'{self.outdir}/spectra'
        self.signalPDir = f'{self.mergedResults}/signalP'

        self.indexesDir = f'{self.pipelineDir}/data/STAR_indexes'

        self.logsDir = f'{self.outdir}/logs'
        # self.check_dirs([self.logsDir])
        # rescore directories
        self.rescoreDir = f'{self.outdir}/rescore'
        self.rescoreDatabaseDir = f'{self.rescoreDir}/databases'
        self.rescoreSearchDir = f'{self.rescoreDir}/peptide_search'
        self.rescorePostProcessDir = f'{self.rescoreDir}/post_processing'
        self.rescoreSummarizedResultsDir = f'{self.rescoreDir}/summarized_results'

        # counts directories
        self.countsDir = f'{self.outdir}/counts'
        self.rawCountsDir = f'{self.countsDir}/raw'
        self.rpkmDir = f'{self.countsDir}/rpkm'
        self.plotsDir = f'{self.countsDir}/plots'
        if self.args.mode == 'ribocov':
            self.check_dirs([self.countsDir, self.rawCountsDir, self.rpkmDir, self.plotsDir])

        # counts files
        self.mappingGroupsRPKMs = f'{self.countsDir}/mapping_groups_rpkm.txt'
        self.mappingGroupsRPKMsFiltered = f'{self.countsDir}/mapping_groups_rpkm_filtered.txt'
        self.mappingGroups = f'{self.countsDir}/mapping_groups.txt'
        self.mappingGroupsUnion = f'{self.countsDir}/mapping_groups_union.txt'
        self.microproteinMappingGroups = f'{self.countsDir}/microprotein_mapping_groups.txt'
        self.microproteinMappingGroupsForPlotsUnion = f'{self.countsDir}/microprotein_mapping_groups_plots_union.txt'
        self.microproteinsMappingGroupsExclusive = f'{self.countsDir}/microprotein_mapping_groups_exclusive.txt'
        self.microproteinsMM = f'{self.countsDir}/multimapping_smorfs.fasta'


        # files
        self.uniqueMicroproteins = f'{self.mergedResults}/microproteins_150.fasta'
        self.uniqueMicroproteinsNRFasta = f'{self.mergedResults}/nr_microproteins_150.fasta'
        self.utpsMicroproteins = f'{self.mergedResults}/microproteins_utps_150.fasta'
        self.utpsMicroproteinsBlast = f'{self.mergedResults}/microproteins_utps_150.fasta_blast_filt.fasta'
        self.microproteinsBlast = f'{self.mergedResults}/microproteins_150.fasta_blast_filt.fasta'
        self.uniqueMicroproteinsGTF = f'{self.mergedResults}/microproteins_150.gtf'
        self.metricsDir = f'{self.outdir}/metrics'
        self.rescoredMicroproteinsFasta = f'{self.rescoreSummarizedResultsDir}/filtered_rescored_microproteins_150.fasta'
        self.rescoredMicroproteinsGTF = f'{self.rescoreSummarizedResultsDir}/filtered_rescored_microproteins_150.gtf'
        self.mergedFullGTF = f'{self.mergedResults}/merged_predicted_microproteins.gtf'


        # transcriptomics
        self.rnaDir = f'{self.outdir}/transcriptomics'
        self.assemblyDir = f'{self.rnaDir}/assembly'
        self.rnaAlnDir = f'{self.rnaDir}/alignments'
        # self.rnaSortedBamDir = f'{self.rnaAlnDir}/sorted_alignments'
        # if self.mode == 'rna':
        self.genomeIndexDir = f'{self.rnaDir}/index'
        self.genomeIndex = f'{self.genomeIndexDir}/hisat_genome_index'
        self.spliceSites = f'{self.genomeIndexDir}/genome.ss'
        self.exons = f'{self.genomeIndexDir}/genome.exon'
        self.transcriptomeFasta = f'{self.rnaDir}/transcriptome.fasta'
        self.check_dirs([self.rnaDir])

        # ribo-seq
        self.riboSeqDir = f'{self.outdir}/ribo-seq'
        self.riboSeqTrimmedDir = f'{self.riboSeqDir}/trimmed_reads'
        self.riboSeqContaminantAlnDir = f'{self.riboSeqDir}/contaminant_alns'
        self.riboSeqAlnDir = f'{self.riboSeqDir}/alignments'
        self.riboSeqCovDir = f'{self.riboSeqDir}/coverage'
        self.riboSeqSortedBamDir = f'{self.riboSeqDir}/sorted_bam_alignments'
        if self.args.mode == 'ribocov':
            self.check_dirs([self.riboSeqDir, self.riboSeqCovDir, self.riboSeqSortedBamDir])
        # conservation
        self.phyloDir = f'{self.outdir}/phylogenetics'
        self.blastDir = f'{self.phyloDir}/blast'

        self.configFile = f'{self.pipelineDir}/config.txt'

        self.toolPaths = {}
        self.read_config_file()
        # print(self.toolPaths)

        # duo proteogenomics
        self.__set_duo_mode()
        self.__set_orf_class_params()
        self.expCutoffsSuffix = self.__check_suffix()

        # MSBooster
        self.boosterDir = f'{self.outdir}/MSBooster'
        self.boosterPinDir = f'{self.boosterDir}/pin_files'
        self.mergedBoosterPin = f'{self.boosterDir}/merged_booster.pin'

        self.splitMicroproteins = f'{self.rescoreSummarizedResultsDir}/microproteins_splitNuc.fasta'

        # orf class
        self.orfClassDir = f'{self.outdir}/orf_class'

        # mapping class
        # if self.args.paralogy:
        self.mappingClassDir = f'{self.outdir}/mapping_classification'
        self.mappingHomologyDir = f'{self.mappingClassDir}/homology'
        self.check_dirs([self.mappingClassDir, self.mappingHomologyDir])


    # def __set__postms_mode(self):
    #     self.percInputSingle = f'{db_path}/percolator_input_single'

    def __check_suffix(self):
        suffix = ''
        if self.args.mode == 'ribocov':
            suffix += f'RPKM-{self.args.rpkm}_rawCounts-{self.args.minRawCounts}'
        modes = ['ribocov', 'search', 'rescore', 'postms']
        if self.args.mode in modes:
            if self.args.proteinFDR:
                fdr = 'protein'
            else:
                fdr = 'peptide'
            suffix += f'_{fdr}FDR-{self.args.qvalue}'
        return suffix

    def __set_duo_mode(self):
        if self.args.mode == 'duo':
            self.duoDir = f'{self.outdir}/duo_proteogenomics'
            self.check_dirs([self.duoDir])
            self.hostMicroproteins = f'{self.duoDir}/{self.args.hostPattern}_microproteins.fasta'
            self.rescoredHostMicroproteins = f'{self.duoDir}/{self.args.hostPattern}_rescored_microproteins.fasta'
            self.pathogenMicroproteins = f'{self.duoDir}/{self.args.pathoPattern}_microproteins.fasta'
            self.rescoredPathogenMicroproteins = f'{self.duoDir}/{self.args.pathoPattern}_rescored_microproteins.fasta'

    def __set_demo_args(self):
        self.args.genome = f'{self.testDir}/genome.fasta'
        self.args.proteome = f'{self.testDir}/proteome.fasta'
        self.args.refseq = f'{self.testDir}/refseq.fasta'
        self.args.gtf = f'{self.testDir}/demo.gtf'
        self.args.mzml = f'{self.testDir}/search'
        # self.args.outdir = f'{sys.path[0]}/demo_outdir/'
        self.outdir = self.args.outdir
        self.args.gtf_folder = f'{self.testDir}/database'
        self.args.external_database = None
        self.args.external_gtf = None
        self.args.uniprotAnnotation = None
        self.args.highHomologyDB = False


        # search
        self.args.cat = True
        self.args.quantify = False
        self.args.mod = None
        self.args.groups = None
        self.args.tmt_mod = None
        self.args.fragment_mass_tolerance = 20
        self.args.std_proteomics = False
        self.args.msPattern = 'mzML'
        self.args.postms_mode = 'cat'
        self.args.hlaPeptidomics = False
        self.args.memory = 16
        self.args.recalculateFDR = False
        self.args.qvalue = 0.01
        self.args.msBooster = False
        self.args.smorfUTPs = False
        self.args.includeLowAnnotation = False
        self.args.proteinFDR = False
        self.args.minReplicates = 1
        self.args.rescore = True
        self.args.rescored = True

        # ribocov
        self.args.adapter = "AGATCGGAAGAGCACACGTCT"
        self.args.genome_index = f'{sys.path[0]}/STAR_indexes/hg19.star'
        self.args.cont_index = f'{sys.path[0]}/STAR_indexes/hg19cont.star'
        self.args.aln = None
        self.args.fastq = f'{sys.path[0]}/demo_data/ribocov'
        self.args.fastx_clipper_path = f"{sys.path[0]}/dependencies/fastx_clipper"
        self.args.fastx_trimmer_path = f"{sys.path[0]}/dependencies/fastx_trimmer"
        self.args.multimappings = 99
        self.args.rpkm = 1
        self.args.plots = True

        # rna
        self.args.strandness = 'rf'
        self.args.index = f'{self.pipelineDir}/data/hisat2_indexes/hg19/hisat_genome_index'
        self.args.reads_folder = f'{self.testDir}/rna'
        self.args.lib_type = 'paired'
        self.args.skip_trimming = False

    def exec(self, cmd):
        self.check_dirs([self.logsDir])
        log = f'{self.logsDir}/{self.args.mode}_log.txt'

        with open(log, 'a') as handler:
            process = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE)
            for c in iter(lambda: process.stdout.read(1), b""):
                sys.stdout.buffer.write(c)
                handler.buffer.write(c)

    def read_config_file(self):
        # print("yay")
        with open(self.configFile, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                # if 'MSBooster' in line:
                #     print(line)
                if not line.startswith("#"):
                    splat = line.split("\t")
                    if len(splat) > 1:
                        tool = splat[0]
                        path = splat[1].rstrip()
                        self.toolPaths[tool] = path
        return self.toolPaths

    def set_translation_attrs(self):
        self.multiMappingDir = f'{self.outdir}/multimappers'
        self.countsDir = f'{self.outdir}/counts'
        self.rawCountsDir = f'{self.countsDir}/raw'
        self.rpkmDir = f'{self.countsDir}/rpkm'
        self.plotsDir = f'{self.countsDir}/plots'
        # check_multiple_dirs([self.multiMappingDir])
        self.check_dirs([self.multiMappingDir, self.countsDir, self.rawCountsDir, self.rpkmDir, self.plotsDir])

    def set_ribocov_params(self):
        self.check_dirs([self.riboSeqDir, self.riboSeqTrimmedDir, self.riboSeqContaminantAlnDir, self.riboSeqAlnDir])


    @staticmethod
    def check_dirs(folders):
        if type(folders) == list:
            for folder in folders:
                if not os.path.exists(folder):
                    os.mkdir(folder)
        else:
            if not os.path.exists(folders):
                os.mkdir(folders)

    @staticmethod
    def get_content(main_dir):
        groups = os.listdir(main_dir)
        for group in groups:
            group_dir = f'{main_dir}/{group}'
            if os.path.isdir(group_dir):
                databases = os.listdir(group_dir)
                for db in databases:
                    if db.endswith("target_database.fasta") or db.endswith("_target_decoy_database.fasta"):
                        if os.path.isdir(f'{group_dir}/{db}'):
                            db_dir = f'{group_dir}/{db}'
                            files = os.listdir(db_dir)
                            for file in files:
                                fullfile = f'{db_dir}/{file}'
                                content = Content(file=file, fullfile=fullfile, db=db, group=group, main_dir=main_dir)
                                yield content

    def __set_orf_class_params(self):
        if self.args.mode == 'anno':
            if self.args.orfClass:
                self.orfClassDir = f'{self.outdir}/orf_class'

                self.check_dirs([self.orfClassDir])

    def select_fasta(self):
        if os.path.exists(self.rescoredMicroproteinsFasta):
            fasta = self.rescoredMicroproteinsFasta
        else:
            fasta = self.uniqueMicroproteins
        return fasta

    def verify_checkpoint(self, outfile, step):
        if os.path.exists(outfile):
            print(f"(!) Found output file {outfile}. Skipping {step}...")
            run = False
        else:
            run = True
        return run


class Content:
    def __init__(self, file, fullfile, group, db, main_dir):
        self.mainDir = main_dir

        self.group = group
        self.groupDir = f'{main_dir}/{group}'

        self.db = db
        self.dbDir = f'{self.groupDir}/{db}'

        self.file = file
        self.fullFile = fullfile


