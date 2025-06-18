import os
import sys
import subprocess

import pandas as pd

from ..utils import check_multiple_dirs



class PipelineStructure:
    def __init__(self, args):
        self.args = args
        self.pipelineDir = f'{sys.path[0]}'
        # self.check_dirs(self.outdir)///
        if self.args.mode == 'demo':
            self.testDir = f'{sys.path[0]}/demo_data'
            self.set_demo_args()
        self.dataFolder = f'{sys.path[0]}/data'
        self.refAnnoDir = f'{self.dataFolder}/reference_annotations'
        self.__define_genomes()

        self.outdir = self.args.outdir
        if self.args.mode == 'rphub':
            self.outdir = self.args.rpHubDir
        self.check_dirs([self.outdir])
        self.translationDir = f'{self.outdir}/translation'

        # db
        self.databaseDir = f'{self.outdir}/databases'
        self.homologyDBDir = f'{self.outdir}/homology_database'
        self.repeatsDBDir = f'{self.outdir}/repeats_database'
        self.repeatsProteogenomicsDBs = f'{self.outdir}/proteogenomics_databases'
        self.refProteome = f'{self.databaseDir}/proteome.fasta'
        self.refProteomeWithDecoy = f'{self.databaseDir}/decoy_proteome.fasta'

        # split db
        self.splitDbDir = f'{self.databaseDir}/split_databases'
        self.splitDbProteomeDir = f'{self.splitDbDir}/proteome'
        self.splitDbProteogenomicsDir = f'{self.splitDbDir}/proteogenomics'
        self.fullContaminantsDb = f'{self.splitDbDir}/contaminants.fasta'

        # folders
        self.resultsDir = f'{self.outdir}/results'
        self.searchDir = f'{self.outdir}/peptide_search'
        self.postProcessDir = f'{self.outdir}/post_processing'
        self.quantDir = f'{self.outdir}/quantification'
        self.summarizedResultsDir = f'{self.outdir}/summarized_results'
        self.mergedResults = f'{self.summarizedResultsDir}/merged'
        self.spectrumDir = f'{self.outdir}/spectra'
        self.signalPDir = f'{self.outdir}/signalP'
        self.signalPstandardDir = f'{self.signalPDir}/standardSizes_proteins'
        self.signalPMicroproteinDir = f'{self.signalPDir}/microproteins'
        self.signalPAnnoMPDir = f'{self.signalPDir}/annotated_microproteins'

        self.splitFastaDir = f'{self.outdir}/split_fasta'

        self.indexesDir = f'{self.pipelineDir}/data/STAR_indexes'

        self.logsDir = f'{self.outdir}/logs'
        # self.check_dirs([self.logsDir])

        # rescore directories
        self.rescoreDir = f'{self.outdir}/rescore'
        self.rescoreDatabaseDir = f'{self.rescoreDir}/databases'
        self.rescoreSearchDir = f'{self.rescoreDir}/peptide_search'
        self.rescorePostProcessDir = f'{self.rescoreDir}/post_processing'
        self.rescoreSummarizedResultsDir = f'{self.rescoreDir}/summarized_results'
        self.rescoreGroupFDRDir = f'{self.rescoreDir}/group_FDR'
        self.rescoreGroupPostProcessDir = f'{self.rescoreGroupFDRDir}/post_processing'
        self.rescoreMPGroupdir = f'{self.rescoreGroupPostProcessDir}/mp'
        self.rescoreAnnoGroupDir = f'{self.rescoreGroupPostProcessDir}/anno'
        self.mpPinFile = f'{self.rescoreGroupFDRDir}/mp.pin'
        self.annoPinFile = f'{self.rescoreGroupFDRDir}/anno.pin'

        self.cascadeDir = f'{self.outdir}/cascade'
        self.cascadeZeroPassDir = f'{self.cascadeDir}/zero_pass'
        self.cascadeZeroPassMzmlDir = f'{self.cascadeZeroPassDir}/mzml_files'
        self.cascadeFirstPassDir = f'{self.cascadeDir}/first_pass'
        self.cascadeSecondPassDir = f'{self.cascadeDir}/second_pass'
        self.cascadeMzmlDir = f'{self.cascadeDir}/mzml_files'

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

        # quant
        self.quantificationDir = f'{self.outdir}/quantification'

        self.flashLFQInput = f'{self.quantificationDir}/flash_lfq_input.tsv'

        # transcriptomics
        self.rnaDir = f'{self.outdir}/transcriptomics'
        self.assemblyDir = f'{self.rnaDir}/assembly'
        self.interAssemblyDir = f'{self.assemblyDir}/intermediate_assemblies'
        self.rnaAlnDir = f'{self.rnaDir}/alignments'
        self.rnaTrimmedDir = f'{self.rnaDir}/trimmed'
        self.rnaQCDir = f'{self.rnaDir}/QC'
        self.rnaNoContDir = f'{self.rnaTrimmedDir}/no_contaminants'
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


        # PGC
        self.pgContextDir = f'{self.outdir}/pg_context'
        self.intermediatePGCFiles = f'{self.pgContextDir}/intermediate_files'
        self.contextFiguresDir = f'{self.pgContextDir}/context_figures'

        # protein groups
        self.proteinGroupsDir = f'{self.outdir}/protein_groups'
        self.proteinGroups = f'{self.proteinGroupsDir}/protein_groups.csv'

        # shiny
        self.shinyRScript = f'{sys.path[0]}/src/shiny/shinyapp.R'

        # wgs
        self.wgsDir = f'{self.outdir}/WGS'
        self.variantsDir = f'{self.wgsDir}/variants'
        self.wgsTrimmeddir = f'{self.wgsDir}/trimmed_reads'
        self.wgsAlignDir = f'{self.wgsDir}/alignments'
        self.wgsSamDir = f'{self.wgsAlignDir}/sam'
        self.wgsBamDir = f'{self.wgsAlignDir}/bam'
        self.wgsDeduplicatedBamDir = f'{self.wgsAlignDir}/deduplicated_bam'
        self.wgsRecalDir = f'{self.variantsDir}/recalibration'
        self.wgsMutectDir = f'{self.variantsDir}/mutect2'

    def check_wgs_dirs(self):
        self.check_dirs([self.wgsDir, self.variantsDir, self.wgsTrimmeddir,
                         self.wgsAlignDir, self.wgsSamDir, self.wgsBamDir, self.wgsDeduplicatedBamDir,
                         self.wgsRecalDir, self.wgsMutectDir])

    # def __set__postms_mode(self):
    #     self.percInputSingle = f'{db_path}/percolator_input_single'

    def __define_genomes(self):
        files = {'hg38': {'genome': f'{self.refAnnoDir}/hg38/gencode.GRCh38.primary_assembly.genome.shortHead.fa',
                          'gtf': f'{self.refAnnoDir}/hg38/gencode.v38.annotation.gtf',
                          'ensembl_gtf': f'{self.refAnnoDir}/hg38/ensembl_hg38_chromRenamed.gtf',
                          'genome_index': f'{self.refAnnoDir}/hg38/STARindex_GencodeGTF',
                          'cont_index': f'{self.refAnnoDir}/hg38/STARindex_RNAcont'},
                'mm39': {'genome': f'{self.refAnnoDir}/mm39/mm39.fa',
                         'gtf': f'{self.refAnnoDir}/mm39/mm39.knownGene.gtf',
                         'proteome': f'{self.refAnnoDir}/mm39/mm_swissProt-trembl-Isoforms_uniprotkb_proteome_UP000000589_2025_03_09.fasta',
                         'ensembl_gtf': f'{self.refAnnoDir}/mm39/ensembl_mm39.gtf',
                         'genome_index': f'{self.refAnnoDir}/mm39/STARindex_GencodeGTF',
                         'cont_index': f'{self.refAnnoDir}/mm39/STARindex_RNAcont'},
        }
        if self.args.genomeAssembly in files:
            self.args.genome = files[self.args.genomeAssembly].get('genome', None)
            self.args.gtf = files[self.args.genomeAssembly].get('gtf', None)
            self.args.ensembleGTF = files[self.args.genomeAssembly].get('ensembl_gtf', None)
            self.args.genome_index = files[self.args.genomeAssembly].get('genome_index', None)
            self.args.cont_index = files[self.args.genomeAssembly].get('cont_index', None)
            if getattr(self.args, 'proteome', None) is None:
                self.args.proteome = files[self.args.genomeAssembly].get('proteome', None)
            if self.args.mode == 'anno':
                self.args.gtf = files[self.args.genomeAssembly]['ensembl_gtf']
            if self.args.mode == 'wgs':
                variant_data = f'{self.refAnnoDir}/hg38/variant'
                self.args.knownSites = [f'{variant_data}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf',
                                        f'{variant_data}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz']
                self.args.germlineResource = f'{variant_data}/hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz'
                self.args.panelOfNormals = f'{variant_data}/somatic-hg38_1000g_pon.hg38.vcf.gz'


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

    def set_demo_args(self):
        # rna
        # self.args.pairedPattern = ''
        # database

        self.args.genome = f'{self.testDir}/genome.fasta'
        self.args.proteome = f'{self.testDir}/proteome.fasta'
        self.args.refseq = f'{self.testDir}/refseq.fasta'
        self.args.gtf = f'{self.testDir}/demo.gtf'
        self.args.mzml = f'{self.testDir}/search/group'
        # self.args.outdir = f'{sys.path[0]}/demo_outdir/'
        self.outdir = self.args.outdir
        self.args.gtf_folder = f'{self.testDir}/database'
        self.args.external_database = None
        self.args.external_gtf = None
        self.args.uniprotAnnotation = None
        self.args.highHomologyDB = False


        # search
        self.args.quantifyOnly = False
        self.args.quantify = False
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
        self.args.keepAnnotated = False
        self.args.maxORFLength = 150
        self.args.phosphorylation = False
        self.args.groupedFDR = False

        # ribocov
        self.args.manualPlot = False
        self.args.adapter = "AGATCGGAAGAGCACACGTCT"
        self.args.genome_index = f'{sys.path[0]}/STAR_indexes/hg19.star'
        self.args.cont_index = f'{sys.path[0]}/STAR_indexes/hg19cont.star'
        self.args.aln = None
        self.args.fastq = f'{sys.path[0]}/demo_data/ribocov'
        self.args.fastx_clipper_path = f"{sys.path[0]}/dependencies/fastx_clipper"
        self.args.fastx_trimmer_path = f"{sys.path[0]}/dependencies/fastx_trimmer"
        self.args.multimappings = 9999
        self.args.rpkm = 0
        self.args.plots = True
        self.args.minRawCounts = 10
        self.args.grouping_method = 'union'

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
        if hasattr(self.args, 'externalFasta'):
            if self.args.externalFasta:
                fasta = self.args.externalFasta
        return fasta

    def split_big_fasta(self, fasta, sequences=10000):
        from Bio import SeqIO
        self.check_dirs([self.splitFastaDir])
        records = SeqIO.parse(fasta, 'fasta')
        
        for i, record in enumerate(records):
            if i % sequences == 0:
                if i != 0:
                    out_fasta.close()
                out_fasta = open(f'{self.splitFastaDir}/split_{i // sequences}.fasta', 'w')
            SeqIO.write(record, out_fasta, 'fasta')
            

    def select_psm_df(self):
        if os.path.exists(self.rescoredMicroproteinsFasta):
            psm = f'{self.rescorePostProcessDir}/group/psm_fixed.txt'
        else:
            psm = f'{self.postProcessDir}/group/db/psm_fixed.txt'
        return psm

    def select_peptides_df(self):
        if os.path.exists(self.rescoredMicroproteinsFasta):
            pep = f'{self.rescorePostProcessDir}/group/peptides_fixed.txt'
        else:
            pep = f'{self.postProcessDir}/group/db/peptides_fixed.txt'
        return pep
    
    def select_database(self, decoy=False, proteome=False, split_db=False):
        """
        Returns full path to the proper database. It will check if rescored. 
        If rescored, it will return the rescored database. If not, it will return the original database.
        if split_db, it will return a list with the paths to every db in the split databases folder."""

        if split_db:
            dbs = []
            if proteome:
                files = os.listdir(self.splitDbProteomeDir)
                for file in files:
                    if not file.endswith(".idx"):
                        dbs.append(os.path.join(self.splitDbProteomeDir, file))
            else:
                files = os.listdir(self.splitDbProteogenomicsDir)
                for file in files:
                    if not file.endswith(".idx"):
                        dbs.append(os.path.join(self.splitDbProteogenomicsDir, file))
            return dbs

        if proteome:
            # if self.args.splitDatabase is not None:
            db = self.refProteome
            if decoy:
                db = self.refProteomeWithDecoy
            return db
        
        if os.path.exists(self.rescoredMicroproteinsFasta):
            db = f'{self.rescoreDatabaseDir}/rescore_target_database.fasta'
            if decoy:
                db = f'{self.rescoreDatabaseDir}/rescore_target_decoy_database.fasta'

        else:
            dbs = os.listdir(self.databaseDir)
            db = None
            for file in dbs:
                if decoy:
                    if file.endswith("target_decoy_database.fasta"):
                        db = f'{self.databaseDir}/{file}'
                        break
                else:
                    if file.endswith("target_database.fasta"):
                        db = f'{self.databaseDir}/{file}'
                        break
        return db


    def is_rescored(self):
        if os.path.exists(self.rescoredMicroproteinsFasta):
            rescored = True
        else:
            rescored = False
        return rescored
    
    def select_search_dir(self):
        if os.path.exists(self.rescoreSearchDir):
            return self.rescoreSearchDir
        else:
            return self.searchDir

    def select_gtf(self):
        if os.path.exists(self.rescoredMicroproteinsGTF):
            gtf = self.rescoredMicroproteinsGTF
        else:
            gtf = self.uniqueMicroproteinsGTF
        return gtf

    # def select_database(self, decoy=False):
    #     db = None
    #     if self.is_rescored():
    #         if not decoy:
    #             db = f'{self.rescoreDatabaseDir}/rescore_target_database.fasta'
    #         else:
    #             db = f'{self.rescoreDatabaseDir}/rescore_target_decoy_database.fasta'
    #     else:
    #         dbs = os.listdir(self.databaseDir)
    #         for fasta in dbs:
    #             if decoy:
    #                 if fasta.endswith("target_decoy_database.fasta"):
    #                     db = fasta
    #             else:
    #                 if fasta.endswith("target_database.fasta"):
    #                     db = fasta
    #     return db


    def verify_checkpoint(self, outfile, step, mute=False):
        if os.path.exists(outfile) and not self.args.overwrite:
            if not mute:
                print(f"(!) Found output file {outfile}. Skipping {step}...")
            run = False
        else:
            run = True
        return run

    def print_row(self, n=80, word='', character='='):
        char_n = (n-(len(word))) / 2
        first = character * int(char_n)
        row = f'{first}{word}{first}'
        print(row)

    def get_microprotein_mapping_groups(self):
        mapping_groups = {}
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            mapping_groups[smorf] = group
        return mapping_groups

    def associate_groups_to_files(self):
        groups_dict = {}
        if self.args.groupsFile is not None:
            df = pd.read_csv(self.args.groupsFile, sep='\t')
            groups = df["group"].tolist()
            files = df["file"].tolist()
            groups_dict = {}
            for group, file in zip(groups, files):
                groups_dict[file.replace(".mzML", "")] = group
        return groups_dict


class Content:
    def __init__(self, file, fullfile, group, db, main_dir):
        self.mainDir = main_dir

        self.group = group
        self.groupDir = f'{main_dir}/{group}'

        self.db = db
        self.dbDir = f'{self.groupDir}/{db}'

        self.file = file
        self.fullFile = fullfile


