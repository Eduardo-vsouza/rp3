import os
import sys
from ..utils import check_multiple_dirs




class PipelineStructure:
    def __init__(self, args):
        self.args = args

        self.outdir = args.outdir
        self.translationDir = f'{self.outdir}/translation'
        self.databaseDir = f'{self.outdir}/databases'
        self.resultsDir = f'{self.outdir}/results'
        self.searchDir = f'{self.outdir}/peptide_search'
        self.postProcessDir = f'{self.outdir}/post_processing'
        self.quantDir = f'{self.outdir}/quantification'
        self.summarizedResultsDir = f'{self.outdir}/summarized_results'
        self.mergedResults = f'{self.summarizedResultsDir}/merged'
        self.spectrumDir = f'{self.outdir}/spectra'
        self.signalPDir = f'{self.mergedResults}/signalP'

        self.pipelineDir = f'{sys.path[0]}'
        self.indexesDir = f'{self.pipelineDir}/data/STAR_indexes'

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


        # files
        self.uniqueMicroproteins = f'{self.mergedResults}/microproteins_150.fasta'
        self.utpsMicroproteins = f'{self.mergedResults}/microproteins_utps_150.fasta'
        self.utpsMicroproteinsBlast = f'{self.mergedResults}/microproteins_utps_150_blast_filt.fasta'
        self.microproteinsBlast = f'{self.mergedResults}/microproteins_150.fasta_blast_filt.fasta'
        self.metricsDir = f'{self.outdir}/metrics'
        self.rescoredMicroproteinsFasta = f'{self.rescoreSummarizedResultsDir}/filtered_rescored_microproteins_150.fasta'
        self.rescoredMicroproteinsGTF = f'{self.rescoreSummarizedResultsDir}/filtered_rescored_microproteins_150.gtf'
        self.mergedFullGTF = f'{self.mergedResults}/merged_predicted_microproteins.gtf'


        # ribo-seq
        self.riboSeqDir = f'{self.outdir}/ribo-seq'
        self.riboSeqTrimmedDir = f'{self.riboSeqDir}/trimmed_reads'
        self.riboSeqContaminantAlnDir = f'{self.riboSeqDir}/contaminant_alns'
        self.riboSeqAlnDir = f'{self.riboSeqDir}/alignments'

        # conservation
        self.phyloDir = f'{self.outdir}/phylogenetics'
        self.blastDir = f'{self.phyloDir}/blast'

        self.configFile = f'{self.pipelineDir}/config.txt'

        self.toolPaths = {}

    def read_config_file(self):
        with open(self.configFile, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                splat = line.split("\t")
                tool = splat[0]
                path = splat[1]
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


class Content:
    def __init__(self, file, fullfile, group, db, main_dir):
        self.mainDir = main_dir

        self.group = group
        self.groupDir = f'{main_dir}/{group}'

        self.db = db
        self.dbDir = f'{self.groupDir}/{db}'

        self.file = file
        self.fullFile = fullfile


