import os
import sys

from src.pipeline import Pipeline
from ..pipeline_config import PipelineStructure


class Demo(PipelineStructure):
    def __init__(self, args):

        super().__init__(args=args)

        self.testedModes = {'translation': None, 'database': None, 'search': None, 'postms': None,
                            'ribocov': None}
        self.issues = {}

    # self.check_dirs([self.demoO])
    def test(self):
        self.__test_database_mode()
        self.__test_search_mode()
        self.__test_ribocov_mode()
        # self.__test_rna_mode()

    def __test_rna_mode(self):
        if not self.args.skip_rna:
            pipe = Pipeline(args=self.args)
            pipe.assemble_transcriptomes()
            skip = False
        else:
            skip = True
        # self.__check_files(file=self.spliceSites, mode='rna', tool='hisat2_extract_splice_sites.py', skipped=skip)
        # self.__check_files(file=self.exons, mode='rna', tool='hisat2_extract_exons.py')
        # self.__check_files(file=f'{self.args.index}/hisat_genome_index.8.ht2', mode='rna', tool='hisat2-build')
        self.__check_files(file=f'{self.assemblyDir}/group/merged_assembly.gtf', mode='rna', tool='stringtie',
                           skipped=skip)


    def __test_database_mode(self):
        if not self.args.skip_database:
            skip = False
            pipe = Pipeline(args=self.args)

            pipe.genome = self.args.genome
            pipe.proteome = self.args.proteome

            pipe.translate_in_silico()
            pipe.generate_databases()
        else:
            skip = True
        # if not skip:
        translation_outfiles = ['GSE125218_K562_Cufflinks_transcript_assembly.nuc',
                                'GSE125218_K562_Cufflinks_transcript_assembly_ORFs.gtf',
                                'GSE125218_K562_Cufflinks_transcript_assembly.pep',
                                'GSE125218_K562_Cufflinks_transcript_assembly.split_nuc']
        database_outfiles = ['GSE125218_K562_Cufflinks_transcript_assembly_target_decoy_database.fasta']
        for file in translation_outfiles:
            fullfile = f'{self.translationDir}/GSE125218_K562_Cufflinks_transcript_assembly/{file}'
            self.__check_files(file=fullfile, mode="translation", tool="GTFtoFasta", skipped=skip,
                               message="Additionally, check if the provided GTF files are correct. If "
                                       "the pipeline was used to assemble the transcriptome, check if it was "
                                       "generated correctly during the 'rna' mode.\n")
            # if os.path.exists(fullfile):
            #     if os.path.getsize(fullfile) > 1000:
            #         self.testedModes['translation'] = "OK"
            #     else:
            #         self.testedModes['translation'] = "Error"
        for file in database_outfiles:
            fullfile = f'{self.databaseDir}/{file}'
            self.__check_files(file=fullfile, mode='database', tool="GTFtoFasta", skipped=skip,
                               message="Check if the provided GTF and proteome fasta files are correct.\n")
                # if os.path.exists(fullfile):
                #     if os.path.getsize(fullfile) > 1000:
                #         self.testedModes['database'] = "OK"
                #     else:
                #         self.testedModes['database'] = "Error"
        # else:
        #     self.testedModes['translation'] = 'Skipped'
        #     self.testedModes['database'] = 'Skipped'

    def __test_search_mode(self):
        pipe = Pipeline(args=self.args)
        if not self.args.skip_search:
            skip = False
            pipe.search_mass_spec()
        else:
            skip = True
        file = f'{self.searchDir}/group/GSE125218_K562_Cufflinks_transcript_assembly_target_decoy_database.fasta/demo_target.pin'
        self.__check_files(file=file, mode='search', skipped=skip, tool="MSFragger", check_size=False)

        if not self.args.skip_postms:
            skip = False
            pipe.post_process_searches()
            pipe.rescore()
        else:
            skip = True
        file = f'{self.postProcessDir}/group/db/peptides_filtered.txt'
        self.__check_files(file=file, mode='postms', skipped=skip, tool='Percolator')
        rescore_file = f'{self.rescoreSummarizedResultsDir}/filtered_rescored_microproteins_150.fasta'
        self.__check_files(file=rescore_file, mode="postms", skipped=skip, tool='MSFragger and Percolator')


    def __test_ribocov_mode(self):
        pipe = Pipeline(args=self.args)
        if not self.args.skip_ribocov:
            skip = False
            pipe.check_riboseq_coverage()
        else:
            skip = True
        trimmed_file = f'{self.riboSeqTrimmedDir}/trimmed_SRR8449580.fastq'
        self.__check_files(file=trimmed_file, mode='ribocov', tool="FastX Toolkit", skipped=skip)
        aln_file = f'{self.riboSeqAlnDir}/aligned_to_genome_no_contaminant_trimmed_SRR8449580.fastqUnmapped.out.mate1Aligned.out.sam'
        self.__check_files(file=aln_file, mode='ribocov', tool='STAR', skipped=skip)
        raw_counts_file = f'{self.countsDir}/raw/rescored_smorfs_plus_reference_annotation_ambiguous_counts.txt'
        self.__check_files(file=raw_counts_file, mode='ribocov', tool='featureCounts', skipped=skip,
                           message='Additionally, featureCounts supports a max of 64 threads. Make sure you have not '
                                   'specified a number of threads exceeding the limit.')
        self.__check_files(file=f'{self.countsDir}/heatmap_.png', tool='Python package nheatmap', skipped=skip,
                           mode='ribocov', check_size=False)

    def __check_files(self, file, mode, tool, skipped=False, message=None, check_size=True):
        if self.testedModes[mode] != 'Failed':
            # if not skipped:
            if os.path.exists(file):
                if check_size:
                    if os.path.getsize(file) > 100:
                        self.testedModes[mode] = "OK"
                    else:
                        self.testedModes[mode] = "Failed"
                else:
                    self.testedModes[mode] = 'OK'
            else:
                self.testedModes[mode] = "Failed"
            # else:
            #     self.testedModes[mode] = 'Skipped'
            if self.testedModes[mode] == 'Failed':
                if mode not in self.issues:
                    self.issues[mode] = []
                self.issues[mode].append(f"{file} was not generated properly during '{mode}' mode."
                                         f" Please check your {tool} installation.\n")
                if tool in self.toolPaths:
                    self.issues[mode].append(f"The provided path to this tool is {self.toolPaths[tool]}.\n")
                if message is not None:
                    self.issues[mode].append(message)

    def print_status(self):
        print("\nWelcome to demo mode.\nWe are currently testing which Rp3 modes are working properly.")
        status = (f'\n{"="*11}Modes{"="*11}\n'
                  # f' rna \t{" "*9}   {self.testedModes["rna"]}\n'
                  f' translation\t    {self.testedModes["translation"]}\n'
                  f' database\t    {self.testedModes["database"]} \n'
                  f' search{" "*9}    {self.testedModes["search"]}\n'
                  f' postms{" "*9}    {self.testedModes["postms"]}\n'
                  f' ribocov{" "*8}    {self.testedModes["ribocov"]}\n'
                  f'{"="*27}\n')
        print(status)

        for mode in self.testedModes:
            if self.testedModes[mode] == 'Failed':
                print(f"- {mode} returned a 'Failed' status. Possible reasons are:")
                issue = self.issues[mode]
                for message in issue:
                    print(f"-- {message}")