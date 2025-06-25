import os
import sys
import re

import pandas as pd
from Bio import SeqIO

from .utils import group_folder_generator, check_multiple_dirs
from .pipeline_config import PipelineStructure

import concurrent.futures

from tqdm import tqdm


class PercolatorPostProcessing(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.args = args
        self.threads = args.threads

        # directories
        self.outdir = args.outdir
        self.translationDir = f'{self.outdir}/translation'
        self.peptideSearchDir = f'{self.outdir}/peptide_search'
        self.postProcessDir = f'{self.outdir}/post_processing'
        self.resultsDir = f'{self.outdir}/results'
        self.summarizedResultsDir = f'{self.outdir}/summarized_results'
        self.mergedDir = f'{self.summarizedResultsDir}/merged'
        self.__check_dirs([self.postProcessDir, self.resultsDir, self.summarizedResultsDir, self.mergedDir])

        self.mode = 'postms'
        self.params = []

    @staticmethod
    def __check_dirs(folders):
        if type(folders) == list:
            for folder in folders:
                if not os.path.exists(folder):
                    os.mkdir(folder)
        else:
            if not os.path.exists(folders):
                os.mkdir(folders)
        # if not os.path.exists(self.postProcessDir):
        #     os.mkdir(self.postProcessDir)
        # if not os.path.exists(self.resultsDir):
        #     os.mkdir(self.resultsDir)
        # if not os.path.exists(self.summarizedResultsDir):
        #     os.mkdir(self.summarizedResultsDir)

    #

    # deprecated
    def convert_to_pin(self):
        groups = os.listdir(self.peptideSearchDir)
        for group in groups:
            dbs = os.listdir(f'{self.peptideSearchDir}/{group}')
            files_to_convert = {}
            for db in dbs:

                if db.endswith(".fasta"):

                    # GSE125218_HeLa-S3_Cufflinks_transcript_assembly_target_database.fasta
                    # would look like
                    # GSE125218_HeLa-S3_Cufflinks_transcript_assembly
                    db_prefix = '_'.join(db.split(".fasta")[0].split("_")[:-2])

                    if db_prefix not in files_to_convert:
                        files_to_convert[db_prefix] = []  # this is to separate different DBs for a same group
                        # for instance, multiple transcriptome assemblies but the same mass spec data

                    pin_dir = f'{self.peptideSearchDir}/{group}/pin_files'
                    if not os.path.exists(pin_dir):
                        os.mkdir(pin_dir)

                    pin_dir_db = f'{pin_dir}/{db_prefix}'  # the pin folder should contain a different folder
                    if not os.path.exists(pin_dir_db):  # for each transcriptome assembly
                        os.mkdir(pin_dir_db)

                    db_folder = f'{self.peptideSearchDir}/{group}/{db}'
                    protein_db_path = f'{self.outdir}/databases/{db}'
                    pepxml_files = os.listdir(f'{self.peptideSearchDir}/{group}/{db}')

                    for pepxml in pepxml_files:
                        # if pepxml.endswith(".pepXML"):
                        #     cmd_convert = f'/home/eduardo/programs/crux/crux-4.1.Linux.x86_64/bin/crux psm-convert ' \
                        #                   f'--output-dir {pin_dir} --protein-database {protein_db_path} --overwrite T' \
                        #                   f' --input-format pepxml {self.peptideSearchDir}/{group}/{db}/{pepxml} ' \
                        #                   f'pin'
                        if pepxml.endswith(".pepXML"):
                            files_to_convert[db_prefix].append(f' {db_folder}/{pepxml}')
            # print(files_to_convert)
            # for db in files_to_convert:
            #     print("\n")
            #     print(db)
            #     for file in files_to_convert[db]:
            #         print(file)

            for db in files_to_convert:
                cmd_convert = f'/home/eduardo/programs/crux/crux-4.1.Linux.x86_64/bin/crux make-pin ' \
                              f'--overwrite T --output-file {self.peptideSearchDir}/{group}/pin_files/{group}_{db}.pin --decoy_prefix rev_ ' \
                              f'{" ".join(files_to_convert[db])}'
                os.system(cmd_convert)

    def merge_pin_replicates(self):
        groups = os.listdir(self.peptideSearchDir)
        for group in groups:
            group_dir = f'{self.peptideSearchDir}/{group}'
            dbs = os.listdir(group_dir)
            for db_dir in dbs:
                if 'target' in db_dir and db_dir.endswith(".fasta"):
                    db_path = f'{group_dir}/{db_dir}'  # the directories are named after the database that was used
                    pin_files = os.listdir(db_path)  # for the search
                    merged_pin = []
                    for pin in pin_files:
                        if pin.endswith(".pin"):
                            with open(f'{db_path}/{pin}', 'r') as handler:
                                lines = handler.readlines()
                                for i, line in enumerate(lines):
                                    if len(merged_pin) < 1:
                                        if i == 0:
                                            merged_pin.append(line)
                                    else:
                                        if i > 0:
                                            # if 'decoy' in pin:
                                            #     cols = line.split("\t")

                                            merged_pin.append(line)
                    percolator_input = f'{db_path}/percolator_input'
                    if not os.path.exists(percolator_input):
                        os.mkdir(percolator_input)
                    with open(f'{percolator_input}/{db_dir}.pin', 'w') as outfile:
                        outfile.writelines(merged_pin)

    def merge_all_pins(self):
        groups = os.listdir(self.peptideSearchDir)
        merged_pin = []
        for group in groups:
            group_dir = f'{self.peptideSearchDir}/{group}'
            if os.path.isdir(group_dir):
                dbs = os.listdir(group_dir)
                for db_dir in dbs:
                    if 'target_decoy' in db_dir and db_dir.endswith(".fasta"):
                        db_path = f'{group_dir}/{db_dir}'  # the directories are named after the database that was used
                        pin_files = os.listdir(db_path)  # for the search
                        for pin in pin_files:
                            if pin.endswith(".pin"):
                                with open(f'{db_path}/{pin}', 'r') as handler:
                                    lines = handler.readlines()
                                    for i, line in enumerate(lines):
                                        if len(merged_pin) < 1:
                                            if i == 0:
                                                merged_pin.append(line)
                                        else:
                                            if i > 0:
                                                merged_pin.append(line)
        with open(f'{self.outdir}/all_pins.pin', 'w') as out:
            out.writelines(merged_pin)

    def percolate(self):
        groups = os.listdir(self.peptideSearchDir)
        for group in groups:
            group_dir = f'{self.peptideSearchDir}/{group}'
            if os.path.isdir(group_dir):
                dbs = os.listdir(group_dir)
                for db_dir in dbs:
                    if 'target' in db_dir and db_dir.endswith(".fasta"):
                        db_path = f'{group_dir}/{db_dir}'  # the directories are named after the database that was used
                        percolator_input = f'{db_path}/percolator_input'
                        pin = f'{percolator_input}/{db_dir}.pin'

                        group_outdir = f'{self.postProcessDir}/{group}'
                        if not os.path.exists(group_outdir):
                            os.mkdir(group_outdir)
                        outdir = f'{self.postProcessDir}/{group}/{db_dir}'
                        if not os.path.exists(outdir):
                            os.mkdir(outdir)

                        cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                                         f'--post-processing-tdc --results-psms {outdir}/psm.txt --results-peptides ' \
                                         f'{outdir}/peptides.txt --no-terminate --num-threads {self.threads} ' \
                                         f'-X {outdir}/pout.xml {pin}'
                        # cmd_percolator = f'percolator --protein-decoy-pattern rev_ ' \
                        #                  f'--post-processing-mix-max --results-psms {outdir}/psm.txt --results-peptides ' \
                        #                  f'{outdir}/peptides.txt -f auto --test-each-iteration -I separate  ' \
                        #                  f'-X {outdir}/pout.xml {pin}'
                        os.system(cmd_percolator)

    # def percolate_single(self):
    #     """
    #     On 'sep' mode, this runs percolator on single .pin files, one by one, instead of using a merged pin file
    #      containing all of these from a group. The intermediate output files generated by percolator will be stored in
    #     percolator_output_single, still inside the group_dir for the search. Afterwards, these pin files are
    #     merged into a single one for the cutoffs to be applied, located in the self.postProcessDir'
    #     """
    #     groups = os.listdir(self.peptideSearchDir)
    #     for group in groups:
    #         group_dir = f'{self.peptideSearchDir}/{group}'
    #         if os.path.isdir(group_dir):
    #             dbs = os.listdir(group_dir)
    #             for db_dir in dbs:
    #                 db_path = f'{group_dir}/{db_dir}'  # the directories are named after the database that was used
    #                 if 'target_decoy' in db_dir and db_dir.endswith(".fasta"):
    #                     pin_files = os.listdir(f'{db_path}')
    #                     outdir = f'{db_path}/percolator_output_single'
    #                     self.__check_dirs([outdir])
    #                     for file in pin_files:
    #                         if file.endswith("target.pin"):
    #                             fullfile = f'{db_path}/{file}'
    #                             outfile = file.replace(".pin", "")
    #                             cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
    #                                              f'--post-processing-tdc --results-psms {outdir}/{outfile}_psm.txt --results-peptides ' \
    #                                              f'{outdir}/{outfile}_peptides.txt --no-terminate --num-threads {self.threads} ' \
    #                                              f'-X {outdir}/{outfile}_pout.xml {fullfile}'
    #                             os.system(cmd_percolator)
    #                     group_outdir = f'{self.postProcessDir}/{group}'
    #                     if not os.path.exists(group_outdir):
    #                         os.mkdir(group_outdir)
    #                     outdir = f'{self.postProcessDir}/{group}/{db_dir}'
    #                     if not os.path.exists(outdir):
    #                         os.mkdir(outdir)

    def percolate_cascade(self):
        # first pass FDR assessment
        db = self.select_database(decoy=True, proteome=True)
        pin = f'{self.searchDir}/cascade_first_pass_final.pin'
        self.print_state(message=f"Running Percolator on {pin} with {db}", color='yellow')
        cmd = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                f'--post-processing-tdc --results-psms {self.cascadePostmsDir}/cascade_first_pass_psm.txt ' \
                f'--results-peptides {self.cascadePostmsDir}/cascade_first_pass_peptides.txt --no-terminate ' \
                f'--num-threads {self.args.threads} -X {self.cascadePostmsDir}/cascade_first_pass_pout.xml ' \
                f'--picked-protein {db} --results-proteins {self.cascadePostmsDir}/cascade_first_pass_proteins.txt {pin}'
        os.system(cmd)
        self.print_state(message=f"Percolator completed for {pin}", color='green')

        # second pass FDR assessment
        db = self.select_database(decoy=True, proteome=False)
        pin = f'{self.searchDir}/cascade_second_pass_final.pin'
        self.print_state(message=f"Running Percolator on {pin} with {db}", color='yellow')
        cmd = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                f'--post-processing-tdc --results-psms {self.cascadePostmsDir}/cascade_second_pass_psm.txt ' \
                f'--results-peptides {self.cascadePostmsDir}/cascade_second_pass_peptides.txt --no-terminate ' \
                f'--num-threads {self.args.threads} -X {self.cascadePostmsDir}/cascade_second_pass_pout.xml ' \
                f'--picked-protein {db} --results-proteins {self.cascadePostmsDir}/cascade_second_pass_proteins.txt {pin}'
        os.system(cmd)
        self.print_state(message=f"Percolator completed for {pin}", color='green')

        # concatenate results for filtering
        group_outdir = f'{self.postProcessDir}/group/db'
        cmd = f'cat {self.cascadePostmsDir}/*_peptides.txt > {group_outdir}/peptides.txt'
        os.system(cmd)
        cmd = f'cat {self.cascadePostmsDir}/*_psm.txt > {group_outdir}/psm.txt'
        os.system(cmd)
        cmd = f'cat {self.cascadePostmsDir}/*_proteins.txt > {group_outdir}/proteins.txt'
        os.system(cmd)

        # cmd = (f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ '
        #         f'--post-processing-tdc --results-psms {group_outdir}/{pin}_psm.txt --results-peptides '
        #         f'{group_outdir}/{pin}_peptides.txt --no-terminate --num-threads {self.args.threads} '
        #         f'-X {group_outdir}/{pin}_pout.xml --picked-protein {self.databaseDir}/{dbss} '
        #         f'--results-proteins {group_outdir}/{pin}_proteins.txt '
        #         f'{self.searchDir}/group/{db}/{pin}')

    def percolate_single(self):
        databases = os.listdir(self.databaseDir)

        for db in databases:
            if db.endswith("target_database.fasta"):
                dbss = db

        group_outdir = f'{self.postProcessDir}/group/db'
        self.check_dirs([f'{self.postProcessDir}/group', group_outdir])
        dbs = os.listdir(f'{self.searchDir}/group')
        for db in dbs:
            if db.endswith("target_decoy_database.fasta"):
                pins = os.listdir(f'{self.searchDir}/group/{db}')
                for pin in pins:
                    run = False
                    if self.args.msBooster:
                        if '_edited' in pin:
                            run = True
                    else:
                        if pin.endswith(".pin"):
                            run = True
                    if run:
                        cmd = (f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ '
                               f'--post-processing-tdc --results-psms {group_outdir}/{pin}_psm.txt --results-peptides '
                               f'{group_outdir}/{pin}_peptides.txt --no-terminate --num-threads {self.args.threads} '
                               f'-X {group_outdir}/{pin}_pout.xml --picked-protein {self.databaseDir}/{dbss} '
                               f'--results-proteins {group_outdir}/{pin}_proteins.txt '
                               f'{self.searchDir}/group/{db}/{pin}')
                        os.system(cmd)
        cmd = f'cat {group_outdir}/*_peptides.txt > {group_outdir}/peptides.txt'
        os.system(cmd)
        cmd = f'cat {group_outdir}/*_psm.txt > {group_outdir}/psm.txt'
        os.system(cmd)
        cmd = f'cat {group_outdir}/*_proteins.txt > {group_outdir}/proteins.txt'
        os.system(cmd)

    def percolate_all_pins(self):
        pin = f'{self.outdir}/all_pins.pin'
        if self.args.cascade:
            pin = os.path.join(self.peptideSearchDir, 'cascade_search.pin')
        # self.flatten_orf_duplicates_pin(pin=pin, output=f'{self.outdir}/all_pins_flattened.pin')
        g_outdir = f'{self.postProcessDir}/group'
        dbs = os.listdir(self.databaseDir)
        for db in dbs:
            if db.endswith("target_decoy_database.fasta"):
                outdir = f'{g_outdir}/db'
                self.check_dirs([g_outdir, outdir])
        databases = os.listdir(self.databaseDir)
        for db in databases:
            if db.endswith("target_database.fasta"):
                dbss = db
        # print(dbss)
        best_positive = ''
        if self.args.splitDatabase is not None:
            best_positive = '--train-best-positive'
        enzyme = self.args.enzyme
        
        if self.args.hlaPeptidomics:
            enzyme = 'no_enzyme'
        cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                         f'--post-processing-tdc --results-psms {outdir}/_psm.txt --results-peptides ' \
                         f'{outdir}/_peptides.txt --protein-enzyme {enzyme} --no-terminate --picked-protein {self.databaseDir}/{dbss} --results-proteins' \
                         f' {outdir}/proteins.txt --num-threads {self.threads} {best_positive} ' \
                         f'-X {outdir}/pout.xml {pin}'
        os.system(cmd_percolator)

    def percolate_multi(self):
        groups = os.listdir(self.peptideSearchDir)
        for group in groups:
            group_dir = f'{self.peptideSearchDir}/{group}'
            dbs = os.listdir(group_dir)
            for db_dir in dbs:
                db_path = f'{group_dir}/{db_dir}'  # the directories are named after the database that was used

                if 'target' in db_dir and db_dir.endswith(".fasta"):
                    pin = ''
                    pin_files = os.listdir(f'{db_path}')
                    for file in pin_files:

                        if file.endswith("pin"):
                            pin += f' {db_path}/{file}'
                            # percolator_input = f'{db_path}/percolator_input'
                            # pin = f'{percolator_input}/{db_dir}.pin'
                        """
                            group_outdir = f'{self.postProcessDir}/{group}'
                            if not os.path.exists(group_outdir):
                                os.mkdir(group_outdir)
                            outdir = f'{self.postProcessDir}/{group}/{db_dir}'
                            if not os.path.exists(outdir):
                                os.mkdir(outdir)

                            cmd_percolator = f'percolator --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                                             f'--post-processing-tdc --results-psms {outdir}/{file}_psm.txt --results-peptides ' \
                                             f'{outdir}/{file}_peptides.txt --no-terminate --num-threads {self.threads} ' \
                                             f'-X {outdir}/pout.xml {db_path}/{file}'
                            # cmd_percolator = f'percolator --protein-decoy-pattern rev_ ' \
                            #                  f'--post-processing-mix-max --results-psms {outdir}/psm.txt --results-peptides ' \
                            #                  f'{outdir}/peptides.txt -f auto --test-each-iteration -I separate  ' \
                            #                  f'-X {outdir}/pout.xml {pin}'
                            os.system(cmd_percolator)
                        """
                    # if file.endswith("pin"):
                    # pin += f' {db_path}/{file}'
                    # percolator_input = f'{db_path}/percolator_input'
                    # pin = f'{percolator_input}/{db_dir}.pin'

                    group_outdir = f'{self.postProcessDir}/{group}'
                    if not os.path.exists(group_outdir):
                        os.mkdir(group_outdir)
                    outdir = f'{self.postProcessDir}/{group}/{db_dir}'
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                    dbss = os.listdir(self.databaseDir)[0]

                    cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                                     f'--post-processing-tdc --results-psms {outdir}/{file}_psm.txt --results-peptides ' \
                                     f'{outdir}/{file}_peptides.txt --no-terminate --num-threads {self.threads} ' \
                                     f'-X {outdir}/pout.xml --picked-protein {dbss} --results-proteins' \
                                     f' {outdir}/proteins.txt {db_path}/percolator_input/{db_dir}.pin'
                    # cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                    #                  f'--post-processing-tdc --results-psms {outdir}/{pin}_psm.txt --results-peptides ' \
                    #                  f'{outdir}/{file}_peptides.txt --no-terminate --num-threads {self.threads} ' \
                    #                  f'-X {outdir}/pout.xml {pin}'
                    # cmd_percolator = f'percolator --protein-decoy-pattern rev_ ' \
                    #                  f'--post-processing-mix-max --results-psms {outdir}/psm.txt --results-peptides ' \
                    #                  f'{outdir}/peptides.txt -f auto --test-each-iteration -I separate  ' \
                    #                  f'-X {outdir}/pout.xml {pin}'
                    os.system(cmd_percolator)

    def fix_multiple_columns(self):
        # if self.args.postms_mode == 'sep':
        #     groups = os.listdir(self.searchDir)
        #     for group in groups:
        #         if os.path.isdir(f'{self.searchDir}/{group}'):
        #             dbs = os.listdir(f'{self.searchDir}/{group}')
        #             for db in dbs:
        #                 if db.endswith("target_decoy_database.fasta") or db == 'db':
        #                     postms_dir = f'{self.postProcessDir}/{group}/{db}'
        #
        #                     db_dir = f'{self.searchDir}/{group}/{db}/percolator_output_single'
        #                     cmd = f'cat {db_dir}/*_peptides.txt > {postms_dir}/peptides.txt'
        #                     os.system(cmd)
        #                     cmd = f'cat {db_dir}/*_psm.txt > {postms_dir}/psm.txt'
        #                     os.system(cmd)
        #                     self.fix(file=f'{postms_dir}/peptides.txt', output=f'{postms_dir}/peptides_fixed.txt')
        #                     self.fix(file=f'{postms_dir}/psm.txt', output=f'{postms_dir}/psm_fixed.txt')
        # else:
        groups = os.listdir(self.postProcessDir)

        for group in groups:
            dbs = os.listdir(f'{self.postProcessDir}/{group}')
            for db in dbs:
                db_dir = f'{self.postProcessDir}/{group}/{db}'
                if self.args.postms_mode == 'cat':
                    cmd = f'cat {db_dir}/*_peptides.txt > {db_dir}/peptides.txt'
                    os.system(cmd)
                    cmd = f'cat {db_dir}/*_psm.txt > {db_dir}/psm.txt'
                    os.system(cmd)
                self.fix(file=f'{db_dir}/peptides.txt', output=f'{db_dir}/peptides_fixed.txt')
                self.fix(file=f'{db_dir}/psm.txt', output=f'{db_dir}/psm_fixed.txt')

    @staticmethod
    def fix(file, output):
        fixed = []
        with open(file, 'r') as pep_handler:
            lines = pep_handler.readlines()
            for i, line in enumerate(lines):
                if i > 0:
                    if 'q-value' not in line:
                        line = line.rstrip()
                        splat = line.split("\t")
                        # print(splat)
                        new_line = ['\t'.join(splat[:5]), ','.join(splat[5:])]
                        fixed.append('\t'.join(new_line))
                        fixed.append("\n")
                else:
                    fixed.append(line)
        with open(output, 'w') as outfile:
            outfile.writelines(fixed)

    def remove_annotated(self, pattern='_ANNO'):
        # pattern = '_MOUSE'
        groups = os.listdir(self.postProcessDir)
        for group in groups:
            group_dir = f'{self.postProcessDir}/{group}'

            dbs = os.listdir(group_dir)
            for db in dbs:
                if db.endswith("target_decoy_database.fasta") or db == 'db':

                    db_dir = f'{group_dir}/{db}'
                    if not os.path.exists(f'{self.resultsDir}/{group}'):
                        os.mkdir(f'{self.resultsDir}/{group}')
                    results_dir = f'{self.resultsDir}/{group}/{db}'
                    if not os.path.exists(results_dir):
                        os.mkdir(results_dir)
                    peptide_df = f'{db_dir}/peptides_fixed.txt'
                    if os.path.getsize(peptide_df) > 0:

                        peptide_filtered = f'{db_dir}/peptides_filtered.txt'
                        self.filter_table(file=peptide_df, output=peptide_filtered, pattern=pattern)
                        #
                        psm_df = f'{db_dir}/psm_fixed.txt'
                        psm_filtered = f'{db_dir}/psm_filtered.txt'
                        self.filter_table(file=psm_df, output=psm_filtered, pattern=pattern)

                        # this generates a table containing only Unique Tryptic Peptides (UTPs), even for predicted MPs
                        peptide_filtered_utps_df = f'{db_dir}/peptides_filtered_utps.txt'
                        self.filter_table(file=peptide_df, output=peptide_filtered_utps_df, pattern=pattern,
                                          utps_only=True)

                        db_path = f'{self.outdir}/databases/{db}'
                        if db == 'db':
                            alldbs = os.listdir(self.databaseDir)
                            for d in alldbs:
                                if d.endswith("target_decoy_database.fasta"):
                                    db_path = f'{self.databaseDir}/{d}'
                                    break
                        fasta = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                                  database=db_path, results_dir=results_dir,
                                                  smorfs=True)

                        fasta_utps = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                                       database=db_path, results_dir=results_dir,
                                                       smorfs=True, utps=True)
                        fasta_proteins = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                                           database=db_path, results_dir=results_dir,
                                                           smorfs=False, utps=True)

                        assembly = db.replace("_target_decoy_database.fasta", "")
                        # if self.args.create_gtf == True:
                        #     self.create_gtf(microproteins_fasta=fasta, results_dir=results_dir,
                        #                     translation_gtf=f'{self.translationDir}/{assembly}/{assembly}_ORFs.gtf')

    @staticmethod
    def get_assembly_name(db):
        assembly = db.replace("_target_decoy_database.fasta", "")
        return assembly

    def filter_table(self, file, output, pattern, utps_only=False):
        df = pd.read_csv(file, sep='\t')
        df = df[df["q-value"] < self.args.qvalue]
        if not self.args.keepAnnotated:
            df = df[df["proteinIds"].str.contains(pattern, regex=False) == False]
        df = df[df["proteinIds"].str.contains("contaminant_", regex=False) == False]
        df = df[df["proteinIds"].str.contains("rev_") == False]
        if utps_only:
            df = df[df["proteinIds"].str.contains(",", regex=False) == False]
        df = self.__filter_min_replicates(df)
        df.to_csv(output, sep='\t', index=False)

    # def flatten_orf_duplicates_pin(self, pin, output):
    #     flattened_entries = []
    #     with open(pin, 'r') as handler, open(output, 'w') as outfile:
    #         lines = handler.readlines()
    #         new_lines = []
    #         for line in lines:
    #             cols = line.rstrip().split('\t')[20:]
    #             coords_checker = []
    #             proteins = []
    #             decoy_checker = []
    #             for col in cols:
    #                 if col not in flattened_entries:
    #                     flattened_entries.append(f'{col}\n')
    #                 if col != 'Proteins':
    #                     if '_F:' in col:
    #                         splat = col.split(":")
    #                         first = splat[0]
    #                         if '+' in first:
    #                             chrom = f'+{first.split("+")[1]}'
    #                         else:
    #                             chrom = f'-{first.split("-")[1]}'
    #                         coords = splat[1].split("_")[0]
    #                         locus = f'{chrom}:{coords}'
    #                         if 'rev_' not in col:
    #
    #                             if locus not in coords_checker:
    #                                 coords_checker.append(locus)
    #                                 proteins.append(col)
    #                         else:
    #                             if locus not in decoy_checker:
    #                                 proteins.append(col)
    #                                 decoy_checker.append(locus)
    #                     else:
    #                         proteins.append(col)
    #                 else:
    #                     proteins.append(col)
    #
    #
    #             prot_list = '\t'.join(proteins)
    #             other_cols = '\t'.join(line.rstrip().split('\t')[:20])
    #             new_line = f'{other_cols}\t{prot_list}\n'
    #             new_lines.append(new_line)
    #             # if any(p.startswith('rev_') for p in cols) and any(p.startswith("MSTRG") for p in cols):
    #             #     print('\n')
    #             #     print(cols)
    #             #     print(coords_checker)
    #             #     print(prot_list)
    #             #     print(other_cols)
    #             #     print(new_line)
    #         outfile.writelines(new_lines)
    #     with open(f'{self.outdir}/flattened_entries.txt', 'w') as out_flat:
    #         out_flat.writelines(flattened_entries)
    # def process_line(self, line, flattened_entries):
    #     cols = line.rstrip().split('\t')[20:]
    #     coords_checker = []
    #     proteins = []
    #     decoy_checker = []
    #     for col in cols:
    #         # print(col)
    #         if col not in flattened_entries:
    #             flattened_entries.append(f'{col}\n')
    #         if col != 'Proteins':
    #             if '_F:' in col:
    #                 splat = col.split(":")
    #                 first = splat[0]
    #                 if '+' in first:
    #                     chrom = f'+{first.split("+")[1]}'
    #                 else:
    #                     chrom = f'-{first.split("-")[1]}'
    #                 coords = splat[1].split("_")[0]
    #                 locus = f'{chrom}:{coords}'
    #                 if 'rev_' not in col:
    #                     if locus not in coords_checker:
    #                         coords_checker.append(locus)
    #                         proteins.append(col)
    #                 else:
    #                     if locus not in decoy_checker:
    #                         proteins.append(col)
    #                         decoy_checker.append(locus)
    #             else:
    #                 proteins.append(col)
    #         else:
    #             proteins.append(col)
    #     prot_list = '\t'.join(proteins)
    #     other_cols = '\t'.join(line.rstrip().split('\t')[:20])
    #     new_line = f'{other_cols}\t{prot_list}\n'
    #     return new_line
    #
    # def flatten_orf_duplicates_pin(self, pin, output):
    #
    #     flattened_entries = []
    #     with open(pin, 'r') as handler, open(output, 'w') as outfile:
    #         lines = handler.readlines()
    #         with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.threads) as executor:
    #             # Process each line in parallel with tqdm progress bar
    #             with tqdm(total=len(lines)) as pbar:
    #                 new_lines = []
    #                 for line in executor.map(lambda line: self.process_line(line, flattened_entries), lines):
    #                     new_lines.append(line)
    #                     pbar.update(1)
    #         outfile.writelines(new_lines)
    #
    #     with open(f'{self.outdir}/flattened_entries.txt', 'w') as out_flat:
    #         out_flat.writelines(flattened_entries)

    @staticmethod
    def flatten_orf_duplicates(proteins):
        """
        for a given percolator protein list, separated by comma, choose a single smORF in the case of multiple
        with the same genome coordinates (but different transcript isoforms, possibly) sharing the same peptide.
        """
        flattened = []
        for prot_list in proteins:
            if ',' not in prot_list:
                flattened.append(prot_list)
            else:
                coords_checker = []
                proteins = prot_list.split(",")
                for prot in proteins:
                    splat = prot.split(":")
                    first = splat[0]
                    if '+' in first:
                        chrom = f'+{first.split("+")[1]}'
                    else:
                        chrom = f'-{first.split("-")[1]}'
                    coords = splat[1].split("_")[0]
                    locus = f'{chrom}:{coords}'
                    if locus not in coords_checker:
                        coords_checker.append(locus)
                        flattened.append(prot)
        return flattened

    def __filter_min_replicates(self, df):
        """ Filter the peptides based on a minimum number of replicates (r) they could be identified in. """
        r = self.args.minReplicates
        replicates = {}
        modified = {}
        files, peptides = df["PSMId"].tolist(), df["peptide"].tolist()
        for file, peptide in zip(files, peptides):

            # I want to consider all modified peptides as the same peptide. This way, if we have only two
            # identifications for a peptide, one with +57 and the other without it, it would be considered as
            # being found in two replicates
            mod_peptide = peptide
            peptide = re.sub(r'[^A-Z]', '', peptide)
            if peptide not in modified:
                modified[peptide] = []
            modified[peptide].append(mod_peptide)

            rep = file.split(".")[0]
            if peptide not in replicates:
                replicates[peptide] = []
            if rep not in replicates[peptide]:
                replicates[peptide].append(rep)
        # print("replicates", replicates)
        # print("modified", modified)
        selected = []
        for peptide in replicates:
            if len(replicates[peptide]) >= int(r):
                for mod in modified[peptide]:
                    selected.append(mod)
        # print(selected)
        df = df[df["peptide"].isin(selected)]
        # print(df)
        return df

    def create_merged_gtf(self):
        utps = f'{self.summarizedResultsDir}/merged/microproteins_utps_150.fasta_blast_filt.fasta'
        mps = f'{self.summarizedResultsDir}/merged/microproteins_150.fasta_blast_filt.fasta'

        def get_entries(fasta):
            entries = []
            seqs = []
            records = SeqIO.parse(fasta, 'fasta')
            for record in records:
                seq = str(record.seq)
                if seq not in seqs:
                    entries.append(str(record.description))
                    seqs.append(seq)
            return entries

        utps_mps = get_entries(utps)
        all_mps = get_entries(mps)

        def filter_gtf(entries, gtf_list, output):
            out_gtf = []
            for gtf in gtf_list:
                with open(gtf, 'r') as handler:
                    lines = handler.readlines()
                    for line in lines:
                        for entry in entries:
                            if entry in line:
                                out_gtf.append(line)
            with open(output, 'w') as outfile:
                outfile.writelines(out_gtf)
            return output

        genn = group_folder_generator(self.resultsDir)
        gtfs = []
        for content in genn:
            gtf = f'{content.dbDir}/microproteins.gtf'
            gtfs.append(gtf)
        utps_gtf = filter_gtf(entries=utps_mps, gtf_list=gtfs,
                              output=f'{self.summarizedResultsDir}/merged/microproteins_utps.gtf')
        all_gtf = filter_gtf(entries=all_mps, gtf_list=gtfs,
                             output=f'{self.summarizedResultsDir}/merged/microproteins.gtf')

    def __protein_fdr(self, rescore=False):
        def count_substring_occurrences(string, substring):
            count = 0
            start_index = 0

            # Loop through the string and find occurrences of the substring
            while True:
                # Find the index of the next occurrence of the substring
                index = string.find(substring, start_index)

                # If no further occurrences are found, break out of the loop
                if index == -1:
                    break

                # Increment the count of occurrences
                count += 1

                # Update the start index for the next iteration
                start_index = index + 1

            return count

        if rescore:
            df = pd.read_csv(f'{self.rescorePostProcessDir}/group/db/proteins.txt', sep='\t')
        else:
            df = pd.read_csv(f'{self.postProcessDir}/group/db/proteins.txt', sep='\t')
        df = df[df["q-value"] != "q-value"]
        df["q-value"] = df["q-value"].astype(float)
        df = df[df["q-value"] < 0.01]
        # if not self.args.keepAnnotated:
        df = df[df["ProteinId"].str.contains("ANNO") == False]
        df = df[df["ProteinId"].str.contains("MOUSE") == False]
        df = df[df["ProteinId"].str.contains("contaminant") == False]
        df = df[df["ProteinId"].str.contains("rev_") == False]
        filtered_proteins = []
        proteins = df["ProteinId"].tolist()
        for prot in proteins:
            prot_list = prot.split(",")
            for protein in prot_list:
                # if 'ANNO' not in protein:
                if '_ANNO' not in prot:
                    filtered_proteins.append(protein)
                else:
                    if self.args.keepAnnotated:
                        matches = count_substring_occurrences(prot, 'ANNO')
                        if 'ANNO' in protein and matches <= 1:
                            filtered_proteins.append(protein)

        return filtered_proteins

    def create_fasta(self, peptides, psms, database, results_dir, smorfs=False, utps=False):
        proteins = {}
        records = SeqIO.parse(database, 'fasta')
        for record in records:
            if smorfs:
                if len(str(record.seq)) <= self.args.maxORFLength:
                    proteins[str(record.id)] = str(record.seq)
            else:
                proteins[str(record.id)] = str(record.seq)
        df = pd.read_csv(peptides, sep='\t')
        to_write = []
        unch = []
        microproteins = df["proteinIds"].tolist()
        if self.args.proteinFDR:
            filtered_proteins = self.__protein_fdr()

        for mp in microproteins:
            add = False
            prot_list = mp.split(",")
            for prot in prot_list:
                if utps:
                    if ',' not in prot:
                        # prot_list = mp.split(",")
                        # for prot in prot_list:
                        # if smorfs:
                        if prot in proteins:
                            add = True
                            # to_write.append(f'>{prot}\n{proteins[prot]}\n')
                else:
                    prot_list = prot.split(",")
                    for protein in prot_list:
                        if protein in proteins:
                            add = True
                            # to_write.append(f'>{protein}\n{proteins[protein]}\n')

                    # for prot_list in proteins:
                    #     proteins_splat = prot_list.split(",")
                    #     for protein in proteins_splat:
                if self.args.proteinFDR:
                    if prot in filtered_proteins:
                        add = True
                    else:
                        add = False
                # else:
                #     add = True

                if add:
                    if prot in proteins:
                        if prot.endswith('_UNCH'):
                            unch.append(f'>{prot}\n{proteins[prot]}\n')
                        else:
                            to_write.append(f'>{prot}\n{proteins[prot]}\n')

        if smorfs and utps:
            fasta = f'{results_dir}/microproteins_utps_150.fasta'
            unch_fasta = f'{results_dir}/uncharacterized_microproteins_utps_150.fasta'
        elif smorfs and not utps:
            fasta = f'{results_dir}/microproteins_150.fasta'
            unch_fasta = f'{results_dir}/uncharacterized_microproteins_150.fasta'
        elif utps and not smorfs:
            fasta = f'{results_dir}/proteins.fasta'

        with open(fasta, 'w') as out:
            out.writelines(to_write)
        if smorfs:
            if self.args.includeLowAnnotation:
                with open(unch_fasta, 'w') as outfile:
                    outfile.writelines(unch)
        return fasta

    def create_gtf(self, microproteins_fasta, results_dir, translation_gtf):
        records = SeqIO.parse(microproteins_fasta, 'fasta')
        entries = []
        for record in records:
            entries.append(str(record.description))

        filtered_gtf = []
        with open(translation_gtf, 'r') as gtf:
            lines = gtf.readlines()
            for line in lines:
                if any(mp in line for mp in entries):
                    filtered_gtf.append(line)
        with open(f'{results_dir}/microproteins.gtf', 'w') as out:
            out.writelines(filtered_gtf)

    def filter_microproteins(self):
        groups = os.listdir(self.resultsDir)
        cat_proteins = []
        cat_microproteins = []
        cat_uniq_proteins = []
        cat_uniq_proteins_seqs = []
        cat_uniq_mps_seqs = []
        cat_uniq_microproteins = []
        for group in groups:
            group_dir = f'{self.resultsDir}/{group}'
            if os.path.isdir(group_dir):
                dbs = os.listdir(group_dir)
                for db in dbs:

                    to_write = []
                    db_dir = f'{group_dir}/{db}'
                    if os.path.exists(f'{db_dir}/microproteins.fasta'):
                        records = SeqIO.parse(f'{db_dir}/microproteins.fasta', 'fasta')
                        for record in records:
                            prot = f'>{str(record.description)}\n{str(record.seq)}\n'
                            cat_proteins.append(prot)
                            seq = str(record.seq)
                            if seq not in cat_uniq_proteins_seqs:
                                cat_uniq_proteins_seqs.append(seq)
                                cat_uniq_proteins.append(prot)
                            if len(str(record.seq)) <= 150:
                                if seq not in cat_uniq_mps_seqs:
                                    cat_uniq_mps_seqs.append(seq)
                                    cat_uniq_microproteins.append(prot)
                                to_write.append(f'>{str(record.description)}\n{str(record.seq)}\n')
                                cat_microproteins.append(prot)
                        with open(f'{db_dir}/microproteins_150.fasta', 'w') as out:
                            out.writelines(to_write)
        with open(f'{self.resultsDir}/cat_microproteins.fasta', 'w') as out_mp, \
                open(f'{self.resultsDir}/cat_proteins.fasta', 'w') as out_prots:
            out_mp.writelines(cat_microproteins)
            out_prots.writelines(cat_proteins)
        with open(f'{self.resultsDir}/cat_unique_microproteins.fasta', 'w') as out_mp, open(
                f'{self.resultsDir}/cat_unique_proteins.fasta', 'w') as out_prot:
            out_mp.writelines(cat_uniq_microproteins)
            out_prot.writelines(cat_uniq_proteins)

    def merge_gtf_results(self):
        groups = os.listdir(self.resultsDir)
        group_gtf = {}
        for group in groups:
            group_dir = f'{self.resultsDir}/{group}'
            if os.path.isdir(group_dir):
                dbs = os.listdir(group_dir)
                for db in dbs:
                    db_dir = f'{group_dir}/{db}'
                    gtf = f'{db_dir}/microproteins.gtf'
                    if os.path.exists(gtf):
                        assembly = self.get_assembly_name(db)
                        if assembly not in group_gtf:
                            group_gtf[assembly] = []
                        with open(gtf, 'r') as handler:
                            lines = handler.readlines()
                            for line in lines:
                                if line not in group_gtf[assembly]:
                                    group_gtf[assembly].append(line)

        for assembly in group_gtf:
            # print(assembly)
            with open(f'{self.resultsDir}/{assembly}.gtf', 'w') as outfile:
                outfile.writelines(group_gtf[assembly])

    # def create_cds_gtf(self):

    def merge_fasta(self):
        genn = group_folder_generator(self.resultsDir)
        # print(genn)
        merged = {
            'merged': {'microproteins_150.fasta': [], 'microproteins_utps_150.fasta': [], 'proteins_utp.fasta': []}}
        for content in genn:
            # if content.
            # print(content.file)
            # print(content.dbDir)

            def add_file(subset, file, folder):
                if folder not in merged:
                    merged[folder] = {}
                if subset not in merged[folder]:
                    merged[folder][subset] = []

                check_multiple_dirs([f'{self.summarizedResultsDir}/{folder}', f'{self.summarizedResultsDir}/merged'])
                records = SeqIO.parse(file, 'fasta')
                for record in records:
                    entry = f'>{str(record.description)}\n{str(record.seq)}\n'
                    if entry not in merged[folder][subset]:
                        merged[folder][subset].append(entry)
                    if entry not in merged['merged'][subset]:
                        merged['merged'][subset].append(entry)

            add_file(subset='microproteins_150.fasta', file=f'{content.dbDir}/microproteins_150.fasta',
                     folder=content.db)
            add_file(subset='microproteins_utps_150.fasta', file=f'{content.dbDir}/microproteins_utps_150.fasta',
                     folder=content.db)
            add_file(subset='proteins_utp.fasta', file=f'{content.dbDir}/proteins.fasta', folder=content.db)

        for folder in merged:
            for file in merged[folder]:
                with open(f'{self.summarizedResultsDir}/{folder}/{file}', 'w') as outfile:
                    outfile.writelines(merged[folder][file])

    # def merge_fasta_subset(self, subset):

    # cat_fasta = {}
    # checker = {}
    # results = os.listdir(self.resultsDir)
    # for subdir in results:
    #     if os.path.isdir(f'{self.resultsDir}/{subdir}'):
    #         assemblies = os.listdir(f'{self.resultsDir}/{subdir}')
    #         for assembly in assemblies:
    #             if assembly not in cat_fasta:
    #                 cat_fasta[assembly] = []
    #                 checker[assembly] = []
    #             assembly_dir = f'{self.resultsDir}/{subdir}/{assembly}'
    #             records = SeqIO.parse(f'{assembly_dir}/proteins.fasta', 'fasta')
    #             for record in records:
    #                 seq = str(record.seq)
    #                 if seq not in checker[assembly]:
    #                     cat_fasta[assembly].append(f'>{str(record.description)}\n{seq}\n')
    #                     checker[assembly].append(seq)
    # for assembly in cat_fasta:
    #     with open(f'{self.summarizedResultsDir}/{assembly}_proteins.fasta', 'w') as out:
    #         out.writelines(cat_fasta[assembly])

    # def merge_utp_fasta(self):
    #     generator = group_folder_generator(self.resultsDir)

    def create_cds_gtf(self):
        fasta = os.listdir(f'{self.summarizedResultsDir}')
        for file in fasta:
            entries = []

            if file.endswith("microproteins.fasta"):
                records = SeqIO.parse(f'{self.summarizedResultsDir}/{file}', 'fasta')
                for record in records:
                    entries.append(str(record.description))
            assembly = '_'.join(file.split("_")[:-4])
            nuc_fasta = f'{self.translationDir}/{assembly}/{assembly}.nuc'
            cds_fasta = []
            records = SeqIO.parse(nuc_fasta, 'fasta')
            for record in records:
                if str(record.description) in entries:
                    cds_fasta.append(f'>{str(record.description)}\n{str(record.seq)}\n')

            with open(f'{self.summarizedResultsDir}/{assembly}_nuc.fasta', 'w') as outfile:
                outfile.writelines(cds_fasta)

    def get_utps(self):
        subdirs = os.listdir(self.postProcessDir)
        peptides = {}
        for folder in subdirs:
            subdir = f'{self.postProcessDir}/{folder}'
            asses = os.listdir(subdir)
            for ass in asses:
                if ass not in peptides:
                    peptides[ass] = []
                peptides_df = f'{subdir}/{ass}/peptides_filtered.txt'
                df = pd.read_csv(peptides_df, sep='\t')
                df = df[df["proteinIds"].str.contains(",") == False]
                proteins = df["proteinIds"].tolist()
                for prot in proteins:
                    line = f'{prot}\n'
                    if line not in peptides[ass]:
                        peptides[ass].append(line)
        for ass in peptides:
            with open(f'{self.summarizedResultsDir}/{ass}_UTPs.txt', 'w') as outfile:
                outfile.writelines(peptides[ass])
