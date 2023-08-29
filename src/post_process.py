import os
import sys

import pandas as pd
from Bio import SeqIO

from .utils import group_folder_generator, check_multiple_dirs
from .pipeline_config import PipelineStructure


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
                        files_to_convert[db_prefix] = []   # this is to separate different DBs for a same group
                        # for instance, multiple transcriptome assemblies but the same mass spec data

                    pin_dir = f'{self.peptideSearchDir}/{group}/pin_files'
                    if not os.path.exists(pin_dir):
                        os.mkdir(pin_dir)

                    pin_dir_db = f'{pin_dir}/{db_prefix}'   # the pin folder should contain a different folder
                    if not os.path.exists(pin_dir_db):    # for each transcriptome assembly
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
                    db_path = f'{group_dir}/{db_dir}'   # the directories are named after the database that was used
                    pin_files = os.listdir(db_path)     # for the search
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

    def percolate(self):
        groups = os.listdir(self.peptideSearchDir)
        for group in groups:
            group_dir = f'{self.peptideSearchDir}/{group}'
            dbs = os.listdir(group_dir)
            for db_dir in dbs:
                if 'target' in db_dir and db_dir.endswith(".fasta"):
                    db_path = f'{group_dir}/{db_dir}'   # the directories are named after the database that was used
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
                        """
                        if file.endswith("pin"):
                            # pin += f' {db_path}/{file}'
                            # percolator_input = f'{db_path}/percolator_input'
                            # pin = f'{percolator_input}/{db_dir}.pin'

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

                    cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                                     f'--post-processing-tdc --results-psms {outdir}/{file}_psm.txt --results-peptides ' \
                                     f'{outdir}/{file}_peptides.txt --no-terminate --num-threads {self.threads} ' \
                                     f'-X {outdir}/pout.xml {db_path}/percolator_input/{db_dir}.pin'
                    # cmd_percolator = f'percolator --protein-decoy-pattern rev_ ' \
                    #                  f'--post-processing-mix-max --results-psms {outdir}/psm.txt --results-peptides ' \
                    #                  f'{outdir}/peptides.txt -f auto --test-each-iteration -I separate  ' \
                    #                  f'-X {outdir}/pout.xml {pin}'
                    os.system(cmd_percolator)
    def fix_multiple_columns(self):
        groups = os.listdir(self.postProcessDir)
        for group in groups:
            dbs = os.listdir(f'{self.postProcessDir}/{group}')
            for db in dbs:
                db_dir = f'{self.postProcessDir}/{group}/{db}'
                if self.args.cat:
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


                    fasta = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                              database=f'{self.outdir}/databases/{db}', results_dir=results_dir,
                                              smorfs=True)

                    fasta_utps = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                                   database=f'{self.outdir}/databases/{db}', results_dir=results_dir,
                                                   smorfs=True, utps=True)
                    fasta_proteins = self.create_fasta(peptides=peptide_filtered, psms=psm_filtered,
                                                   database=f'{self.outdir}/databases/{db}', results_dir=results_dir,
                                                   smorfs=False, utps=True)


                    assembly = db.replace("_target_decoy_database.fasta", "")
                    if self.args.create_gtf == True:
                        self.create_gtf(microproteins_fasta=fasta, results_dir=results_dir,
                                        translation_gtf=f'{self.translationDir}/{assembly}/{assembly}_ORFs.gtf')

    @staticmethod
    def get_assembly_name(db):
        assembly = db.replace("_target_decoy_database.fasta", "")
        return assembly

    @staticmethod
    def filter_table(file, output, pattern, utps_only=False):
        df = pd.read_csv(file, sep='\t')
        df = df[df["q-value"] < 0.01]
        df = df[df["proteinIds"].str.contains(pattern, regex=False) == False]
        df = df[df["proteinIds"].str.contains("contaminant_", regex=False) == False]
        df = df[df["proteinIds"].str.contains("rev_") == False]
        if utps_only:
            df = df[df["proteinIds"].str.contains(",", regex=False) == False]
        df.to_csv(output, sep='\t', index=False)

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
        utps_gtf = filter_gtf(entries=utps_mps, gtf_list=gtfs, output=f'{self.summarizedResultsDir}/merged/microproteins_utps.gtf')
        all_gtf = filter_gtf(entries=all_mps, gtf_list=gtfs, output=f'{self.summarizedResultsDir}/merged/microproteins.gtf')


    def create_fasta(self, peptides, psms, database, results_dir, smorfs=False, utps=False):
        proteins = {}
        records = SeqIO.parse(database, 'fasta')
        for record in records:
            if smorfs:
                if len(str(record.seq)) <= 150:
                    proteins[str(record.id)] = str(record.seq)
            else:
                proteins[str(record.id)] = str(record.seq)
        df = pd.read_csv(peptides, sep='\t')
        to_write = []
        microproteins = df["proteinIds"].tolist()

        for prot in microproteins:
                if utps:
                    if ',' not in prot:
                        # prot_list = mp.split(",")
                        # for prot in prot_list:
                        # if smorfs:
                        if prot in proteins:
                            to_write.append(f'>{prot}\n{proteins[prot]}\n')
                else:
                    prot_list = prot.split(",")
                    for protein in prot_list:
                        if protein in proteins:
                            to_write.append(f'>{protein}\n{proteins[protein]}\n')

        if smorfs and utps:
            fasta = f'{results_dir}/microproteins_utps_150.fasta'
        elif smorfs and not utps:
            fasta = f'{results_dir}/microproteins_150.fasta'
        elif utps and not smorfs:
            fasta = f'{results_dir}/proteins.fasta'
        with open(fasta, 'w') as out:
            out.writelines(to_write)
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
        with open(f'{self.resultsDir}/cat_unique_microproteins.fasta', 'w') as out_mp, open(f'{self.resultsDir}/cat_unique_proteins.fasta', 'w') as out_prot:
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
            print(assembly)
            with open(f'{self.resultsDir}/{assembly}.gtf', 'w') as outfile:
                outfile.writelines(group_gtf[assembly])

    # def create_cds_gtf(self):

    def merge_fasta(self):
        genn = group_folder_generator(self.resultsDir)
        print(genn)
        merged = {'merged': {'microproteins_150.fasta': [], 'microproteins_utps_150.fasta': [], 'proteins_utp.fasta': []}}
        for content in genn:
            # if content.
            print(content.file)
            print(content.dbDir)
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

            add_file(subset='microproteins_150.fasta', file=f'{content.dbDir}/microproteins_150.fasta', folder=content.db)
            add_file(subset='microproteins_utps_150.fasta', file=f'{content.dbDir}/microproteins_utps_150.fasta', folder=content.db)
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





