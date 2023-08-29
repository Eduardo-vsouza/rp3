import os
import sys

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure
from ..decoy import Decoy
from ..post_process import PercolatorPostProcessing
from ..utils import GTFFilter, GTFGatherer


class PeptideReScoring(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.MSFraggerPath = self.toolPaths["MSFragger"]
        self.params = []


        self.check_dirs([self.rescoreDir, self.rescoreDatabaseDir, self.rescoreSearchDir, self.rescorePostProcessDir,
                         self.rescoreSummarizedResultsDir])

        self.rescoreDatabase = f'{self.rescoreDatabaseDir}/rescore_target_decoy_database.fasta'

    def generate_databases(self):
        print("Generating databases for rescoring")
        checker = []
        fasta = []
        records = SeqIO.parse(f'{self.microproteinsBlast}', 'fasta')  # microproteins identified with the first search
        for record in records:
            seq = str(record.seq)
            if seq not in checker:
                checker.append(seq)
                fasta.append(f'>{str(record.description).replace(" ", "_")}\n{str(record.seq)}\n')
        reference = SeqIO.parse(self.args.proteome, 'fasta')
        for record in reference:
            fasta.append(f'>{str(record.description).replace(" ", "_")}_ANNO\n{str(record.seq)}\n')
        with open(f'{self.rescoreDatabaseDir}/rescore_target_database.fasta', 'w') as outfile:
            outfile.writelines(fasta)

        decoy = Decoy(db=f'{self.rescoreDatabaseDir}/rescore_target_database.fasta')
        decoy.reverse_sequences()
        decoy.add_contaminants()
        decoy.to_fasta(output=f'{self.rescoreDatabase}')

    def re_search_peptides(self, min_pep_len=7, max_pep_len=50):
        print("Performing peptide search")
        pattern = self.args.msPattern
        groups = os.listdir(self.args.mzml)
        for group in groups:
            group_files = ''
            group_dir = f'{self.args.mzml}/{group}'
            files = os.listdir(group_dir)
            for file in files:
                if file.endswith(pattern):
                    print(file)
                    group_files += f' {group_dir}/{file}'
            cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                  f'--database_name {self.rescoreDatabase} --decoy_prefix rev ' \
                  f'--num_threads {self.args.threads} --digest_min_length {min_pep_len} ' \
                  f'--digest_max_length {max_pep_len}{group_files}'
            os.system(cmd)
            out_group_dir = f'{self.rescoreSearchDir}/{group}'
            self.check_dirs([out_group_dir])
            self.params.append(cmd)
            # print("\n TARGET \n\n")
            for file in group_files.split(" "):
                if file.endswith(pattern):
                    mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                         f'{out_group_dir}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                    os.system(mv)

    def re_percolate(self):
        self.__merge_pin_files()
        groups = os.listdir(self.rescoreSearchDir)
        for group in groups:
            group_dir = f'{self.rescoreSearchDir}/{group}'
            percolator_input = f'{group_dir}/percolator_input'
            pin = f'{percolator_input}/{group}.pin'

            group_outdir = f'{self.rescorePostProcessDir}/{group}'
            if not os.path.exists(group_outdir):
                os.mkdir(group_outdir)

            cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                             f'--post-processing-tdc --results-psms {group_outdir}/psm.txt --results-peptides ' \
                             f'{group_outdir}/peptides.txt --no-terminate --num-threads {self.args.threads} ' \
                             f'-X {group_outdir}/pout.xml {pin}'
            os.system(cmd_percolator)

    def __merge_pin_files(self):
        groups = os.listdir(self.rescoreSearchDir)
        for group in groups:
            group_dir = f'{self.rescoreSearchDir}/{group}'
            pin_files = os.listdir(group_dir)
            merged_pin = []
            for pin in pin_files:
                if pin.endswith(".pin"):
                    with open(f'{group_dir}/{pin}', 'r') as handler:
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
            percolator_input = f'{group_dir}/percolator_input'
            if not os.path.exists(percolator_input):
                os.mkdir(percolator_input)
            with open(f'{percolator_input}/{group}.pin', 'w') as outfile:
                outfile.writelines(merged_pin)


    def re_assess_fdr(self):
        groups = os.listdir(self.rescorePostProcessDir)
        microproteins = {}
        records = SeqIO.parse(self.rescoreDatabase, 'fasta')
        for record in records:
            microproteins[str(record.description).replace(" ", "_")] = str(record.seq)

        perc = PercolatorPostProcessing(args=self.args)

        for group in groups:
            fasta = []
            fixed = f'{self.rescorePostProcessDir}/{group}/peptides_fixed.txt'
            perc.fix(file=f'{self.rescorePostProcessDir}/{group}/peptides.txt',
                     output=fixed)
            df = pd.read_csv(fixed, sep='\t')
            df = df[df["q-value"] != "q-value"]
            df = df[df["q-value"] < 0.01]
            df = df[df["proteinIds"].str.contains("ANNO") == False]
            df = df[df["proteinIds"].str.contains("MOUSE") == False]
            df = df[df["proteinIds"].str.contains("contaminant") == False]
            df = df[df["proteinIds"].str.contains("rev_") == False]
            proteins = df["proteinIds"].tolist()
            for prot_list in proteins:
                proteins_splat = prot_list.split(",")
                for protein in proteins_splat:
                    fasta.append(f'>{protein}\n{microproteins[protein]}\n')
            with open(f'{self.rescorePostProcessDir}/{group}/filtered_rescored_smorfs.fasta', 'w') as handler:
                handler.writelines(fasta)

    def merge_results(self):
        groups = os.listdir(self.rescorePostProcessDir)
        fasta = []
        for group in groups:
            records = SeqIO.parse(f'{self.rescorePostProcessDir}/{group}/filtered_rescored_smorfs.fasta', 'fasta')
            for record in records:
                if len(str(record.seq)) <= 150:
                    fasta.append(f'>{str(record.description)}\n{str(record.seq)}\n')
        with open(self.rescoredMicroproteinsFasta, 'w') as outfile:
            outfile.writelines(fasta)

    def filter_gtf(self):
        gatherer = GTFGatherer(args=self.args)
        gatherer.cat_gtfs()
        data = GTFFilter(args=self.args, gtf=f'{self.mergedFullGTF}',
                         fasta=self.rescoredMicroproteinsFasta)
        data.get_entries()
        # data.filter_gtf(output=self.rescoredMicroproteinsGTF)
        data.filter_gtf_tmp()