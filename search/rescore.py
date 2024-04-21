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
        fasta_file = self.microproteinsBlast
        if not os.path.exists(self.microproteinsBlast):
            fasta_file = self.uniqueMicroproteins
        if self.args.smorfUTPs:
            fasta_file = self.utpsMicroproteinsBlast
        records = SeqIO.parse(f'{fasta_file}', 'fasta')  # microproteins identified with the first search
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
                    group_files += f' {group_dir}/{file}'
            if self.args.hlaPeptidomics:
                cmd = f'java -Xmx256g -jar {self.MSFraggerPath} --output_format pin ' \
                      f'--database_name {self.rescoreDatabase} --decoy_prefix rev --search_enzyme_name nonspecific ' \
                      f'--num_threads {self.args.threads} --fragment_mass_tolerance 20 --num_enzyme_termini 0 ' \
                      f'--precursor_true_tolerance 6 --digest_mass_range 500.0_1500.0 ' \
                      f'--max_fragment_charge 3 --search_enzyme_cutafter ARNDCQEGHILKMFPSTWYV ' \
                      f'--digest_min_length 8 --digest_max_length 25{group_files}'
            else:
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
                             f'-X {group_outdir}/pout.xml --picked-protein {self.rescoreDatabase} --results-proteins {group_outdir}/proteins.txt {pin}'
            os.system(cmd_percolator)

    def re_percolate_all_pins(self):
        self.__merge_all_pin_files()
        pin = f'{self.rescoreDir}/all_pins.pin'
        if self.args.msBooster:
            pin = self.mergedBoosterPin
        group_outdir = f'{self.rescorePostProcessDir}/group'
        self.check_dirs([group_outdir])
        cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                         f'--post-processing-tdc --results-psms {group_outdir}/psm.txt --results-peptides ' \
                         f'{group_outdir}/peptides.txt --no-terminate --num-threads {self.args.threads} ' \
                         f'-X {group_outdir}/pout.xml --picked-protein {self.rescoreDatabase} --results-proteins {group_outdir}/proteins.txt {pin}'
        os.system(cmd_percolator)

    def __merge_all_pin_files(self):
        groups = os.listdir(self.rescoreSearchDir)
        merged_pin = []
        for group in groups:
            group_dir = f'{self.rescoreSearchDir}/{group}'
            pin_files = os.listdir(group_dir)
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
        with open(f'{self.rescoreDir}/all_pins.pin', 'w') as out:
            out.writelines(merged_pin)

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
        # print(microproteins)
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
            if self.args.smorfUTPs:
                df = df[df["proteinIds"].str.contains(",") == False]
            proteins = df["proteinIds"].tolist()

            if self.args.proteinFDR:
                filtered_proteins = self.__protein_fdr(group=group)

            for prot_list in proteins:
                # print("protlist", prot_list)
                proteins_splat = prot_list.split(",")
                for protein in proteins_splat:
                    if self.args.proteinFDR:
                        if protein in filtered_proteins:
                            add = True
                        else:
                            add = False
                    else:
                        add = True
                    if add:
                        fasta.append(f'>{protein}\n{microproteins[protein]}\n')
            with open(f'{self.rescorePostProcessDir}/{group}/filtered_rescored_smorfs.fasta', 'w') as handler:
                handler.writelines(fasta)

    def __protein_fdr(self, group):
        df = pd.read_csv(f'{self.rescorePostProcessDir}/{group}/proteins.txt', sep='\t')
        df = df[df["q-value"] != "q-value"]
        df = df[df["q-value"] < 0.01]
        df = df[df["ProteinId"].str.contains("ANNO") == False]
        df = df[df["ProteinId"].str.contains("MOUSE") == False]
        df = df[df["ProteinId"].str.contains("contaminant") == False]
        df = df[df["ProteinId"].str.contains("rev_") == False]
        filtered_proteins = []
        proteins = df["ProteinId"].tolist()
        for prot in proteins:
            prot_list = prot.split(",")
            for protein in prot_list:
                filtered_proteins.append(protein)
        return filtered_proteins

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

    def filter_gtf(self, rescored=True):
        gatherer = GTFGatherer(args=self.args)
        gatherer.cat_gtfs()
        fasta = self.select_fasta()
        records = SeqIO.parse(fasta, 'fasta')
        checker = []
        nr_fasta = []
        for record in records:
            seq = str(record.seq)
            if seq not in checker:
                checker.append(seq)
                nr_fasta.append(f'>{str(record.description)}\n{seq}\n')
        # if rescored:
        #     fasta = self.rescoredMicroproteinsFasta
        # else:
        #     fasta = self.uniqueMicroproteinsNRFasta
        with open(self.uniqueMicroproteinsNRFasta, 'w') as outfile:
            outfile.writelines(nr_fasta)
        data = GTFFilter(args=self.args, gtf=f'{self.mergedFullGTF}',
                         fasta=fasta)
        data.get_entries()
        # data.filter_gtf(output=self.rescoredMicroproteinsGTF)
        data.filter_gtf_tmp()
