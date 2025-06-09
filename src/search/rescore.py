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
        self.args.rescored = True

        self.check_dirs([self.rescoreDir, self.rescoreDatabaseDir, self.rescoreSearchDir, self.rescorePostProcessDir,
                         self.rescoreSummarizedResultsDir])

        self.rescoreDatabase = f'{self.rescoreDatabaseDir}/rescore_target_decoy_database.fasta'

    def generate_databases(self):
        print("Generating databases for rescoring")
        checker = []
        fasta = []
        # fasta_file = self.select_fasta()
        fasta_file = self.microproteinsBlast
        if not os.path.exists(self.microproteinsBlast):
            fasta_file = self.uniqueMicroproteins
        # if self.args.smorfUTPs:
        #     fasta_file = self.utpsMicroproteinsBlast
        records = SeqIO.parse(f'{fasta_file}', 'fasta')  # microproteins identified with the first search
        for record in records:
            seq = str(record.seq)
            if seq not in checker:
                checker.append(seq)
                fasta.append(f'>{str(record.description).replace(" ", "_")}\n{str(record.seq)}\n')
        reference = SeqIO.parse(self.args.proteome, 'fasta')
        for record in reference:
            if self.args.keepAnnotated:
                fasta.append(f'>{str(record.description).replace(" ", "_").replace(",", "_")}\n{str(record.seq)}\n')
            else:
                fasta.append(f'>{str(record.description).replace(" ", "_").replace(",", "_")}_ANNO\n{str(record.seq)}\n')
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

             # = os.listdir(group_dir)
            if not os.path.isdir(groups[0]):
                files = os.listdir(self.args.mzml)
                group = 'group'
                single_group = True
                group_dir = f'{self.args.mzml}'

            else:
                single_group = False
                files = os.listdir(f'{self.args.mzML}/{group}')
                group_dir = f'{self.args.mzml}/{group}'

            for file in files:
                if file.endswith(pattern):
                    group_files += f' {group_dir}/{file}'
            phospho = ''
            mod = ''
            i = 3
            if self.args.amidation:
                amida = f' --variable_mod_0{i} -0.9840_c*_1'
                i += 1
            else:
                amida = ''

            if self.args.pyroGlu:
                pyroglu = f' --variable_mod_0{i} -17.0265_nQ_1'
                i += 1
            else:
                pyroglu = ''
            if self.args.mod is not None:
                mod = f' --variable_mod_0{i} {self.args.mod}'
                i += 1
            if self.args.phosphorylation:
                phospho = f' --variable_mod_0{i} 79.9663_STY_3'
            if self.args.hlaPeptidomics:
                cmd = f'java -Xmx256g -jar {self.MSFraggerPath} --output_format pin ' \
                      f'--database_name {self.rescoreDatabase} --decoy_prefix rev --search_enzyme_name nonspecific ' \
                      f'--num_threads {self.args.threads}{phospho}{mod}{amida}{pyroglu} --fragment_mass_tolerance 20 --num_enzyme_termini 0 ' \
                      f'--precursor_true_tolerance 6 --digest_mass_range 500.0_1500.0 ' \
                      f'--max_fragment_charge 3 --search_enzyme_cutafter ARNDCQEGHILKMFPSTWYV ' \
                      f'--use_all_mods_in_first_search 1 --digest_min_length 8 --digest_max_length 25{group_files}'
            else:
                if self.args.quantifyOnly or self.args.quantify:
                    cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format tsv ' \
                          f'--database_name {self.rescoreDatabase} --decoy_prefix rev ' \
                          f'--num_threads {self.args.threads}{phospho}{mod}{amida}{pyroglu} --digest_min_length {min_pep_len} ' \
                          f'--use_all_mods_in_first_search 1 --digest_max_length {max_pep_len}{group_files}'
                    os.system(cmd)

                if not self.args.quantifyOnly:
                    cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format pin ' \
                          f'--database_name {self.rescoreDatabase} --decoy_prefix rev ' \
                          f'--num_threads {self.args.threads}{phospho}{mod}{amida} --digest_min_length {min_pep_len} ' \
                          f'--use_all_mods_in_first_search 1 --digest_max_length {max_pep_len}{group_files}'
                    os.system(cmd)

            out_group_dir = f'{self.rescoreSearchDir}/{group}'
            self.check_dirs([out_group_dir])
            self.params.append(cmd)
            # print("\n TARGET \n\n")
            for file in group_files.split(" "):
                if file.endswith(pattern):
                    if self.args.quantifyOnly:
                        suffix = 'tsv'
                        mv = f'mv {file.replace(f".{pattern}", f".{suffix}")} ' \
                             f'{out_group_dir}/{file.split("/")[-1].replace(f".{pattern}", f"_target.{suffix}")}'
                        os.system(mv)

                    else:
                        suffix = 'pin'
                        mv = f'mv {file.replace(f".{pattern}", f".{suffix}")} ' \
                             f'{out_group_dir}/{file.split("/")[-1].replace(f".{pattern}", f"_target.{suffix}")}'
                        os.system(mv)
                        if self.args.quantify:
                            suffix = 'tsv'
                            mv = f'mv {file.replace(f".{pattern}", f".{suffix}")} ' \
                                 f'{out_group_dir}/{file.split("/")[-1].replace(f".{pattern}", f"_target.{suffix}")}'
                            os.system(mv)
            if single_group:
                break

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
        if self.args.groupedFDR:
            self.check_dirs([self.rescoreGroupFDRDir, self.rescoreGroupPostProcessDir,
                             self.rescoreMPGroupdir, self.rescoreAnnoGroupDir])
            self.__assess_group_fdr(pin)
        else:
            group_outdir = f'{self.rescorePostProcessDir}/group'
            self.check_dirs([group_outdir])
            cmd_percolator = f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ ' \
                             f'--post-processing-tdc --results-psms {group_outdir}/psm.txt --results-peptides ' \
                             f'{group_outdir}/peptides.txt --no-terminate --num-threads {self.args.threads} ' \
                             f'-X {group_outdir}/pout.xml --picked-protein {self.rescoreDatabase} --results-proteins {group_outdir}/proteins.txt {pin}'
            os.system(cmd_percolator)

    def percolate_single(self):
        group_outdir = f'{self.rescorePostProcessDir}/group'
        self.check_dirs([group_outdir])
        pins = os.listdir(f'{self.rescoreSearchDir}/group')
        for pin in pins:
            if pin.endswith("edited.pin"):
                cmd = (f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ '
                       f'--post-processing-tdc --results-psms {group_outdir}/{pin}_psm.txt --results-peptides '
                       f'{group_outdir}/{pin}_peptides.txt --no-terminate --num-threads {self.args.threads} '
                       f'-X {group_outdir}/{pin}_pout.xml --picked-protein {self.rescoreDatabase} '
                       f'--results-proteins {group_outdir}/{pin}_proteins.txt {self.rescoreSearchDir}/group/{pin}')
                os.system(cmd)
        cmd = f'cat {group_outdir}/*_peptides.txt > {group_outdir}/peptides.txt'
        os.system(cmd)
        cmd = f'cat {group_outdir}/*_psm.txt > {group_outdir}/psm.txt'
        os.system(cmd)
        cmd = f'cat {group_outdir}/*_proteins.txt > {group_outdir}/proteins.txt'
        os.system(cmd)


    def __assess_group_fdr(self, pin):
        """
        Divides the .pin files into two subsets: one containing canonical proteins, and the other containing predicted
        microproteins, then calculates the FDR for each and merges the results together
        """

        self.__generate_group_pin(pin)
        print(f"--Assessing FDR for microproteins separately")
        cmd_mp = (f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ '
                  f'--post-processing-tdc --results-psms {self.rescoreMPGroupdir}/psm.txt --results-peptides '
                  f'{self.rescoreMPGroupdir}/peptides.txt --no-terminate --num-threads {self.args.threads} '
                  f'-X {self.rescoreMPGroupdir}/pout.xml --picked-protein {self.rescoreDatabase} '
                  f'--results-proteins {self.rescoreMPGroupdir}/proteins.txt {self.mpPinFile}')
        os.system(cmd_mp)
        print(f"--Assessing FDR for canonical proteins separately")
        cmd_anno = (f'{self.toolPaths["percolator"]} --protein-report-duplicates --protein-decoy-pattern rev_ '
                    f'--post-processing-tdc --results-psms {self.rescoreAnnoGroupDir}/psm.txt --results-peptides '
                    f'{self.rescoreAnnoGroupDir}/peptides.txt --no-terminate --num-threads {self.args.threads} '
                    f'-X {self.rescoreAnnoGroupDir}/pout.xml --picked-protein {self.rescoreDatabase} '
                    f'--results-proteins {self.rescoreAnnoGroupDir}/proteins.txt {self.annoPinFile}')
        os.system(cmd_anno)

    def __generate_group_pin(self, pin):
        mp_pin = []
        anno_pin = []
        with open(pin, 'r') as handler:
            lines = handler.readlines()
            mp_pin.append(lines[0])
            anno_pin.append(lines[0])
            proteins_col = len(lines[0].rstrip().split('\t'))-1
            print(lines[0].split('\t')[proteins_col])

            for line in lines:
                line = line.rstrip()
                cols = line.split('\t')
                anno = False
                mp = False
                proteins = cols[proteins_col:]
                for protein in proteins:
                    if 'ANNO' in protein:
                        anno = True
                    if '_F:' in protein:
                        mp = True
                if mp and not anno:
                    mp_pin.append(f'{line}\n')
                if anno:
                    anno_pin.append(f'{line}\n')
        with open(self.mpPinFile, 'w') as out:
            out.writelines(mp_pin)
        with open(self.annoPinFile, 'w') as out:
            out.writelines(anno_pin)

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

    def re_assess_fdr_grouped(self):
        peptides_cat = f'{self.rescoreGroupFDRDir}/peptides_anno_MP.txt'
        proteins_cat = f'{self.rescoreGroupFDRDir}/proteins_anno_MP.txt'
        cat_pep = (f'cat {self.rescoreMPGroupdir}/peptides.txt {self.rescoreAnnoGroupDir}/peptides.txt > '
                   f'{peptides_cat}')
        os.system(cat_pep)
        cat_prot = (f'cat {self.rescoreMPGroupdir}/proteins.txt {self.rescoreAnnoGroupDir}/proteins.txt > '
                    f'{proteins_cat}')
        os.system(cat_prot)

        microproteins = {}
        records = SeqIO.parse(self.rescoreDatabase, 'fasta')
        for record in records:
            microproteins[str(record.description).replace(" ", "_")] = str(record.seq)
        perc = PercolatorPostProcessing(args=self.args)
        peptides_cat_fixed = f'{self.rescoreGroupFDRDir}/peptides_anno_MP_fixed.txt'
        perc.fix(file=peptides_cat, output=peptides_cat_fixed)

        df = pd.read_csv(peptides_cat_fixed, sep='\t')
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
            filtered_proteins = self.__protein_fdr(file=proteins_cat)
        fasta = []
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
        with open(f'{self.rescoreGroupFDRDir}/filtered_rescored_smorfs.fasta', 'w') as handler:
            handler.writelines(fasta)


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
            perc.fix(file=f'{self.rescorePostProcessDir}/{group}/psm.txt',
                     output=f'{self.rescorePostProcessDir}/{group}/psm_fixed.txt')
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
            file = f'{self.rescorePostProcessDir}/{group}/proteins.txt'
            if self.args.proteinFDR:
                filtered_proteins = self.__protein_fdr(file=file)

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

    def __protein_fdr(self, file):
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

        df = pd.read_csv(file, sep='\t')
        df = df[df["q-value"] != "q-value"]
        # df["q-value"].astype(float)
        df['q-value'] = df['q-value'].astype('float64')

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
                if '_ANNO' not in prot:
                    filtered_proteins.append(protein)
                else:
                    if self.args.keepAnnotated:
                        matches = count_substring_occurrences(prot, 'ANNO')
                        if 'ANNO' in protein and matches <= 1:
                            filtered_proteins.append(protein)
        return filtered_proteins

    def merge_results(self):
        groups = os.listdir(self.rescorePostProcessDir)
        fasta = []
        for group in groups:
            records = SeqIO.parse(f'{self.rescorePostProcessDir}/{group}/filtered_rescored_smorfs.fasta', 'fasta')
            for record in records:
                if len(str(record.seq)) <= int(self.args.maxORFLength):
                    fasta.append(f'>{str(record.description)}\n{str(record.seq)}\n')
        with open(self.rescoredMicroproteinsFasta, 'w') as outfile:
            outfile.writelines(fasta)

    def filter_gtf(self, rescored=True):
        self.print_row()
        print(f"--Filtering GTF files.")
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

    # def generate_nucleotide_fasta(self):
    #     records = SeqIO.parse(self.rescoredMicroproteinsFasta, 'fasta')
    #     for record in records:

