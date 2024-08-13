import os
import sys
import json
import inspect

import pandas as pd
from Bio import SeqIO

from .utils import group_folder_generator, check_multiple_dirs
from .pipeline_config import PipelineStructure


class MOFF(PipelineStructure):
    def __init__(self, args, outdir):
        super().__init__(args=args)
        self.args = args

        # directories
        self.outdir = outdir  # this is redefined so we can re-use the class on the compare mode
        self.searchDir = f'{self.outdir}/peptide_search'
        self.quantificationDir = f'{self.outdir}/quantification'
        self.moffOutdir = f'{self.outdir}/moFF'
        self.postProcessDir = f'{self.outdir}/post_processing'
        self.summarizedDir = f'{self.outdir}/summarized_results'
        self.flashLFQDir = f'{self.outdir}/flashLFQ'

        # for running moFF
        self.databases = {}
        # self.flashLFQPath = self.args.flash_lfq_path

        self.__check_dirs()

        # params
        self.mode = 'quant'
        self.params = []
        # print(self.outdir)

    def __check_dirs(self):
        folders = [self.quantificationDir, self.moffOutdir]
        for folder in folders:
            if not os.path.exists(folder):
                os.mkdir(folder)
        if not os.path.exists(self.searchDir):
            print("Peptide search folder not found. Please run this pipeline_config in 'search' mode.")

    @staticmethod
    def __check_dir(folder):
        if not os.path.exists(folder):
            os.mkdir(folder)

    def get_fdr_peptides(self):
        print(f"--Gathering peptides")
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        # gen = group_folder_generator(self.postProcessDir)
        proteins_by_file = {}
        peptides = self.select_peptides_df()
        # for content in gen:
        entries = []
        fasta = self.select_fasta()
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            entries.append(str(record.description))
        # peptides = f'{content.dbDir}/peptides_fixed.txt'
        df = pd.read_csv(peptides, sep='\t')
        df = df[df["q-value"] <= 0.01]
        df = df[~df["proteinIds"].str.contains(",")]
        df = df[(df["proteinIds"].str.contains("rev_") == False) & (df["proteinIds"].str.contains("contaminant") == False)]
        if self.args.no_anno:
            df = df[df["proteinIds"].str.contains("_ANNO") == False]
            df = df[df["proteinIds"].isin(entries)]
        proteins = df["proteinIds"].tolist()
        group_db = 'group'
        if group_db not in proteins_by_file:
            proteins_by_file[group_db] = []

        for prot in proteins:
            prot_list = prot.split(",")
            for protein in prot_list:
                # if not protein.startswith("")
                proteins_by_file[group_db].append(protein)
        with open(f'{self.summarizedDir}/proteins_by_group_db.json', 'w') as outfile:
            json.dump(proteins_by_file, outfile)

    def generate_moff_input(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        groups = os.listdir(self.searchDir)
        with open(f'{self.summarizedDir}/proteins_by_group_db.json') as json_file:
            proteins_by_group = json.load(json_file)

            for group in groups:
                group_dir = f'{self.searchDir}/{group}'
                databases = os.listdir(group_dir)
                for db in databases:
                    if db.endswith("target_decoy_database.fasta") or db.endswith("_target_decoy.fasta"):
                        self.__check_dir(folder=f'{self.quantificationDir}/{group}')
                        self.__check_dir(folder=f'{self.quantificationDir}/{group}/{db}')

                        db_dir = f'{group_dir}/{db}'
                        files = os.listdir(db_dir)
                        for file in files:
                            print(file)
                            if file.endswith(".tsv"):
                                df = self.convert_tsv_to_moff(file=f'{db_dir}/{file}')
                                group_db = f'{group}/{db}'
                                proteins = proteins_by_group[group_db]
                                df = df[df["prot"].isin(proteins)]
                                df.to_csv(f'{self.quantificationDir}//{group}/{db}/{file}_moff_input.tsv', sep='\t', index=False)

    def generate_flash_lfq_input(self):
        print(f"--Generating FlashLFQ input")
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        rescored = self.is_rescored()
        search_dir = self.select_search_dir()
        groups = os.listdir(search_dir)

        with open(f'{self.summarizedDir}/proteins_by_group_db.json') as json_file:
            proteins_by_group = json.load(json_file)

            # for group in groups:
            group_dir = f'{search_dir}/group'

            databases = os.listdir(group_dir)
            df = pd.DataFrame()
            # if db.endswith("target_decoy_database.fasta") or db.endswith("_target_decoy.fasta"):
            generate = True
            self.__check_dir(folder=f'{self.quantificationDir}/group')
            # self.__check_dir(folder=f'{self.quantificationDir}/group/{db}')

            db_dir = f'{group_dir}'
            files = os.listdir(db_dir)
            for file in files:
                # print(file)
                if file.endswith(".tsv") and 'spectraRT' not in file:
                    ndf = self.convert_tsv_to_flash_lfq(file=f'{db_dir}/{file}')
                    df = pd.concat([df, ndf])
                    # group_db = f'{group}/{db}'
                    group_db = 'group'
                    proteins = proteins_by_group[group_db]
                    df = df[df["Protein Accession"].isin(proteins)]
            if generate:
                df.to_csv(f'{self.quantificationDir}/flash_lfq_input.tsv', sep='\t', index=False)

    def convert_tsv_to_moff(self, file):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        df = pd.read_csv(file, sep='\t')
        # df = df[['a', 'b']]
        df = df[['peptide', 'protein', 'modification_info', 'retention_time', 'precursor_neutral_mass',
                 'calc_neutral_pep_mass', 'charge']]
        mz = self.__get_mz(df)
        mods = self.__reformat_modifications(df)
        df = df.drop(columns=["precursor_neutral_mass", "modification_info"])

        df.insert(2, "mod_peptide", mods)
        df.insert(3, "mz", mz)
        data = {"peptide": df["peptide"].tolist(), 'prot': df["protein"].tolist(),
                'mod_peptide': df["mod_peptide"].tolist(), 'rt': df["retention_time"].tolist(), 'mz': df["mz"].tolist(),
                "mass": df["calc_neutral_pep_mass"].tolist(), 'charge': df["charge"].tolist()}
        ndf = pd.DataFrame(data=data)


        return ndf

    def __add_mzml_files_to_df(self, df, file_name):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        # mzml = f'{os.path.abspath(self.args.mzml)}/{file_name.replace(".tsv", ".mzML").replace("_target", "")}'
        mzml = f'{file_name.replace(".tsv", ".mzML").replace("_target", "")}'
        mzml_column = []
        length = len(df["peptide"].tolist())
        for i in range(length):
            mzml_column.append(mzml)
        df.insert(2, "File Name", mzml_column)
        return df


    def convert_tsv_to_flash_lfq(self, file):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        df = pd.read_csv(file, sep='\t')
        # print(df.columns)
        print('└──', file)
        # df = df[['a', 'b']]
        # if 'modification_info' in df.columns:
        df = df[['peptide', 'protein', 'modification_info', 'retention_time', 'calc_neutral_pep_mass', 'charge']]
        mods = self.__reformat_modifications_flash_lfq(df)
        proteins = self.__reformat_protein_accession(df)
        df = df.drop(columns=['protein'])
        df.insert(1, "protein", proteins)
        df.insert(2, "mod_peptide", mods)
        #df = df.drop(columns=["precursor_neutral_mass", "modification_info"])
        df = df.drop(columns=["modification_info"])
        df = self.__add_mzml_files_to_df(df, file.split("/")[-1])
        #df.insert(3, "mz", mz)
        data = {"File Name" : df["File Name"].tolist(), "Base Sequence": df["peptide"].tolist(),
                'Full Sequence': df["mod_peptide"].tolist(),
                "Peptide Monoisotopic Mass": df["calc_neutral_pep_mass"].tolist(),
                'Scan Retention Time': df["retention_time"].tolist(),
                'Precursor Charge': df["charge"].tolist(), "Protein Accession": df["protein"].tolist()}
        ndf = pd.DataFrame(data=data)


        return ndf

    @staticmethod
    def __reformat_protein_accession(df):
        fixed = []
        proteins = df["protein"].tolist()
        for protein in proteins:
            fixed.append(protein.split(" ")[0])
        return fixed

    @staticmethod
    def __reformat_modifications(df):
        """
        This is for moFF input. See __reformat_modifications_flashlfq for FlashLFQ input.
        {
        "<cmm>": {"deltaChem":[3,2,1,1],"desc":"Carboxyamidomethylation C  unimod:4"},
        "<ox>": {"deltaChem":[0,0,0,1],"desc":"oxidation oxidation unimod:35" } ,
        "ace-":  {"deltaChem":[2,2,0,1],"desc":"Acetylation N-term unimod:1" },
        }
        :param df:
        :return:
        """
        moff_modifications = {'57.021465': '<cmm>', '15.9949': '<ox>', '42.0106': 'ace-'}

        mods = df["modification_info"].tolist()
        peps = df["peptide"].tolist()
        reformatted_mods = []
        for mod, pep in zip(mods, peps):
            modified_peptide = ''
            if type(mod) != float:
                mod_list = mod.split(", ")
                peptide_raw = [aa for aa in pep]
                # print("mod", mod)
                # print("mod_list", mod_list)
                # print("pep", pep)
                # print("pep_raw", peptide_raw)
                number_of_mods = 0  # this is important because everytime we add <ox>, for instance, the index
                                    # increases by 1, thus changing the remaining indexes.
                for ptm in mod_list:
                    mass = ptm.split("(")[-1][:-1]
                    # print("ptm", ptm)
                    # print("mass", mass)
                    modified_aa = ptm.split("(")[0]
                    if modified_aa.startswith("N-term"):
                        pos = 0
                    else:
                        pos = int(''.join([i for i in modified_aa if i.isdigit()]))
                    # print("pos", pos)
                    moff_mod = moff_modifications[mass]
                    # print("moff_mod", moff_mod)
                    peptide_raw.insert(pos+number_of_mods, moff_mod)
                    number_of_mods += 1
                modified_peptide = f'NH2-{"".join(peptide_raw)}-COOH'
                if 'ace-' in modified_peptide:
                    modified_peptide = modified_peptide.replace("NH2-", "")
            reformatted_mods.append(modified_peptide)
                # print("final_mod_pep", modified_peptide)
                # print("\n")
        return reformatted_mods

    @staticmethod
    def __reformat_modifications_flash_lfq(df):
        """

        {
        "<cmm>": {"deltaChem":[3,2,1,1],"desc":"Carboxyamidomethylation C  unimod:4"},
        "<ox>": {"deltaChem":[0,0,0,1],"desc":"oxidation oxidation unimod:35" } ,
        "ace-":  {"deltaChem":[2,2,0,1],"desc":"Acetylation N-term unimod:1" },
        }
        :param df:
        :return:
        """
        moff_modifications = {'57.021465': '[Carbamidomethylation]', '15.9949': '[Oxidation]',
                              '42.0106': '[Acetylation]'}

        mods = df["modification_info"].tolist()
        peps = df["peptide"].tolist()
        reformatted_mods = []
        for mod, pep in zip(mods, peps):
            modified_peptide = ''
            if type(mod) != float:
                mod_list = mod.split(", ")
                peptide_raw = [aa for aa in pep]
                # print("mod", mod)
                # print("mod_list", mod_list)
                # print("pep", pep)
                # print("pep_raw", peptide_raw)
                number_of_mods = 0  # this is important because everytime we add <ox>, for instance, the index
                # increases by 1, thus changing the remaining indexes.
                for ptm in mod_list:
                    mass = ptm.split("(")[-1][:-1]
                    # print("ptm", ptm)
                    # print("mass", mass)
                    modified_aa = ptm.split("(")[0]
                    if modified_aa.startswith("N-term"):
                        pos = 0
                    else:
                        pos = int(''.join([i for i in modified_aa if i.isdigit()]))
                    # print("pos", pos)
                    moff_mod = moff_modifications[mass]
                    # print("moff_mod", moff_mod)
                    peptide_raw.insert(pos + number_of_mods, moff_mod)
                    number_of_mods += 1
                modified_peptide = f'{"".join(peptide_raw)}'
                if 'ace-' in modified_peptide:
                    modified_peptide = modified_peptide.replace("NH2-", "")
            if modified_peptide == '':
                modified_peptide = pep
            reformatted_mods.append(modified_peptide)
            # print("final_mod_pep", modified_peptide)
            # print("\n")
        return reformatted_mods

    @staticmethod
    def __get_mz(df):
        mz = []
        masses = df["precursor_neutral_mass"].tolist()
        charges = df["charge"].tolist()
        for mass, charge in zip(masses, charges):
            mz.append(float(mass)/float(charge))
        return mz


    def separate_by_database(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = {}
        groups = os.listdir(self.quantificationDir)
        for group in groups:
            groupdir = f'{self.quantificationDir}/{group}'
            dbs = os.listdir(groupdir)
            for db in dbs:
                if db not in databases:
                    databases[db] = {}
                dbdir = f'{groupdir}/{db}'
                files = os.listdir(dbdir)
                for file in files:
                    mzml = f'{self.args.mzml}/{group}/{file.split(".")[0]}.mzML'
                    databases[db][f'{dbdir}/{file}'] = mzml
        self.databases = databases

    def run_moff(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        for db in self.databases:
            moff_db_dir = f'{self.moffOutdir}/{db}'
            if not os.path.exists(moff_db_dir):
                os.mkdir(moff_db_dir)

            tsv_list = []
            mzml_list = []
            for tsv in self.databases[db]:
                print(os.path.getsize(tsv))
                if os.path.getsize(tsv) > 100:
                    if '20120322_EXQ1_MiBa_SA_HCC1143_1_A.tsv_moff_input.tsv' not in tsv and '20120322_EXQ1_MiBa_SA_HCC1937_1_A.tsv_moff_input.tsv' not in tsv:
                        if '20120322_EXQ1_MiBa_SA_HCC1937_2_A.tsv_moff_input.tsv' not in tsv:
                            tsv_list.append(tsv)
                            mzml_list.append(self.databases[db][tsv])


            cmd = f'/home/eduardo/programs/moFF/moff_env_py3.6/bin/python3 {self.args.moff_path} --tsv_list {" ".join(tsv_list)} --raw_list {" ".join(mzml_list)} ' \
                  f'--tol 20 --loc_out {moff_db_dir} --peptide_summary'
            # print(cmd)
            self.params.append(cmd)
            os.system(cmd)

    def run_flashlfq(self):
        print(f"--Running FlashLFQ on protein identifications")
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        # gen = group_folder_generator(self.quantificationDir)
        # check_multiple_dirs([self.flashLFQDir])
        # for content in gen:

        # group_dir = f'{self.flashLFQDir}/{content.group}'
        # db_dir = f'{group_dir}/{content.db}'
        # check_multiple_dirs([group_dir, db_dir])
        # if content.file.endswith(".tsv"):
        file = f'{self.quantificationDir}/flash_lfq_input.tsv'
        cmd = f'{self.flashLFQPath} --idt {file} --ppm 20 --rep {self.args.mzml} ' \
              f'--out {self.quantificationDir} --mbr true --thr {self.args.threads}'
        self.params.append(cmd)
        os.system(cmd)

    def count_spectra(self):
        dbs = os.listdir(self.summarizedDir)
        for db in dbs:
            if db != 'merged':
                db_dir = f'{self.summarizedDir}/{db}'
                mps = f'{db_dir}/microproteins_150.fasta'

                counts = {}

                records = SeqIO.parse(mps, 'fasta')
                for record in records:
                    counts[str(record.description)] = {'counts': 0, 'samples': [], 'n_replicates': 0}
                print(counts)
                perc_groups = os.listdir(self.postProcessDir)
                for group in perc_groups:
                    group_dir = f'{self.postProcessDir}/{group}'
                    perc = f'{group_dir}/{db}/peptides_filtered.txt'
                    df = pd.read_csv(perc, sep='\t')
                    print(list(counts.keys()))
                    df = df[df['proteinIds'].str.contains('|'.join(counts.keys()))]
                    print(df)

                    prots = df["proteinIds"].tolist()
                    files = df["PSMId"].tolist()
                    for prot_list, file in zip(prots, files):
                        # if ',' in prot:
                        proteins = prot_list.split(",")
                        for protein in proteins:
                            if protein in counts:
                                counts[protein]['counts'] += 1
                                file = file.split(".")[0]
                                if file not in counts[protein]['samples']:
                                    counts[protein]['samples'].append(file)
                                    counts[protein]['n_replicates'] = len(counts[protein]['samples'])
                data = {'protein': [], 'spectral_counts': [], 'samples': [], 'n_replicates': []}
                for smorf in counts:
                    data['protein'].append(smorf)
                    data['spectral_counts'].append(counts[smorf]['counts'])
                    data['samples'].append(','.join(counts[smorf]['samples']))
                    data['n_replicates'].append(counts[smorf]['n_replicates'])
                df = pd.DataFrame(data=data)
                df.to_csv(f'{db_dir}/{db.replace(".fasta", "")}_spectral_counts.txt', sep='\t', index=False)


    def count_spectra_annotated(self):
        # groups = os.listdir(self.postProcessDir)
        genn = group_folder_generator(self.postProcessDir)
        counts = {}
        for content in genn:
            if content.file == 'peptides.txt':
                self.__fix_peptides_anno(file=content.fullFile, output=f"{content.dbDir}/peptides_fixed_anno.txt")
                df = pd.read_csv(f'{content.dbDir}/peptides_fixed_anno.txt', sep="\t")
                df = df[df["q-value"] != "q-value"]
                df = df[df["q-value"] < 0.01]
                df = df[df["proteinIds"].str.contains(",") == False]
                proteins = df["proteinIds"].tolist()
                files = df["PSMId"].tolist()
                for prot_list, file in zip(proteins, files):
                    prots = prot_list.split(",")
                    sample = file.split(".")[0]
                    for prot in prots:
                        if prot not in counts:
                            counts[prot] = {"counts": 0, "samples": []}
                        counts[prot]["counts"] += 1
                        if sample not in counts[prot]["samples"]:
                            counts[prot]["samples"].append(sample)
                smorfs = self.filter_annotated_smorfs(db=f"{self.outdir}/databases/{content.db}")

        data = {"protein": [], "spectral_counts": [], "samples": [], "n_replicates": []}
        for protein in counts:
            data["protein"].append(protein)
            data["spectral_counts"].append(counts[protein]["counts"])
            data["samples"].append(",".join(counts[protein]["samples"]))
            data["n_replicates"].append(len(counts[protein]["samples"]))
        df = pd.DataFrame(data=data)
        df.to_csv(f"{self.quantificationDir}/annotated_spec_counts.txt", sep="\t", index=False)
        summ_df = df.drop(columns=["samples"])
        summ_df.to_csv(f"{self.quantificationDir}/annotated_spec_counts_summarized.txt", sep="\t", index=False)
        smorfs_df = summ_df[summ_df["protein"].isin(smorfs)]
        smorfs_df.to_csv(f"{self.quantificationDir}/annotated_smorfs_spec_counts_summarized.txt", sep="\t", index=False)

    def __fix_peptides_anno(self, file, output):
        with open(file, 'r') as handler, open(output, 'w') as outfile:
            lines = handler.readlines()
            fixed = []
            for line in lines:
                if not line.startswith("PSMId"):
                    cols = line.split("\t")
                    prots = cols[5:]
                    prot_list = ''
                    for prot in prots:
                        prot_list += f'{prot.replace(",", "")},'
                    first_cols = "\t".join(cols[:5])
                    newline = f'{first_cols}\t{prot_list[:-1]}\n'
                    fixed.append(newline)
                else:
                    fixed.append(line)
            outfile.writelines(fixed)

    @staticmethod
    def filter_annotated_smorfs(db):
        records = SeqIO.parse(db, 'fasta')
        smorfs = []
        for record in records:
            if "ANNO" in str(record.description):
                if len(str(record.seq)) <= 150:
                    smorfs.append(str(record.description))
        return smorfs

