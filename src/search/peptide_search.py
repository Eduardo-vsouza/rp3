import os
import sys
import inspect

import pandas as pd

from ..pipeline_config import PipelineStructure


class PeptideSearch(PipelineStructure):
    def __init__(self, mzml_folder, outdir, threads, mod, args, quantify=False):
        super().__init__(args=args)
        self.args = args
        self.mod = mod
        self.quantify = quantify
        self.mzMLFolder = mzml_folder # this should be organized in subfolders, each containing mzML for each group
        self.databaseDir = f'{outdir}/databases'
        self.outdir = outdir
        if self.args.groups is not None:
            self.groupsPerFile = self.read_groups(groups_df=self.args.groups)
        self.__check_dirs()
        self.MSFraggerPath = self.toolPaths["MSFragger"]
        self.threads = threads
        # It's been a while since the last time I had joy working on this shit. What am I doing with my life?
        self.mode = 'search'
        self.params = []

    def __check_dirs(self):
        out_search = f'{self.outdir}/peptide_search'
        if not os.path.exists(out_search):
            os.mkdir(out_search)
        dbs = os.listdir(self.databaseDir)
        groups = os.listdir(self.mzMLFolder)
        if not os.path.isdir(groups[0]):
            groups = [self.mzMLFolder]
            single_group = True
        else:
            single_group = False
            # groups = os.listdir(f'{self.mzMLFolder}/{group}')
        for group in groups:
            if not single_group:
                if not os.path.exists(f'{out_search}/{group}'):
                    os.mkdir(f'{out_search}/{group}')
            else:
                self.check_dirs([f'{out_search}/group'])

            for db in dbs:
                # if self.args.groups is not None:
                #     if '_'.join(db.split("_")[:2]) == group or '_'.join(db.split("_")[:1]) == group.split("_")[0] or '_'.join(db.split("_")[:1]) == group:
                #         create = True
                #     else:
                #         print(group)
                #         print(db)
                #         create = False
                # else:
                create = True
                if create:
                    out_db = f'{out_search}/{group}/{db}'
                    if single_group:
                        out_db = f'{out_search}/group/{db}'
                    if not os.path.exists(out_db):
                        os.mkdir(out_db)

    def read_groups(self, groups_df, filegroup=False):
        groups_per_file = {}
        df = pd.read_csv(groups_df, sep='\t')
        files, groups = df["file"].tolist(), df["group"].tolist()
        for file, group in zip(files, groups):
            # file_group = file.split("_")[1]
            if not file.endswith("mzML"):
                file = f'{file}.mzML'
            if filegroup:
                groups_per_file[file].append(file)
            else:
                if group not in groups_per_file:
                    groups_per_file[group] = []
                groups_per_file[group].append(file)
        return groups_per_file

    def iterate_searches(self, min_pep_len=7, max_pep_len=50):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = os.listdir(self.databaseDir)
        mod = ''
        mzml = os.listdir(self.mzMLFolder)
        if self.args.groups is not None:
            groups_per_file = self.read_groups(groups_df=self.args.groups)
        for group in mzml:
            files = os.listdir(f'{self.mzMLFolder}/{group}')
            for file in files:
                fullfile = f'{self.mzMLFolder}/{group}/{file}'
                if fullfile.endswith("mzML"):
                    pattern = "mzML"
                elif fullfile.endswith("mzXML"):
                    pattern = 'mzXML'
                elif fullfile.endswith(".d"):
                    pattern = "d"
                else:
                    pattern = 'mzML'
                for db in databases:
                    search_files = ''
                    if db.endswith(".fasta") and 'target' in db:
                        if self.args.groups is not None:
                            print(file)
                            # this_Db = '_'.join(db.split("_")[:2])
                            # print(this_Db)
                            if file in groups_per_file:
                                if group == groups_per_file[file] and '_'.join(db.split("_")[:2]) == group:
                                    run = True
                                else:
                                    run = False
                            else:
                                run = False
                        else:
                            run = True
                        if run:
                            search_files += f' {fullfile}'

                            # if mod is not None:
                            #
                            #     cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --decoy_prefix rev --output_format pin ' \
                            #           f'--database_name {self.databaseDir}/{db} ' \
                            #           f'--num_threads {self.threads} --digest_min_length {min_pep_len} --clip_nTerm_M 0 --allow_multiple_variable_mods_on_residue 1 ' \
                            #           f'--digest_max_length {max_pep_len}{mod} {fullfile}'
                            # else:
                            if self.args.comet:
                                cmd = self.__run_comet(mzml=fullfile, db=f'{self.databaseDir}/{db}')
                            else:
                                # print(file)
                                cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                                      f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                                      f'--num_threads {self.threads} --digest_min_length {min_pep_len} ' \
                                      f'--digest_max_length {max_pep_len} {fullfile}'
                            self.params.append(cmd)
                            if not os.path.exists(f'{self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_target.pin")}'):
                                # print("\n TARGET \n\n")
                                os.system(cmd)
                                cmd_mv_search = f'mv {fullfile.replace(f".{pattern}", ".pin", )} {self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_target.pin")}'
                                os.system(cmd_mv_search)
                            if self.quantify:
                                cmd_xml = f'java -Xmx32g -jar {self.MSFraggerPath} --decoy_prefix rev --output_format tsv ' \
                                      f'--database_name {self.databaseDir}/{db} --digest_min_length {min_pep_len} ' \
                                      f'--digest_max_length {max_pep_len}{mod} ' \
                                      f'--num_threads {self.threads} {fullfile}'
                                os.system(cmd_xml)
                                self.params.append(cmd_xml)

                            # else:
                            #     raise FileNotFoundError
                            if not os.path.exists(f'{self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_decoy.pin")}'):
                                # if mod is None:
                                # print("\n DECOY \n\n")
                                decoy_db = f'{self.databaseDir}/{db.replace("_target_", "_decoy_")}'
                                if self.args.comet:
                                    cmd = self.__run_comet(mzml=fullfile, db=decoy_db)
                                else:
                                    cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                                          f'--database_name {decoy_db} --decoy_prefix rev ' \
                                          f'--num_threads {self.threads} --digest_min_length {min_pep_len} ' \
                                          f'--digest_max_length {max_pep_len} {fullfile}'
                                os.system(cmd)
                                cmd_mv_decoy_search = f'mv {fullfile.replace(f".{pattern}", ".pin")} {self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_decoy.pin")}'
                                os.system(cmd_mv_decoy_search)
                            else:
                                print(f'{self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_decoy.pin")}')
                            if self.quantify:
                                cmd_mv_search_xml = f'mv {fullfile.replace(f".{pattern}", ".tsv")} {self.outdir}/peptide_search/{group}/{db}/.'
                                self.params.append(cmd_mv_search_xml)
                                os.system(cmd_mv_search_xml)
                            # os.system(cmd_mv_search)

    def __check_ptms(self):
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
        tmt_mod = ''

        mod = ''
        if self.args.mod is not None:
            mod = f' --variable_mod_0{i} {self.mod}'
            i += 1
        if self.args.tmt_mod is not None:
            tmt_mod = f' --variable_mod_03 {self.args.tmt_mod}_K_3 --variable_mod_04 {self.args.tmt_mod}_n*_3 '
        else:
            tmt_mod = ''
        return tmt_mod, mod, amida, pyroglu

    def search_files(self):
        db = self.select_database(decoy=True)
        tmt_mod, mod, amida, pyroglu = self.__check_ptms()

        if self.args.groupsFile is not None:
            groups = self.read_groups(groups_df=self.args.groupsFile)
            # print(groups)
            print(f"Metadata was provided. Files will be searched as:\n")
            print(f"Group\tFiles")
            for group in groups:
                print(f"{group}\t{', '.join(groups[group])}")
            for group in groups:
                filepaths = ''
                for file in groups[group]:
                    filepaths += f' {self.args.mzml}/{file}'
                print(filepaths)                
                if self.args.hlaPeptidomics:
                    self.__search_hla_peptidomics(db=db, search_files=filepaths, full_paths=True)
                else:
                    if self.args.engine == 'comet':
                        self.__run_comet(files=filepaths)
                    elif self.args.engine == 'msfragger':
                        cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format pin ' \
                        f'--database_name {db} --decoy_prefix rev ' \
                        f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                        f'--use_all_mods_in_first_search 1 --digest_min_length {self.args.digest_min_length}{tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {self.args.digest_min_length}{filepaths}'
                        os.system(cmd)
                db_relative = db.split("/")[-1]
                files = os.listdir(self.args.mzml)
                for file in files:
                    if file.endswith(".pin"):
                        cmd_mv = (f'mv {self.args.mzml}/{file} '
                        f'{self.outdir}/peptide_search/group/{db_relative}/{file.replace(f".pin", "_target.pin")}')
                        # os.system(cmd_mv)
                        self.exec(cmd_mv)

        else:
            files = os.listdir(self.args.mzml)

            for file in files:
                if self.args.engine == 'comet':
                    self.__run_comet(files=file)
                elif self.args.engine == 'msfragger':
                    # if self.quantify: # --output_format tsv
                    #     cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format tsv ' \
                    #         f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                    #         f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                    #         f'--use_all_mods_in_first_search 1 --digest_min_length {min_pep_len} {tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {max_pep_len}{search_files[db]}'
                    #     os.system(cmd)

                    cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format pin ' \
                            f'--database_name {db} --decoy_prefix rev ' \
                            f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                            f'--use_all_mods_in_first_search 1 --digest_min_length {self.args.digest_min_length}{tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {self.args.digest_min_length} {self.args.mzml}/{file}'
                    os.system(cmd)
                    db_relative = db.split("/")[-1]
                    cmd_mv = (f'mv {self.args.mzml}/{file.replace(f".mzML", ".pin")} '
                    f'{self.outdir}/peptide_search/group/{db_relative}/{file.replace(f".mzML", "_target.pin")}')
                    # os.system(cmd_mv)
                    self.exec(cmd_mv)

    def __run_comet(self, files):
        params = self._generate_params()
        db = self.select_database(decoy=True)
        print(f"--Running comet on {self.args.mzml} with database {db}")

        cmd = f'{self.toolPaths["comet"]} -D{db} -P{params} {files}'
        print(cmd)

        return cmd

    def __define_comet_params(self):


    def iterate_searches_cat(self, min_pep_len=7, max_pep_len=50):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = os.listdir(self.databaseDir)
        mzml = os.listdir(self.mzMLFolder)
        if self.args.groups is not None:
            groups_per_file = self.read_groups(groups_df=self.args.groups)
        # print("groups per file \n\n")
        # print(groups_per_file)

        if self.args.tmt_mod is not None:
            tmt_mod = f' --variable_mod_03 {self.args.tmt_mod}_K_3 --variable_mod_04 {self.args.tmt_mod}_n*_3 '
        else:
            tmt_mod = ''
        for group in mzml:
            search_files = {}
            database = None
            if not os.path.isdir(mzml[0]):
                files = os.listdir(self.mzMLFolder)
                single_group = True
            else:
                single_group = False
                files = os.listdir(f'{self.mzMLFolder}/{group}')

            for file in files:
                if not file.endswith(".pin"):
                    for db in databases:


                        if db.endswith(".fasta") and 'target_decoy' in db:
                            if self.args.groups is not None:
                                if file in groups_per_file:
                                    if groups_per_file[file] == db.split("_")[0]:
                                        run = True
                                    else:
                                        run = False
                                else:
                                    run = False
                            else:
                                run = True
                            if file.endswith(".xml"):
                                run = False

                            if run:
                                if single_group:
                                    fullfile = f'{self.mzMLFolder}/{file}'
                                else:
                                    fullfile = f'{self.mzMLFolder}/{group}/{file}'
                                if fullfile.endswith("mzML"):
                                    pattern = "mzML"
                                elif fullfile.endswith("mzXML"):
                                    pattern = 'mzXML'
                                elif fullfile.endswith(".d"):
                                    pattern = "d"
                                else:
                                    pattern = 'mzML'

                                database = db
                                if db not in search_files:
                                    search_files[db] = ''
                                    # search_files[decoy] = ''
                                if fullfile.endswith("mzML"):
                                    search_files[db] += f' {fullfile}'
            for db in search_files:
                splat = search_files[db].split(" ")
                outfile = splat[1]
                run = True
                i = 3

                if run:
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
                    tmt_mod = ''

                    mod = ''
                    if self.args.mod is not None:
                        mod = f' --variable_mod_0{i} {self.mod}'
                        i += 1
                    print(amida, pyroglu)
                    if self.args.hlaPeptidomics:
                        self.__search_hla_peptidomics(db=db, search_files=search_files)
                    else:
                        if self.quantify: # --output_format tsv
                            cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format tsv ' \
                              f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                              f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                              f'--use_all_mods_in_first_search 1 --digest_min_length {min_pep_len} {tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {max_pep_len}{search_files[db]}'
                            os.system(cmd)
                        cmd = f'java -Xmx{self.args.memory}g -jar {self.MSFraggerPath} --output_format pin ' \
                              f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                              f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                              f'--use_all_mods_in_first_search 1 --digest_min_length {min_pep_len}{tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {max_pep_len}{search_files[db]}'
                        print(cmd)
                        self.params.append(cmd)
                        os.system(cmd)

                    splat = search_files[db].split(" ")
                    for file in splat:
                        # print(pattern)
                        if file.endswith(pattern):
                            # print(file)
                            if not single_group:
                                if not self.args.quantifyOnly:
                                    mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                                         f'{self.outdir}/peptide_search/{group}/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                                    self.exec(mv)
                                if self.quantify:
                                    cmd_mv = (f'mv {file.replace(f".{pattern}", ".tsv")} '
                                              f'{self.outdir}/peptide_search/{group}/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.tsv")}')
                                    # os.system(cmd_mv)
                                    self.exec(cmd_mv)

                            else:
                                mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                                     f'{self.outdir}/peptide_search/group/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                                self.exec(mv)
                            # os.system(mv)
                        else:
                            print(file)
            if single_group:
                break

    def __search_hla_peptidomics(self, db, search_files, full_paths=False):
        if full_paths:
            files = search_files
            database = db
        else:
            files = search_files[db]
            database = f'{self.databaseDir}/{db}'
        cmd = f'java -Xmx256g -jar {self.toolPaths["MSFragger"]} --output_format pin ' \
              f'--database_name {database} --decoy_prefix rev_ --search_enzyme_name nonspecific ' \
              f'--num_threads {self.threads} --fragment_mass_tolerance 20 --num_enzyme_termini 0 ' \
              f'--precursor_true_tolerance 20 --digest_mass_range 600.0_1500.0 ' \
              f'--max_fragment_charge 3 --search_enzyme_cutafter ARNDCQEGHILKMFPSTWYV ' \
              f'--digest_min_length 8 --digest_max_length 12 {files}'
        # print(cmd)
        self.params.append(cmd)
        os.system(cmd)

    def iterate_searches_multi(self, min_pep_len=7, max_pep_len=50, database_dir='default'):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        if database_dir == 'default':
            database_dir = self.databaseDir
        else:
            database_dir = database_dir

        databases = os.listdir(database_dir)
        mzml = os.listdir(self.mzMLFolder)
        if self.args.groups is not None:
            groups_per_file = self.read_groups(groups_df=self.args.groups)
        for group in mzml:
            search_files = {}
            files = os.listdir(f'{self.mzMLFolder}/{group}')
            for file in files:
                # print(file)
                fullfile = f'{self.mzMLFolder}/{group}/{file}'
                if fullfile.endswith("mzML"):
                    pattern = "mzML"
                elif fullfile.endswith("mzXML"):
                    pattern = 'mzXML'
                elif fullfile.endswith(".d"):
                    pattern = "d"
                else:
                    pattern = 'mzML'
                for db in databases:
                    mod = ''
                    if self.mod is not None:
                        mod = f' --variable_mod_03 {self.mod}'
                    if db.endswith(".fasta") and 'target_decoy' in db:
                        # print(file)
                        if self.args.groups is not None:
                            if file in groups_per_file:
                                if group == groups_per_file[file] and '_'.join(db.split("_")[:2]) == group:
                                    run = True
                                else:
                                    run = False
                            else:
                                run = False
                        else:
                            run = True
                        if run:
                            if db not in search_files:
                                search_files[db] = ''
                            search_files[db] += f' {fullfile}'

            for db in search_files:
                cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                      f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                      f'--num_threads {self.threads} --digest_min_length {min_pep_len} ' \
                      f'--digest_max_length {max_pep_len}{search_files[db]}'
                self.params.append(cmd)
                # print("\n TARGET \n\n")
                for file in search_files[db].split(" "):
                    if file.endswith(pattern):
                        mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                             f'{self.outdir}/peptide_search/{group}/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                        os.system(mv)

    # def __run_comet(self, mzml, db):
    #     cmd = f'{self.args.comet_path} -D{db} -P{self.args.comet_params} {mzml}'
    #     print(cmd)
    #     return cmd
