import os
import sys
import inspect

import pandas as pd


class MSFragger:
    def __init__(self, mzml_folder, outdir, threads, mod, args, quantify=False):
        self.args = args
        self.mod = mod
        self.quantify = quantify
        self.mzMLFolder = mzml_folder # this should be organized in subfolders, each containing mzML for each group
        self.databaseDir = f'{outdir}/databases'
        self.outdir = outdir
        if self.args.groups is not None:
            self.groupsPerFile = self.read_groups(groups_df=self.args.groups)
        self.__check_dirs()
        self.MSFraggerPath = '/home/microway/programs/msfragger_3.5/MSFragger-3.5/MSFragger-3.5.jar'
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
        for group in groups:
            if not os.path.exists(f'{out_search}/{group}'):
                os.mkdir(f'{out_search}/{group}')

            for db in dbs:
                if self.args.groups is not None:
                    if '_'.join(db.split("_")[:2]) == group:
                        create = True
                    else:
                        create = False
                else:
                    create = True
                if create:
                    out_db = f'{out_search}/{group}/{db}'
                    if not os.path.exists(out_db):
                        os.mkdir(out_db)

    def read_groups(self, groups_df):
        groups_per_file = {}
        df = pd.read_csv(groups_df, sep='\t')
        files, groups = df["files"].tolist(), df["groups"].tolist()
        for file, group in zip(files, groups):
            # file_group = file.split("_")[1]
            groups_per_file[file] = group
        return groups_per_file

    def iterate_searches(self, min_pep_len=7, max_pep_len=50):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = os.listdir(self.databaseDir)
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
                    mod = ''
                    if self.mod is not None:
                        mod = f' --variable_mod_03 {self.mod}'
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
                                print(file)
                                cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                                      f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                                      f'--num_threads {self.threads} --digest_min_length {min_pep_len} ' \
                                      f'--digest_max_length {max_pep_len} {fullfile}'
                            self.params.append(cmd)
                            if not os.path.exists(f'{self.outdir}/peptide_search/{group}/{db}/{file.replace(f".{pattern}", "_target.pin")}'):
                                print("\n TARGET \n\n")
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
                                print("\n DECOY \n\n")
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

    def iterate_searches_cat(self, min_pep_len=7, max_pep_len=50):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = os.listdir(self.databaseDir)
        mzml = os.listdir(self.mzMLFolder)
        if self.args.groups is not None:
            groups_per_file = self.read_groups(groups_df=self.args.groups)
        if self.args.tmt_mod is not None:
            tmt_mod = f'--variable_mod_03 {self.args.tmt_mod}_K_3 --variable_mod_04 {self.args.tmt_mod}_n*_3 '
        else:
            tmt_mod = ''
        for group in mzml:
            search_files = {}
            database = None
            files = os.listdir(f'{self.mzMLFolder}/{group}')
            for file in files:
                for db in databases:
                    mod = ''
                    if self.mod is not None:
                        mod = f' --variable_mod_03 {self.mod}'

                    if db.endswith(".fasta") and 'target_decoy' in db:

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
                            fullfile = f'{self.mzMLFolder}/{group}/{file}'
                            if fullfile.endswith("mzML"):
                                pattern = "mzML"
                            elif fullfile.endswith("mzXML"):
                                pattern = 'mzXML'
                            elif fullfile.endswith(".d"):
                                pattern = "d"
                            else:
                                pattern = 'mzML'
                            # print(file)

                            database = db
                            # decoy = db.replace("_target_", "_decoy_")
                            if db not in search_files:
                                search_files[db] = ''
                                # search_files[decoy] = ''
                            search_files[db] += f' {fullfile}'
                            # search_files[decoy] += f' {fullfile}'

            for db in search_files:
                cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --output_format pin ' \
                      f'--database_name {self.databaseDir}/{db} --decoy_prefix rev ' \
                      f'--num_threads {self.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
                      f'--digest_min_length {min_pep_len} {tmt_mod} --digest_max_length {max_pep_len}{search_files[db]}'
                self.params.append(cmd)
                os.system(cmd)
                print("\n TARGET \n\n")
                print(db)
                print(search_files)
                print(group)
                splat = search_files[db].split(" ")
                for file in splat:
                    print(pattern)
                    if file.endswith(pattern):
                        print(file)
                        mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                             f'{self.outdir}/peptide_search/{group}/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                        os.system(mv)
                    else:
                        print(file)

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
                print(file)
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
                        print(file)
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
                print("\n TARGET \n\n")
                for file in search_files[db].split(" "):
                    if file.endswith(pattern):
                        mv = f'mv {file.replace(f".{pattern}", ".pin")} ' \
                             f'{self.outdir}/peptide_search/{group}/{db}/{file.split("/")[-1].replace(f".{pattern}", "_target.pin")}'
                        os.system(mv)

    def __run_comet(self, mzml, db):
        cmd = f'{self.args.comet_path} -D{db} -P{self.args.comet_params} {mzml}'
        print(cmd)
        return cmd
