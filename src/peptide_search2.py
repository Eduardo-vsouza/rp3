import os
import sys
import inspect


class MSFragger:
    def __init__(self, mzml_folder, outdir, threads, mod, quantify=False):
        self.mod = mod
        self.quantify = quantify
        self.mzMLFolder = mzml_folder # this should be organized in subfolders, each containing mzML for each group
        self.databaseDir = f'{outdir}/databases'
        self.outdir = outdir
        self.__check_dirs()
        self.MSFraggerPath = '/home/eduardo/programs/msfragger_3.5/MSFragger-3.5/MSFragger-3.5.jar'
        self.threads = threads

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
                out_db = f'{out_search}/{group}/{db}'
                if not os.path.exists(out_db):
                    os.mkdir(out_db)

    def iterate_searches(self, min_pep_len=7, max_pep_len=50):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        databases = os.listdir(self.databaseDir)
        mzml = os.listdir(self.mzMLFolder)
        for group in mzml:
            files = os.listdir(f'{self.mzMLFolder}/{group}')
            for file in files:
                fullfile = f'{self.mzMLFolder}/{group}/{file}'
                for db in databases:
                    mod = ''
                    if self.mod is not None:
                        mod = f' --variable_mod_03 {self.mod}'
                    if db.endswith(".fasta") and 'target_decoy' in db:
                        if mod is not None:
                            cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --decoy_prefix rev --output_format pin ' \
                                  f'--database_name {self.databaseDir}/{db} ' \
                                  f'--num_threads {self.threads} --digest_min_length {min_pep_len} --clip_nTerm_M 0 --allow_multiple_variable_mods_on_residue 1 ' \
                                  f'--digest_max_length {max_pep_len}{mod} {fullfile}'
                        else:
                            cmd = f'java -Xmx32g -jar {self.MSFraggerPath} --decoy_prefix rev --output_format pin ' \
                                  f'--database_name {self.databaseDir}/{db} ' \
                                  f'--num_threads {self.threads} --digest_min_length {min_pep_len} ' \
                                  f'--digest_max_length {max_pep_len} {fullfile}'
                        self.params.append(cmd)
                        os.system(cmd)
                        if self.quantify:
                            cmd_xml = f'java -Xmx32g -jar {self.MSFraggerPath} --decoy_prefix rev --output_format tsv ' \
                                  f'--database_name {self.databaseDir}/{db} --digest_min_length {min_pep_len} ' \
                                  f'--digest_max_length {max_pep_len}{mod} ' \
                                  f'--num_threads {self.threads} {fullfile}'
                            os.system(cmd_xml)
                            self.params.append(cmd_xml)
                        if fullfile.endswith("mzML"):
                            pattern = "mzML"
                        elif fullfile.endswith("mzXML"):
                            pattern = 'mzXML'
                        else:
                            pattern = 'mzML'
                        # else:
                        #     raise FileNotFoundError
                        cmd_mv_search = f'mv {fullfile.replace(f".{pattern}", ".pin",)} {self.outdir}/peptide_search/{group}/{db}/.'
                        self.params.append(cmd_mv_search)
                        if self.quantify:
                            cmd_mv_search_xml = f'mv {fullfile.replace(f".{pattern}", ".tsv",)} {self.outdir}/peptide_search/{group}/{db}/.'
                            self.params.append(cmd_mv_search_xml)
                            os.system(cmd_mv_search_xml)
                        os.system(cmd_mv_search)



