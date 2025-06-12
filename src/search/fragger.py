import os
import sys

from .basesearch import BaseSearch


class MSFragger(BaseSearch):
    def __init__(self, args, outdir):
        """
        outdir: should be the search dir. 
        """
        super().__init__(args)
        # self.cometDir = f'{sys.path[0]}/dependencies/comet'
        self.searchOutdir = outdir
    
    def run(self, filepaths):
        tmt_mod, mod, amida, pyroglu = self.__check_ptms()

        db = self.select_database(decoy=True)
        cmd = f'java -Xmx{self.args.memory}g -jar {self.toolPaths["MSFragger"]} --output_format pin ' \
        f'--database_name {db} --decoy_prefix rev ' \
        f'--num_threads {self.args.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
        f'--use_all_mods_in_first_search 1 --digest_min_length {self.args.digest_min_length}{tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {self.args.digest_min_length} {filepaths}'
        os.system(cmd)
        db_relative = db.split("/")[-1]
        files = os.listdir(self.args.mzml)
        for file in files:
            if file.endswith(".pin"):
                cmd_mv = (f'mv {self.args.mzml}/{file} '
                f'{self.outdir}/peptide_search/group/{db_relative}/{file.replace(f".pin", "_target.pin")}')
                os.system(cmd_mv)
                # self.exec(cmd_mv)


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