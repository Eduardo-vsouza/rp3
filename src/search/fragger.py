import os
import sys

from .basesearch import BaseSearch
from .cascade import Cascade


class MSFragger(BaseSearch):
    def __init__(self, args, outdir):
        """
        outdir: should be the search dir. 
        """
        super().__init__(args)
        # self.cometDir = f'{sys.path[0]}/dependencies/comet'
        self.searchOutdir = outdir
    
    def run(self):
        
        if self.args.cascade:
            self.__run_cascade()
        else:
            self.__run_standard()

    def __run_cascade(self):
        db = self.select_database(decoy=True, proteome=True)
        print(f"--Running first-pass MSFragger on {self.args.mzml} with reference proteome")
        # self.shower_comets(db=db, mzml_dir=self.args.mzml, pattern=self.args.fileFormat)
        if self.args.hlaPeptidomics:        
            cmd = self.__hla_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
        else:
            cmd = self.__std_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
        os.system(cmd)

        self.move_pin_files(outdir=self.cascadeFirstPassDir)

        # SECOND PASS on proteogenomics database
        db = self.select_database(decoy=True, proteome=False)
        print(f"--Running second-pass MSFragger on {self.cascadeMzmlDir} with proteogenomics database")

        # implement Cascade() to filter mzml here
        cascade = Cascade(args=self.args)
        # get scans from reference proteome that passed the first search
        cascade.get_first_pass_scans()
        # remove ref proteome scans from mzml files and store them in cascadeMzmlDir
        cascade.filter_mzml(mzml_dir=self.args.mzml,
                            outdir=self.cascadeMzmlDir)
        # run comet on filtered mzML; files will be stored in the same directory
        if self.args.hlaPeptidomics:        
            cmd = self.__hla_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
        else:
            cmd = self.__std_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
        os.system(cmd)
        
        self.move_pin_files(mzml_dir=self.cascadeMzmlDir, outdir=self.cascadeSecondPassDir) 
        cascade.concatenate_pin_files()

    def __run_standard(self):

        db = self.select_database(decoy=True)

        if self.args.hlaPeptidomics:        
            cmd = self.__hla_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
        else:
            cmd = self.__std_command(db=db, files=self.get_mzml(mzml_dir=self.args.mzml))
            
        os.system(cmd)
        db_relative = db.split("/")[-1]
        self.move_pin_files(mzml_dir=self.args.mzml, outdir=f'{self.outdir}/peptide_search/group/{db_relative}')



    def __hla_command(self, db, files):
        """
        return: command to run MSFragger with parameters optimized for HLA peptidomics searches.
        """
        print(f"--Running MSFragger with parameters optimized for HLA peptidomics searches")

        cmd = f'java -Xmx256g -jar {self.toolPaths["MSFragger"]} --output_format pin ' \
        f'--database_name {db} --decoy_prefix rev_ --search_enzyme_name nonspecific ' \
        f'--num_threads {self.threads} --fragment_mass_tolerance 20 --num_enzyme_termini 0 ' \
        f'--precursor_true_tolerance 20 --digest_mass_range 600.0_1500.0 --allowed_missed_cleavage_1 0 ' \
        f'--max_fragment_charge 3 --search_enzyme_cutafter ARNDCQEGHILKMFPSTWYV ' \
        f'--digest_min_length 8 --digest_max_length 12 {files}'
        print(cmd)
        # os.system(cmd)
        return cmd

    def __std_command(self, db, files):
        tmt_mod, mod, amida, pyroglu = self.__check_ptms()

        cmd = f'java -Xmx{self.args.memory}g -jar {self.toolPaths["MSFragger"]} --output_format pin ' \
            f'--database_name {db} --decoy_prefix rev ' \
            f'--num_threads {self.args.threads} --fragment_mass_tolerance {self.args.fragment_mass_tolerance} ' \
            f'--use_all_mods_in_first_search 1 --digest_min_length {self.args.digest_min_length}{tmt_mod}{mod}{amida}{pyroglu} --digest_max_length {self.args.digest_min_length} {files}'
        return cmd

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