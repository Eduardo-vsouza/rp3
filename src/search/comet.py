import os
import sys

import pandas as pd

from .basesearch import BaseSearch
from .cascade import Cascade
from ..pipeline_config import PipelineStructure


class Comet(BaseSearch):
    def __init__(self, args, outdir):
        """
        outdir: should be the search dir. 
        """
        super().__init__(args)
        self.cometDir = f'{sys.path[0]}/dependencies/comet'
        self.searchOutdir = outdir

    def run(self):
        """
        Step 1 of cascade search.
        """
        self.params = self.__define_comet_params()

        if self.args.cascade:
            self._run_cascade()
        else:
            self._run_standard()

    def _run_standard(self):
        """
        step 1.5, optional
        """
        # params = self.__define_comet_params()
        db = self.select_database(decoy=True)
        print(f"--Running Comet on {self.args.mzml} with database {db}")
        self.shower_comets(db=db, mzml_dir=self.args.mzml, pattern=self.args.fileFormat)
        self.move_pin_files(outdir='standard')

    def _run_cascade(self):
        """
        Step 1.5, optional
        Private function
        To be called by self.run()
        """

        # FIRST PASS on reference proteome
        db = self.select_database(decoy=True, proteome=True)
        print(f"--Running first-pass Comet on {self.args.mzml} with reference proteome")
        self.shower_comets(db=db, mzml_dir=self.args.mzml, pattern=self.args.fileFormat)
        self.move_pin_files(outdir=self.cascadeFirstPassDir)

        # SECOND PASS on proteogenomics database
        db = self.select_database(decoy=True, proteome=False)
        print(f"--Running second-pass Comet on {self.cascadeMzmlDir} with proteogenomics database")

        # implement Cascade() to filter mzml here
        cascade = Cascade(args=self.args)
        # get scans from reference proteome that passed the first search
        cascade.get_first_pass_scans()
        # remove ref proteome scans from mzml files and store them in cascadeMzmlDir
        cascade.filter_mzml(mzml_dir=self.args.mzml,
                            outdir=self.cascadeMzmlDir)
        # run comet on filtered mzML; files will be stored in the same directory
        self.shower_comets(db=db, mzml_dir=self.cascadeMzmlDir, pattern='_filtered.mzML')

        self.move_pin_files(mzml_dir=self.cascadeMzmlDir, outdir=self.cascadeSecondPassDir) 
        self.concatenate_pin_files()

    def shower_comets(self, db, mzml_dir, pattern='.mzML'):
        """
        Iterate comet on the provided folder with mzml files.
        """
        print(f"--Running Comet on {mzml_dir} with database {db}")
        mzml = os.listdir(mzml_dir)
        files = ''
        for file in mzml:
            if file.endswith(pattern):
                # run = self.verify_checkpoint(outfile=os.path.join(mzml_dir, file), step="comet")
                # if run:
                files += f' {mzml_dir}/{file}'
        if not files:
            print(f"--No mzML files found in {mzml_dir} with pattern {pattern}.")
        else:
            cmd = f'{self.toolPaths["comet"]} -D{db} -P{self.params}{files}'
            os.system(cmd)
       
    def concatenate_pin_files(self):
        print(f"--Concatenating pin files from first and second pass of cascade search...")
        cmd = f'cat {self.cascadeFirstPassDir}/*pin {self.cascadeSecondPassDir}/*pin > {self.searchOutdir}/cascade_search_unfixed.pin'
        os.system(cmd)

        cmd = f"grep -v 'SpecId' {self.searchOutdir}/cascade_search_unfixed.pin > {self.searchOutdir}/cascade_search_filtered.pin"
        os.system(cmd)

        cmd = f"awk 'FNR<2' {self.searchOutdir}/cascade_search_unfixed.pin | cat - {self.searchOutdir}/cascade_search_filtered.pin > {self.searchOutdir}/cascade_search.pin"
        os.system(cmd)


    def __define_comet_params(self):
        print(f"--Defining Comet params...")
        if self.args.highRes:
            params = f'{self.cometDir}/cometParams_highRes_rp3.txt'
        elif self.args.lowRes:
            params = f'{self.cometDir}/cometParams_lowRes_rp3.txt'
        elif self.args.hlaPeptidomics:
            params = f'{self.cometDir}/cometParams_hla_rp3.txt'
        else:
            params = f'{self.cometDir}/cometParams.txt'
        cmd = f'cp {params} {self.searchOutdir}/comet_params.txt'
        os.system(cmd)
        return params