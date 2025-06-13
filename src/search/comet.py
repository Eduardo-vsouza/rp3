import os
import sys

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

    def run(self, files):
        """
        Step 1 of cascade search.
        """
        checked = self.check_files(files)
        if self.args.cascade:
            self._run_cascade(files=checked)
        else:
            self._run_standard(files=checked)

    def _run_standard(self, files):
        """
        step 1.5, optional
        """
        params = self.__define_comet_params()
        db = self.select_database(decoy=True)
        print(f"--Running Comet on {self.args.mzml} with database {db}")

        cmd = f'{self.toolPaths["comet"]} -D{db} -P{params}{files}'
        print(cmd)
        os.system(cmd)

        self.move_pin_files()

    def _run_cascade(self, files):
        """
        Step 1.5, optional
        Private function
        To be called by self.run()
        """
        self.params = self.__define_comet_params()
        db = self.select_database(decoy=True, proteome=True)

        # FIRST PASS on reference proteome
        print(f"--Running first-pass Comet on {self.args.mzml} with reference proteome")
        cmd = f'{self.toolPaths["comet"]} -D{db} -P{self.params}{files}'
        print(cmd)
        os.system(cmd)
        self.move_pin_files(outdir='cascade_first')

        # SECOND PASS on proteogenomics database
        db = self.select_database(decoy=True)
        print(f"--Running second-pass Comet on {self.cascadeMzmlDir} with proteogenomics database")

        # implement Cascade() to filter mzml here
        cascade = Cascade(args=self.args)
        cascade.get_first_pass_scans(cascade_dir=self.cascadeFirstPassDir)
        cascade.filter_mzml(mzml_dir=self.cascadeMzmlDir)

        self.shower_comets(db=db, mzml_dir=self.cascadeMzmlDir)

        # mzml = os.listdir(self.cascadeMzmlDir)
        # files = ''
        # for file in mzml:
        #     if file.endswith('_filtered.mzML'):
        #         files += f' {self.cascadeMzmlDir}/{file}'
        # cmd = f'{self.toolPaths["comet"]} -D{db} -P{params}{files}'
        # os.system(cmd)

        self.move_pin_files(outdir='standard')

    def shower_comets(self, db, params, mzml_dir, pattern='.mzML'):
        """
        Iterate comet on the provided folder with mzml files.
        """
        mzml = os.listdir(mzml_dir)
        files = ''
        for file in mzml:
            if file.endswith(pattern):
                files += f' {mzml_dir}/{file}'
        cmd = f'{self.toolPaths["comet"]} -D{db} -P{params}{files}'
        os.system(cmd)

    def __remove_reference_scans(self):
        """
        Step 1.6, following the end of both cascade searches (reference proteome and proteogenomics database)
        Remove scans from the first-pass search that were already identified in the reference proteome.
        """
        


    def concatenate_pin_files(self):
        # cmd = f'cat {self.cascade}'
        ...

    def __define_comet_params(self):
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