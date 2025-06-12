import os
import sys

from .basesearch import BaseSearch
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
        checked = self.check_files(files)
        params = self.__define_comet_params()
        db = self.select_database(decoy=True)
        print(f"--Running Comet on {self.args.mzml} with database {db}")

        cmd = f'{self.toolPaths["comet"]} -D{db} -P{params}{checked}'
        print(cmd)
        os.system(cmd)

        self.move_pin_files()
        return cmd


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