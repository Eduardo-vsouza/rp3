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
        self.check_dirs([self.cascadeZeroPassDir, self.cascadeZeroPassMzmlDir])
        # ZERO PASS to remove contaminant sequences
        db = self.fullContaminantsDb
        print(f"--Running Comet on {self.args.mzml} with contaminant database {db}")
        self.index_database(db=db)
        self.shower_comets(db=db, mzml_dir=self.args.mzml, pattern=self.args.fileFormat)
        self.move_pin_files(outdir=self.cascadeZeroPassDir)  # copies from self.args.mzml by default
        
        cascade = Cascade(args=self.args)
        cascade.get_zero_pass_scans()  # remove contaminant scans from self.args.mzml
        


        # FIRST PASS on reference proteome
        if self.args.splitDatabase is not None:
            dbs = self.select_database(decoy=True, proteome=True, split_db=True)
            for i, db in enumerate(dbs):
                print(f"--Db {i}/{len(dbs)}")

                self.index_database(db=db)
                print(f"--Running first-pass Comet on {self.cascadeZeroPassMzmlDir} with {db}")
                self.shower_comets(db=db, mzml_dir=self.cascadeZeroPassMzmlDir, pattern=self.args.fileFormat)
                self.move_pin_files(mzml_dir=self.cascadeZeroPassMzmlDir,  # move from zeroPass
                                    outdir=self.cascadeFirstPassDir, split_i=i)  

        else:
            db = self.select_database(decoy=True, proteome=True)
            self.index_database(db=db)
            print(f"--Running first-pass Comet on {self.args.mzml} with reference proteome")
            self.shower_comets(db=db, mzml_dir=self.cascadeZeroPassMzmlDir,
                                pattern=self.args.fileFormat)
            self.move_pin_files(mzml_dir=self.cascadeZeroPassMzmlDir,
                                outdir=self.cascadeFirstPassDir)

        # SECOND PASS on proteogenomics database
        # get scans from reference proteome that passed the first search
        cascade.get_first_pass_scans()
        self.print_row(word="Second-pass Search", color='blue')
        if self.args.splitDatabase is not None:
            dbs = self.select_database(decoy=True, proteome=False, split_db=True)
            for i, db in enumerate(dbs):
                self.print_state(message=f"Running second-pass Comet on {self.cascadeMzmlDir} with {db}",
                                 color='yellow')
                # print(f"")
                self.print_state(message=f"Db {i}/{len(dbs)}", color='yellow')
                # print(f"--Db {i}/{len(dbs)}")
                self.index_database(db=db)
                self.shower_comets(db=db, mzml_dir=self.cascadeMzmlDir, pattern=self.args.fileFormat)
                self.move_pin_files(mzml_dir=self.cascadeMzmlDir, outdir=self.cascadeSecondPassDir, split_i=i)
                self.print_state(message=f"Second-pass Comet completed for {db}", color='green')
        else:
            self.print_state(message=f"Running second-pass Comet on {self.cascadeMzmlDir} with {db}",
                    color='yellow')
            db = self.select_database(decoy=True, proteome=False)
            self.index_database(db=db)
            self.shower_comets(db=db, mzml_dir=self.cascadeMzmlDir, pattern='_filtered.mzML')
            self.move_pin_files(mzml_dir=self.cascadeMzmlDir, outdir=self.cascadeSecondPassDir) 
            self.print_state(message=f"Second-pass Comet completed for {db}", color='green')

            # print(f"--Running second-pass Comet on {self.cascadeMzmlDir} with proteogenomics database")



        cascade.concatenate_pin_files()

        if not self.args.keepIntermediate:
            self.print_state(message=f"Removing intermediate mzML files...", color='yellow')
            cmd = f'rm {self.cascadeMzmlDir}/*.mzML'
            os.system(cmd)

            cmd = f'rm {self.cascadeZeroPassMzmlDir}/*.mzML'
            os.system(cmd)
            self.print_state(message=f"Intermediate mzML files removed.", color='green')

    def index_database(self, db):
        """
        Index the database for comet search.
        """
        # db = self.select_database(decoy=True)
        if not self.args.noCometIndex:
            self.print_state(message=f"Indexing database {db} for Comet search.", color='yellow')
            # print(f"--Indexing database {db} for Comet search.")
            if not os.path.exists(f'{db}.idx') and not db.endswith(".idx"):
                cmd = f'{self.toolPaths["comet"]} -D{db} -i -P{self.params}'
                os.system(cmd)
                self.print_state(message=f"Database {db} indexed successfully.", color='green')
                # print(f"--Database {db} indexed successfully.")
            else:
                self.print_state(message=f"Database {db} already indexed. Skipping indexing step.", color='red')
                # print(f"--Database {db} already indexed. Skipping indexing step.")
        else:
            self.print_state(message=f"Skipping database indexing for Comet search as per user request.", color='red')
            # print(f"--Skipping database indexing for Comet search as per user request.")
            # if not os.path.exists(f'{db}.idx'):
            #     print(f"--Warning: Database {db} is not indexed. This may affect search performance.")

    def shower_comets(self, db, mzml_dir, pattern='.mzML'):
        """
        Iterate comet on the provided folder with mzml files.
        """
        self.print_state(message=f"Running Comet on {mzml_dir} with database {db}", color='blue')
        # print(f"--Running Comet on {mzml_dir} with database {db}")
        mzml = os.listdir(mzml_dir)
        files = ''
        for file in mzml:
            if file.endswith(pattern):
                # run = self.verify_checkpoint(outfile=os.path.join(mzml_dir, file), step="comet")
                # if run:
                files += f' {mzml_dir}/{file}'
        if not files:
            self.print_state(message=f"No mzML files found in {mzml_dir} with pattern {pattern}.", color='red')
            # print(f"--No mzML files found in {mzml_dir} with pattern {pattern}.")
        else:
            if not db.endswith(".idx"):
                if not self.args.noCometIndex:
                    database = f'{db}.idx'
                else:
                    database = db
                cmd = f'{self.toolPaths["comet"]} -D{database} -P{self.params}{files}'
                os.system(cmd)
                self.print_state(message=f"Comet search completed for {len(mzml)} files.", color='green')
    


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