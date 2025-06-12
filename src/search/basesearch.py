import os

from ..pipeline_config import PipelineStructure


class BaseSearch(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
    
        self.args = args
        self.mod = self.args.mod
        self.quantify = self.args.quantify
        self.mzMLFolder = self.args.mzml
        self.databaseDir = self.databaseDir
        self.outdir = self.args.outdir
        if self.args.groups is not None:
            self.groupsPerFile = self.read_groups(groups_df=self.args.groups)
        # self.__check_dirs()
        self.threads = self.args.threads
        self.mode = 'search'
        self.params = []

        if self.args.cascade:
            self.check_dirs([self.cascadeDir, self.cascadeMzmlDir, self.cascadePinDir])

    def check_files(self, files):
        checked = ''
        filelist = files.split(" ")
        for file in filelist:
            if file.endswith(self.args.fileFormat):
                checked += f' {file}'
        return checked


    def move_pin_files(self, outdir="standard"):
        # outdirs = {"standard": }
        if outdir == 'standard':
            db_relative = self.select_database(decoy=True).split("/")[-1]
            output_dir = f'{self.searchDir}/group/{db_relative}'
        elif outdir == 'rescore':
            output_dir = f'{self.outdir}/rescore/peptide_search/group'
        elif outdir == 'cascade':
            output_dir = self.cascadeDir

        files = os.listdir(self.args.mzml)
        for file in files:
            if file.endswith(".pin"):
                cmd_mv = (f'mv {self.args.mzml}/{file} '
                          f'{output_dir}/{file.replace(".pin", "_target.pin")}')
                os.system(cmd_mv)
            elif file.endswith(".txt") or file.endswith(".xml"):
                cmd_mv = (f'mv {self.args.mzml}/{file} {output_dir}/{file}')
                os.system(cmd_mv)