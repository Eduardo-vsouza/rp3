import os

from ..pipeline_config import PipelineStructure


class BaseSearch(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
    
        self.args = args
        self.print_row(word="Search", color='blue')
        self.print_state(message=f"Engine: {self.args.engine}", color='blue', marker='')
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
            self.check_dirs([self.cascadeDir, self.cascadeMzmlDir, self.cascadeFirstPassDir, self.cascadeSecondPassDir])

    def check_files(self, files):
        checked = ''
        filelist = files.split(" ")
        for file in filelist:
            if file.endswith(self.args.fileFormat):
                checked += f' {file}'
        return checked

    def move_pin_files(self, mzml_dir=None, outdir="standard", split_i=None):
        # outdirs = {"standard": }
        if outdir == 'standard':
            db_relative = self.select_database(decoy=True).split("/")[-1]
            output_dir = f'{self.searchDir}/group/{db_relative}'
        elif outdir == 'rescore':
            output_dir = f'{self.outdir}/rescore/peptide_search/group'
        # elif outdir == 'cascade_first':
        #     output_dir = self.cascadeFirstPassDir
        # elif outdir == 'cascade_second':
        #     output_dir = self.cascadeSecondPassDir
        else:
            output_dir = outdir
        if mzml_dir is None:
            mzml_dir = self.args.mzml
        files = os.listdir(mzml_dir)
        for i, file in enumerate(files):
            pattern = '_target.pin'
            if split_i is not None:
                pattern = f'_{split_i}_target.pin'
            if file.endswith(".pin"):
                cmd_mv = (f'mv {mzml_dir}/{file} '
                          f'{output_dir}/{file.replace(".pin", pattern)}')
                os.system(cmd_mv)
            elif file.endswith(".txt") or file.endswith(".xml"):
                cmd_mv = (f'mv {mzml_dir}/{file} {output_dir}/{file}')
                os.system(cmd_mv)

    
    def get_mzml(self, mzml_dir):
        files = ''
        mzml = os.listdir(mzml_dir)
        files = ''
        for file in mzml:
            if file.endswith(self.args.fileFormat):
                files += f' {mzml_dir}/{file}'
        return files