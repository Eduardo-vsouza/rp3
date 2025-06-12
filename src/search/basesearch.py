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
        self.__check_dirs()
        self.threads = self.args.threads
        self.mode = 'search'
        self.params = []
