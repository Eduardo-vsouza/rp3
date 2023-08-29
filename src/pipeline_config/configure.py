import os

from .pipeline_structure import PipelineStructure


class Configure(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.toolPaths = {}

    def read_config_file(self):
        with open(self.configFile, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                splat = line.split("\t")
                tool = splat[0]
                path = splat[1]
                self.toolPaths[tool] = path
        return self.toolPaths

