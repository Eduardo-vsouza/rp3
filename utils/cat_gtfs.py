import os
import sys

from ..pipeline_config import PipelineStructure



class GTFGatherer(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)


    def cat_gtfs(self):
        dbs = os.listdir(self.translationDir)
        gtf = []
        filesa = ''
        for db in dbs:
            fulldir = f'{self.translationDir}/{db}'
            files = os.listdir(fulldir)
            for file in files:
                if file.endswith("_ORFs.gtf"):
                    # with open(f'{fulldir}/{file}', 'r') as handler:
                    #     lines = handler.readlines()
                    #     for line in lines:
                    #         if line not in gtf:
                    #             gtf.append(line)
                    filesa += f' {fulldir}/{file}'
        # with open(self.mergedFullGTF, 'w') as outfile:
        #     outfile.writelines(gtf)
        cmd = f'cat {filesa} > {self.mergedFullGTF}'
        print(cmd)
        os.system(cmd)