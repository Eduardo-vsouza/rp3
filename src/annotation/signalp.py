import os
import sys

from ..pipeline_config import PipelineStructure


class SignalP(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.plotsDir = f'{self.signalPDir}/plots'
        self.check_dirs([self.signalPDir, self.plotsDir])

    def run(self):
        files = os.listdir(self.mergedResults)
        # print(files)
        # if self.microproteinsBlast.split("/")[-1] in files:
        file = self.microproteinsBlast
        if self.args.rescored:
            file = self.rescoredMicroproteinsFasta
        # else:
        #     file = self.uniqueMicroproteins
        # print(file)
        cmd = f'{self.toolPaths["signalP"]} --fastafile {file} --format all ' \
              f'--organism {self.args.organism} --output_dir {self.signalPDir} --mode slow-sequential'
        os.system(cmd)

    def organize_files(self):
        # self.check_dirs(f'{plots}')
        cmd = f'mv {self.signalPDir}/*plot.* {self.plotsDir}/.'
        os.system(cmd)

    def save_metrics(self):
        cmd = f'grep ">" {self.signalPDir}/processed_entries.fasta | wc -l > {self.metricsDir}/signal_peptides.txt'
        os.system(cmd)