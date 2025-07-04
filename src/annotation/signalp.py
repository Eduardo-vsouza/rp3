import os
import sys

from ..pipeline_config import PipelineStructure


class SignalP(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.plotsDir = f'{self.signalPDir}/plots'
        self.check_dirs([self.signalPDir, self.plotsDir])

    def run(self):
        self.check_dirs([self.signalPDir, self.signalPstandardDir,
                         self.signalPAnnoMPDir, self.signalPMicroproteinDir])

        # rp3 microproteins
        # if self.args.rescored:
        #     file = self.rescoredMicroproteinsFasta
        # else:
        #     file = self.uniqueMicroproteins
        file = self.select_fasta()
        self.__run_signalp(file, self.signalPMicroproteinDir)

        if self.args.predictAnnotated:
            # annotated microproteins
            file = f'{self.proteinGroupsDir}/annotated_microproteins.fasta'
            if os.path.exists(file):
                self.__run_signalp(file, self.signalPAnnoMPDir)

            # standard-sized proteins
            file = f'{self.proteinGroupsDir}/standard.fasta'
            if os.path.exists(file):
                self.__run_signalp(file, self.signalPstandardDir)

    def __run_signalp(self, file, outdir):
        mode = 'fast'
        if self.args.signalpMode == 'slow':
            mode = 'slow-sequential'
        cmd = f'{self.toolPaths["signalP"]} --fastafile {file} --format all --organism eukarya ' \
              f'--organism {self.args.organism} --output_dir {outdir} --mode {mode} -wp ' \
              f'{self.args.threads} -tt {self.args.threads}'
        os.system(cmd)

    def organize_files(self):
        # self.check_dirs(f'{plots}')
        cmd = f'mv {self.signalPDir}/*plot.* {self.plotsDir}/.'
        os.system(cmd)

    def save_metrics(self):
        cmd = f'grep ">" {self.signalPDir}/processed_entries.fasta | wc -l > {self.metricsDir}/signal_peptides.txt'
        os.system(cmd)