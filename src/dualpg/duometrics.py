import os
import sys

from Bio import SeqIO

from ..pipeline_config import PipelineStructure


class DuoMetrics(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

    def separate_results(self):
        if os.path.exists(self.rescoreDir):
            self.__divide_fasta(pattern=self.args.hostPattern, fasta=self.rescoredMicroproteinsFasta,
                                outfile=self.rescoredHostMicroproteins)
            self.__divide_fasta(pattern=self.args.pathoPattern, fasta=self.rescoredMicroproteinsFasta,
                                outfile=self.rescoredPathogenMicroproteins)
        self.__divide_fasta(pattern=self.args.hostPattern, fasta=self.uniqueMicroproteins,
                            outfile=self.hostMicroproteins)
        self.__divide_fasta(pattern=self.args.pathoPattern, fasta=self.uniqueMicroproteins,
                            outfile=self.pathogenMicroproteins)

    @staticmethod
    def __divide_fasta(pattern, fasta, outfile):
        records = SeqIO.parse(fasta, 'fasta')
        to_write = []
        for record in records:
            entry = str(record.description)
            if pattern in entry:
                to_write.append(f'>{entry}\n{str(record.seq)}\n')
        with open(outfile, 'w') as handler:
            handler.writelines(to_write)

    def count_orfs(self):
        files = os.listdir(self.duoDir)
        for file in files:
            if file.endswith("fasta"):
                cmd = f'echo "{file}\n" >> {self.duoDir}/microproteins_report.txt'
                os.system(cmd)
                cmd_report = f'mip_report.sh {self.duoDir}/{file} >> {self.duoDir}/microproteins_report.txt'
                os.system(cmd_report)
                cmd = f'echo "---------------\n" >> {self.duoDir}/microproteins_report.txt'
                os.system(cmd)

