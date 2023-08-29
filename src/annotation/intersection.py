import os
import sys

from Bio import SeqIO
from venn import venn
import matplotlib.pyplot as plt

from ..pipeline_config import PipelineStructure


class ResultsIntersection(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.vennSets = {}

    def compare_groups(self):
        groups = os.listdir(self.summarizedResultsDir)
        for group in groups:
            name = '_'.join(group.split("_")[:2])
            if group != 'merged' and name not in self.vennSets:
                self.vennSets[name] = []
            groupdir = f'{self.summarizedResultsDir}/{group}'
            records = SeqIO.parse(f'{groupdir}/microproteins_150.fasta_blast_filt.fasta', 'fasta')
            for record in records:
                seq = str(record.seq)
                if name in self.vennSets:
                    self.vennSets[name].append(seq)
        for subset in self.vennSets:
            self.vennSets[subset] = set(self.vennSets[subset])
        fig = plt.figure()
        venn(self.vennSets)
        plt.savefig(f'{self.metricsDir}/microprotein_blast_filt_venn.png', dpi=600)
