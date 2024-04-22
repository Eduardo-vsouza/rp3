import os
import sys

from Bio import SeqIO
import pandas as pd

from ..pipeline_config import PipelineStructure


class ResultsSummary(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.data = {}

    def gather_smorf_data(self):
        records = SeqIO.parse(self.microproteinsBlast, 'fasta')
        for record in records:
            self.data[str(record.description)] = {'seq': str(record.seq), "signalP": 0, "conservation": [],
                                                  "rescored": "False", "Ribo-Seq coverage": "False"}
        records = SeqIO.parse(self.rescoredMicroproteinsFasta, 'fasta')
        for record in records:
            self.data[str(record.description)]["rescored"] = "True"

    def get_signalp(self):
        records = SeqIO.parse(f'{self.signalPDir}/processed_entries.fasta', 'fasta')
        for record in records:
            self.data[str(record.description)]["signalP"] = str(record.seq)

    def get_conservation(self):
        df = pd.read_csv(f'{self.phyloDir}/smorfs_entries_per_species.xls', sep='\t')
        species, smorfs = df["species"].tolist(), df["smorf"].tolist()
        for sp, smorf in zip(species, smorfs):
            if sp not in self.data[smorf]["conservation"]:
                self.data[smorf]["conservation"].append(sp)

    def get_riboseq_cov(self):
        df = pd.read_csv(self.mappingGroups, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            if group == "No coverage":
                self.data[smorf]["Ribo-Seq coverage"] = "False"
            else:
                self.data[smorf]["Ribo-Seq coverage"] = "True"

    def save(self):
        data = {"smorf": [], "signalP": [], "conservation": [], "rescored": [], "Ribo-Seq coverage": []}
        for smorf in self.data:
            data["smorf"].append(smorf)
            data["signalP"].append(self.data[smorf]["signalP"])
            data["conservation"].append(','.join(self.data[smorf]["conservation"]))
            data["rescored"].append(self.data[smorf]["rescored"])
            data["Ribo-Seq coverage"].append(self.data[smorf]["Ribo-Seq coverage"])
        df = pd.DataFrame(data=data)
        df.to_csv(f'{self.metricsDir}/smorfs_summary.xls', sep='\t', index=False)

