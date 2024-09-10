import os

from Bio import SeqIO
import pandas as pd

from ..pipeline_config import PipelineStructure


class ProtSplit:
    def __init__(self, outdir):
        self.outdir = outdir

        self.db = self.__select_db()

        self.proteinGroupsDir = f'{outdir}/protein_groups'
        PipelineStructure.check_dirs([self.proteinGroupsDir])
        self.sequences = {}

        self.proteinDict = self.__split_results()
        self.proteinGroups = f'{self.proteinGroupsDir}/protein_groups.csv'

    def __select_db(self):
        db = None
        if os.path.exists(f'{self.outdir}/rescore/summarized_results/filtered_rescored_microproteins_150.fasta'):
            db = f'{self.outdir}/rescore/databases/rescore_target_database.fasta'
        else:
            dbs = os.listdir(f'{self.outdir}/databases')
            for fasta in dbs:
                if fasta.endswith("target_database.fasta"):
                    db = f'{self.outdir}/databases/{fasta}'
        return db

    def split_protein_groups(self):
        print(f"--Splitting proteins into standard-sized proteins, unannotated microproteins, and "
              f"annotated microproteins.")
        data = {'protein': [], 'group': []}
        for group in self.proteinDict:
            for prot in self.proteinDict[group]:
                data['protein'].append(prot)
                data['group'].append(group)
        df = pd.DataFrame(data)
        df.to_csv(self.proteinGroups, sep='\t', index=False)
        # self.proteinDict = {}

    def __split_results(self):
        file = f'{self.outdir}/rescore/post_processing/group/peptides_fixed.txt'
        if not os.path.exists(file):
            file = f'{self.outdir}/post_processing/group/db/peptides_fixed.txt'
        df = pd.read_csv(file, sep='\t')
        df = df[df["q-value"] <= 0.01]
        df = df[~df["proteinIds"].str.contains("rev_")]
        df = df[~df["proteinIds"].str.contains("contaminant")]
        prots = df["proteinIds"].tolist()
        filtered_proteins = []
        for prot in prots:
            proto = prot.replace(",_", "_").replace("__", "_")
            if ',' not in prots:
                proteins = proto.split(",")
                for protein in proteins:
                    if 'ANNO' not in prots:
                        filtered_proteins.append(protein)
        protein_dict = self.__get_prot_dict(filtered_proteins)
        return protein_dict

    def __get_prot_dict(self, filtered_proteins, mp_threshold=150):
        protein_dict = {'standard': [], 'annotated_microproteins': [], 'microproteins': []}

        records = SeqIO.parse(self.db, 'fasta')
        for record in records:
            entry = str(record.description)
            seq = str(record.seq)
            # if 'sp|' in entry:
            entry = entry.replace(",_", "_").replace("__", "_")
            self.sequences[entry] = seq
            if entry in filtered_proteins:
                if len(seq) <= mp_threshold:

                        if '_ANNO' not in entry:
                                protein_dict['microproteins'].append(entry)
                        else:
                            protein_dict['annotated_microproteins'].append(entry)
                else:
                    protein_dict['standard'].append(entry)
        return protein_dict

    def get_protein_groups(self, order='group_prot'):
        """
        param order: group_prot returns a dictionary with group, prot as key, value. prot_group returns a
        dictionary with prot, group as key, value
        return: a dictionary containing key, value as protein, group or group, protein, as speficied by order
        """
        protein_groups = {}
        df = pd.read_csv(self.proteinGroups, sep='\t')
        proteins, groups = df["protein"].tolist(), df["group"].tolist()
        for prot, group in zip(proteins, groups):
            if order == 'group_prot':
                if group not in protein_groups:
                    protein_groups[group] = []
                protein_groups[group].append(prot)
            else:
                protein_groups[prot] = group
        return protein_groups

    def split_fasta(self):
        groups = self.get_protein_groups()
        for group in groups:
            fasta = []
            for entry in self.sequences:
                if entry in groups[group]:
                    fasta.append(f'>{entry}\n{self.sequences[entry]}\n')
            with open(f'{self.proteinGroupsDir}/{group}.fasta', 'w') as out:
                out.writelines(fasta)
