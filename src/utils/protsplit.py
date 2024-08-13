from Bio import SeqIO
import pandas as pd

from ..pipeline_config import PipelineStructure


class ProtSplit(PipelineStructure):
    def __init__(self, args, outdir):
        super().__init__(args=args)
        self.outdir = outdir
        self.check_dirs([self.proteinGroupsDir])
        # self.peptidesDF
        self.proteinDict = self.__split_results()

    def split_protein_groups(self):
        run = self.verify_checkpoint(outfile=self.proteinGroups, step="protein group splitting")
        if run:
            data = {'protein': [], 'group': []}
            for group in self.proteinDict:
                for prot in self.proteinDict[group]:
                    data['protein'].append(prot)
                    data['group'].append(group)
            df = pd.DataFrame(data)
            df.to_csv(self.proteinGroups, sep='\t', index=False)
            self.proteinDict = {}

    def __split_results(self):
        run = self.verify_checkpoint(outfile=self.proteinGroups, step="protein group splitting", mute=True)
        if run:
            file = self.select_peptides_df()
            file = f'{self.rescorePostProcessDir}/group/peptides_fixed.txt'
            df = pd.read_csv(file, sep='\t')
            df = df[df["q-value"] <= 0.01]
            df = df[~df["proteinIds"].str.contains("rev_")]
            df = df[~df["proteinIds"].str.contains("contaminant")]
            prots = df["proteinIds"].tolist()
            filtered_proteins = []
            for prot in prots:
                proteins = prot.split(",")
                for protein in proteins:
                    filtered_proteins.append(protein)
            protein_dict = self.__get_prot_dict(filtered_proteins)
            return protein_dict
        else:
            return None

    def __get_prot_dict(self, filtered_proteins, mp_threshold=150):
        protein_dict = {'standard': [], 'annotated_microproteins': [], 'microproteins': []}
        db = self.select_database(decoy=False)
        records = SeqIO.parse(db, 'fasta')
        for record in records:
            entry = str(record.description)
            if entry in filtered_proteins:
                seq = str(record.seq)
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
