import os
import sys

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from ..pipeline_config import PipelineStructure
from ..utils import ProtSplit


class MS2Summ(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        print(f"--Summarizing results...")
        # self.proteinGroups = self.__get_protein_groups()
        self.expGroups = {}
        if self.args.groupsFile is not None:
            self.expGroups = self.associate_groups_to_files()   

        self.masterDF = {'Protein': [], 'Sequence': [], 'Protein Group': [], 'Gene Name': [], 'Total Peptides': [], 'Total spectral counts': []}
        self.masterDict = {}

    def __get_protein_groups(self):
        print(f"--Classifying protein into groups...")
        prot_groups = {}
        protsplit = ProtSplit(outdir=self.args.outdir)
        protsplit.split_protein_groups()
        this_prot_groups = protsplit.get_protein_groups(order='prot_group')  # dict: protein, group
        prot_groups.update(this_prot_groups)

        return prot_groups
    
    def __retrieve_gene_name(self, protein):
        if 'GN=' in protein:
            gene_name = protein.split("GN=")[1].split("_")[0]
            return gene_name
        else:
            return ''

    def gather_peptides(self):
        print(f"--Gathering peptides...")
        filepath = self.select_peptides_df()
        peptides_df = pd.read_csv(filepath, sep='\t')
        peptides_df = peptides_df[peptides_df["q-value"] != 'q-value']
        peptides_df["q-value"] = peptides_df["q-value"].astype(float)
        peptides_df = peptides_df[peptides_df["q-value"] <= 0.01]
        print(f"Total peptides at 1% FDR: {len(peptides_df)}")

        peptides_df = peptides_df[~peptides_df["proteinIds"].str.contains("rev_")]
        peptides_df = peptides_df[~peptides_df["proteinIds"].str.contains("contaminant")]
        psm_dict = peptides_df.set_index("PSMId").to_dict(orient="index")
        for psm in tqdm(psm_dict):
            prot_list = psm_dict[psm]["proteinIds"]
            selected, unique = self.__judge_protein_list(prot_list)
            # print(selected)
            if len(selected) == 0:
                continue
            else:
                for protein in selected:
                    file = psm.split(".")[0]
                    if self.args.groupsFile is not None:
                        exp_group = self.expGroups[file]
                    else:
                        exp_group = 'group'
                    if protein not in self.masterDict:
                        self.masterDict[protein] = {}
                        self.masterDict[protein]['Gene Name'] = self.__retrieve_gene_name(protein)
                        # self.masterDict[protein]['Protein Group'] = self.proteinGroups[protein]

                        self.masterDict[protein]['Total Peptides'] = 0
                        self.masterDict[protein]['Total Spectral counts'] = 0
                        self.masterDict[protein]['Total Unique spectral counts'] = 0
                        self.masterDict[protein]['Total Unique Peptides'] = 0

                    if f'Peptides Seq {exp_group}' not in self.masterDict[protein]:
                        self.masterDict[protein][f'Peptides Seq {exp_group}'] = ''
                        self.masterDict[protein][f'Peptides Count {exp_group}'] = 0

                        self.masterDict[protein][f'Unique Peptides Seq {exp_group}'] = ''
                        self.masterDict[protein][f'Unique Peptides Count {exp_group}'] = 0

                        self.masterDict[protein][f'Spectral Counts {exp_group}'] = 0
                        self.masterDict[protein][f'Unique Spectral Counts {exp_group}'] = 0


                    self.masterDict[protein]['Total Peptides'] += 1
                    self.masterDict[protein][f'Peptides Seq {exp_group}'] += f'{psm_dict[psm]["peptide"]},'
                    self.masterDict[protein][f'Peptides Count {exp_group}'] += 1
                    if unique:
                        self.masterDict[protein][f'Unique Peptides Seq {exp_group}'] += f'{psm_dict[psm]["peptide"]},'
                        self.masterDict[protein][f'Unique Peptides Count {exp_group}'] += 1
                        self.masterDict[protein]['Total Unique Peptides'] += 1
        # print(self.masterDict)

    def gather_spectral_counts(self):
        print(f"--Gathering spectral counts...")
        filepath = self.select_psm_df()
        psm_df = pd.read_csv(filepath, sep='\t')
        peptides_df = psm_df[psm_df["q-value"] != 'q-value']
        peptides_df["q-value"] = peptides_df["q-value"].astype(float)
        peptides_df = peptides_df[peptides_df["q-value"] <= 0.01]
        print(f"Total PSMs at 1% FDR: {len(peptides_df)}")
        peptides_df = peptides_df[~peptides_df["proteinIds"].str.contains("rev_")]
        peptides_df = peptides_df[~peptides_df["proteinIds"].str.contains("contaminant")]
        psm_dict = peptides_df.set_index("PSMId").to_dict(orient="index")
        for psm in tqdm(psm_dict):
            prot_list = psm_dict[psm]["proteinIds"]
            selected, unique = self.__judge_protein_list(prot_list)
            if len(selected) == 0:
                continue
            else:
                for protein in selected:            # since we define proteins in thed ict before with gather_peptides(),
                    if protein in self.masterDict:  # use only the PSMs with at least one accepted peptide 
                        file = psm.split(".")[0]    
                        if self.args.groupsFile is not None:
                            exp_group = self.expGroups[file]
                        else:
                            exp_group = 'group'
                        self.masterDict[protein]['Total Spectral counts'] += 1


                        if f'Spectral Counts {exp_group}' not in self.masterDict[protein]:
                            self.masterDict[protein][f'Spectral Counts {exp_group}'] = 0
                            self.masterDict[protein][f'Unique Spectral Counts {exp_group}'] = 0


                        self.masterDict[protein][f'Spectral Counts {exp_group}'] += 1
                        if unique:
                            self.masterDict[protein][f'Unique Spectral Counts {exp_group}'] += 1
                            self.masterDict[protein]['Total Unique spectral counts'] += 1

    def __judge_protein_list(self, prot_list):
        """
        Accepts a comma-separated list of protein matches for a given MS peptide or PSM.
        Selects based on criteria which proteins to accept
        """
        proteins = prot_list.split(",")
        # print(proteins)
        unique = True
        # if only smORFs are in prot_list, return a list of all proteins
        
        if 'ANNO' not in prot_list and '_F:' in prot_list:
        
            # if more than one protein is in the list, classify as non-unique
            if ',' in prot_list:
                unique = False
            return proteins, unique            

        selected = []
        # if any annotated protein is in the list, do not return any novel
        all_novel = True
        annotated = 0

        for protein in proteins:
            if 'ANNO' in protein or 'GN=' in protein:
                all_novel = False
                annotated += 1
                selected.append(protein)
            else:
                if '_F:' in protein:
                    if 'ANNO' not in prot_list:
                        selected.append(protein)  # this makes sure that we do not include a novel one whose peptide also matches an annotated
        # 3) razor assignment
        if self.args.pepAssign == 'razor' and len(selected) > 1:
            # a) if any annotated, restrict to those
            anno_hits = [p for p in selected if 'ANNO' in p or 'GN=' in p]
            candidates = anno_hits if anno_hits else selected

            # b) pick the one with highest Total Peptides so far
            counts = {
                p: self.masterDict.get(p, {}).get('Total Peptides', 0)
                for p in candidates
            }
            razor_prot = max(counts, key=counts.get)
            return [razor_prot], True

        # 4) fallback to unique/non-unique logic
        if annotated <= 1:
            return selected, unique
        else:
            return selected, False
        
    def gather_protein_sequences(self):
        print(f"--Gathering protein sequences...")
        db = self.select_database()
        records = SeqIO.parse(db, 'fasta')
        db_seqs = {}
        print(f"--Loading database sequences...")
        for record in tqdm(records):
            entry = str(record.description)
            seq = str(record.seq)
            db_seqs[entry] = seq
        if self.args.cascade:
            db = self.select_database(proteome=True, decoy=False)
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                entry = str(record.description)
                seq = str(record.seq)
                if entry not in db_seqs:
                    db_seqs[entry] = seq
        print(f"Total sequences in database: {len(db_seqs)}")
        added = 0
        print("--Inserting protein sequences into summary...")
        for protein in tqdm(self.masterDict):
            if protein in db_seqs:
                self.masterDict[protein]['Sequence'] = db_seqs[protein]
                added += 1
            else:
                self.masterDict[protein]['Sequence'] = ''
                print(f"Warning: {protein} not found in database. This suggest an error when processing the pipeline.")
        print(f"--Total sequences added: {added}")

    def save(self):
        df = pd.DataFrame.from_dict(self.masterDict, orient='index').reset_index().rename(columns={'index': 'Protein'})
        rescored = self.is_rescored()
        if rescored:
            outdir = self.rescoreSummarizedResultsDir
        else:
            outdir = self.summarizedResultsDir
        # Reorder dataframe columns so that the experimental group columns are grouped together
        base_cols = [col for col in df.columns if col in ['Protein', 'Sequence', 'Gene Name',
                                                           'Total Peptides', 'Total Unique Peptides',
                                                           'Total Spectral counts', 'Total Unique spectral counts']]
        peptides_seq = sorted([col for col in df.columns if col.startswith('Peptides Seq ')])
        peptides_count = sorted([col for col in df.columns if col.startswith('Peptides Count ')])
        unique_peptides_seq = sorted([col for col in df.columns if col.startswith('Unique Peptides Seq ')])
        unique_peptides_count = sorted([col for col in df.columns if col.startswith('Unique Peptides Count ')])
        spectral_counts = sorted([col for col in df.columns if col.startswith('Spectral Counts ') and col != 'Total Spectral counts'])
        unique_spectral_counts = sorted([col for col in df.columns if col.startswith('Unique Spectral Counts ') and col != 'Total Unique spectral counts'])

        new_order = (base_cols + peptides_seq + peptides_count +
                     unique_peptides_seq + unique_peptides_count +
                     spectral_counts + unique_spectral_counts)
        df = df[new_order]
        df = df[df["Total Unique Peptides"] > 0]
        df.to_csv(f'{outdir}/ms2summ.csv', sep='\t', index=False)
        print(f"MS2 summary file saved to {outdir}/ms2summ.csv")
        total_unique_sc = df[df["Total Unique spectral counts"] > 0]
        total_unique_pep = df[df["Total Unique Peptides"] > 0]
        print(f"Total proteins with unique spectral counts: {len(total_unique_sc)}")
        print(f"Total proteins with unique peptides: {len(total_unique_pep)}")
        prots = len(df["Protein"].tolist())
        print(f"Total proteins: {prots}")

            

