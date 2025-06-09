import os

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure
from ..pgcontext import PGContext
from ..utils import ProtSplit



class SpecComparison:
    def __init__(self, args):
        self.args = args

        self.filteredMicroproteins = {}
        self.__read_fasta_results()

        self.proteins = {}
        self.peptides = {}
        self.overlaps = {}
        self.proteinSequences = {}
        self.proteinsUnique = {}

        self.mergedProteins = {'protein': [], 'sequence': [], 'overlaps': [], 'protein_group': []}

        self.proteinGroups = self.__get_protein_groups()



    def __get_protein_groups(self):
        prot_groups = {}

        for results in self.args.results:
            protsplit = ProtSplit(outdir=results)
            protsplit.split_protein_groups()
            this_prot_groups = protsplit.get_protein_groups(order='prot_group')  # dict: protein, group
            prot_groups.update(this_prot_groups)
            self.proteinSequences.update(protsplit.sequences)

        return prot_groups

    def get_spec_counts(self):
        outdirs, groups = self.args.results, self.args.groups
        for outdir, group in zip(outdirs, groups):
            # fasta = self.select_fasta()
            # rescored = self.check_rescore()
            if self.args.rescored:
                fasta = f'{outdir}/rescore/summarized_results/filtered_rescored_microproteins_150.fasta'
                pep_df = f'{outdir}/rescore/post_processing/group/peptides_fixed.txt'
            else:
                pep_df = f'{outdir}/post_processing/group/db/peptides_fixed.txt'
            df = pd.read_csv(pep_df, sep='\t')
            df = df[df["q-value"] != "q-value"]
            df = df[df["q-value"] <= 0.01]
            df = df[~df["proteinIds"].str.contains("rev_")]
            df = df[~df["proteinIds"].str.contains("contaminant")]
            proteins = df["proteinIds"].tolist()
            peptides = df["peptide"].tolist()
            # df = df[~df["proteinIds"].str]
            for protein, peptide in zip(proteins, peptides):
                proteina = protein.replace(",_", "_").replace("__", "_")
                protlist = proteina.split(",")
                for prot in protlist:
                    # if 'sp|' in prot:
                    #     prot = '|'.join(protein.split("|")[:2])
                    # if '_F:' in prot and 'ANNO' not in protein:
                    #     add = True
                    #     print(prot)
                    # else:
                    #     add = False
                    # print(prot)
                    add = self.__check_filters(prot, proteina)
                    # print(add)
                    if add:
                        if prot not in self.proteins:
                            self.proteins[prot] = {}
                            self.peptides[prot] = {}
                        if group not in self.proteins[prot]:
                            self.proteins[prot][group] = 0
                            self.peptides[prot][group] = []
                        self.proteins[prot][group] += 1
                        if peptide not in self.peptides[prot][group]:
                            self.peptides[prot][group].append(peptide)
                    if ',' not in proteina:  # check if they are unique
                        if prot not in self.proteinsUnique:
                            self.proteinsUnique[prot] = {}
                        if group not in self.proteinsUnique[prot]:
                            self.proteinsUnique[prot][group] = 0
                        self.proteinsUnique[prot][group] += 1


        for group in self.args.groups:
            for prot in self.proteins:
                if group not in self.proteins[prot]:
                    self.proteins[prot][group] = 0
                    self.peptides[prot][group] = ''
            for prot in self.proteinsUnique:
                if group not in self.proteinsUnique[prot]:
                    self.proteinsUnique[prot][group] = 0

    def add_overlaps(self):
        for outdir in self.args.results:
            pgc = PGContext(args=self.args)
            pgc.outdir = outdir
            pgc.overlappedGTF = f'{outdir}/pg_context/intermediate_files/overlapped_smorfs.gtf'
            if os.path.exists(pgc.overlappedGTF):
                pgc.gather_overlaps()
                for smorf in pgc.overlappedSmorfs:
                    if smorf not in self.overlaps:
                        self.overlaps[smorf] = []
                    for feature in pgc.overlappedSmorfs[smorf]:
                        orf = feature.gene
                        # print(orf)
                        if orf is not None:
                            self.overlaps[smorf].append(orf)

    def create_data_frame(self):
        for group in self.args.groups:
            if group not in self.mergedProteins:
                self.mergedProteins[group] = []
        for group in self.args.groups:
            self.mergedProteins[f'{group}_peptides'] = []
            self.mergedProteins[f'{group}_uniquePeptides'] = []
        for protein in self.proteins:
            # if 'sp|' in protein:
            #     protein = '|'.join(protein.split("|")[:2])
            # if protein in self.proteinGroups:

            self.mergedProteins['protein'].append(protein)
            if protein in self.proteinSequences:
                seq = self.proteinSequences[protein].rstrip()
            else:
                seq = ''
            if protein in self.proteinGroups:
                protgroup = self.proteinGroups[protein]
            else:
                protgroup = ''
            self.mergedProteins['protein_group'].append(protgroup)

            self.mergedProteins['sequence'].append(seq)
            if protein in self.overlaps:
                # if self.overlaps[protein] is not None:
                self.mergedProteins['overlaps'].append(', '.join(list(set(self.overlaps[protein]))))
            else:
                self.mergedProteins['overlaps'].append('Intergenic')
            for group in self.args.groups:
                self.mergedProteins[f'{group}_peptides'].append(','.join(self.peptides[protein][group]))
                if protein in self.proteinsUnique:
                    self.mergedProteins[f'{group}_uniquePeptides'].append(self.proteinsUnique[protein][group])
                else:
                    self.mergedProteins[f'{group}_uniquePeptides'].append(0)
            # if len(self.proteinSequences[protein]) > 150:
            #     print(protein, "longer than 150")
            for group in self.proteins[protein]:
                self.mergedProteins[group].append(self.proteins[protein][group])
            # else:
            #     print(protein)
        PipelineStructure.check_dirs([self.args.outdir])
        df = pd.DataFrame(data=self.mergedProteins)
        # df.to_excel(f'{self.args.outdir}/group_comparison.xlsx', index=False, encoding='utf-8')

        df.to_csv(f'{self.args.outdir}/group_comparison.xlsx', sep='\t', index=False)
        df.to_csv(f'{self.args.outdir}/group_comparison.csv', sep='\t', index=False)
        df.to_csv(f'{self.args.outdir}/group_comparison.tsv', sep='\t', index=False)
        # filter it a little
        # df = df[df[""]]
        groups = []
        for group in self.args.groups:
            groups.append(group)
        df = df[(df[f"{groups[0]}_uniquePeptides"] > 0) | (df[f"{groups[1]}_uniquePeptides"] > 0)]
        df.to_csv(f'{self.args.outdir}/group_comparison_filtered.csv', sep='\t', index=False)
        full_fasta = []
        mp_fasta = []
        anno_mp_fasta = []


        entries, seqs = df["protein"].tolist(), df["sequence"].tolist()
        for entry, seq in zip(entries, seqs):
            full_fasta.append(f'>{entry}\n{seq}\n')
        with open(f'{self.args.outdir}/group_comparison.fasta', 'w') as handler:
            handler.writelines(full_fasta)

        for entry, seq in zip(entries, seqs):
            if '_F:' in entry:
                mp_fasta.append(f'>{entry}\n{seq}\n')
            if 'ANNO' in entry:
                if len(seq) <= 150:
                    anno_mp_fasta.append(f'>{entry}\n{seq}\n')
        with open(f'{self.args.outdir}/group_comparison_mp.fasta', 'w') as handler:
            handler.writelines(mp_fasta)
        with open(f'{self.args.outdir}/group_comparison_anno_mp.fasta', 'w') as handler:
            handler.writelines(anno_mp_fasta)
        
    def __read_fasta_results(self):


        for outdir, group in zip(self.args.results, self.args.groups):
            if self.args.rescored:
                fasta = f'{outdir}/rescore/summarized_results/filtered_rescored_microproteins_150.fasta'
            else:
                fasta = f'{outdir}/summarized_results/merged/microproteins_150.fasta'
            records = SeqIO.parse(fasta, 'fasta')
            for record in records:
                entry = str(record.description)

                if not self.args.microproteins:
                    self.filteredMicroproteins[entry] = str(record.seq)
                else:
                    if len(str(record.seq)) <= 150:
                        self.filteredMicroproteins[entry] = str(record.seq)
                # if 'ALFA' in entry:
                #     self.filteredMicroproteins[entry] = str(record.seq)
        # self.filteredMicroproteins['E4orf6/7-ALFA'] = 'MAANSTSDLHTPGTQLSVADIIVITVYFALNVAVGIWSSCRASRNTVNGYFLAGRDMTWWPIGASLFASSEGSGLFIGLAGSGAAGGLAVAGFEWNATYVLLALAWVFVPIYISSEIVTLPEYIQKRYGGQRIRMYLSVLSLLLSVFTKISLDLYAGALFVHICLGWNFYLSTILTLGITALYTIAGGLAAVIYTDALQTLIMVVGAVILTIKAFDQIGGYGQLEAAYAQAIPSRTIANTTCHLPRTDAMHMFRDPHTGDLPWTGMTFGLTIMATWYWCTDQVIVQRSLSARDLNHAKAGSILASYLKMLPMGLIIMPGMISRALFPDDVGCVVPSECLRACGAEVGCSNIAYPKLVMELMPIGLRGLMIAVMLAALMSSLTSIFNSSSTLFTMDIWRRLRPRSGERELLLVGRLVIVALIGVSVAWIPVLQDSNSGQLFIYMQSVTSSLAPPVTAVFVLGVFWRRANEQGAFWGLIAGLVVGATRLVLEFLNPAPPCGEPDTRPAVLGSIHYLHFAVALFALSGAVVVAGSLLTPPPQSVQIENLTWWTLAQDVPLGTKAGDGQTPQKHAFWARVCGFNAILLMCVNIFFYAYFA'

    def separate_up_regulated(self):
        df = pd.read_csv(f'{self.args.outdir}/group_comparison.xlsx', sep='\t')
        # df = pd.read_csv(f'{self.args.outdir}/group_comparison.csv', sep='\t')

        # for group in self.args.groups:
        #     if group != self.highlightGroup:

        # filtered_df = df[(df["1640"] > df['UI']) & (df['1640'] > df['102'])]
        filtered_df = df[(df[self.args.groups[-1]] > df[self.args.groups[0]])]
        treat = self.args.groups[-1]
        filtered_df.to_csv(f'{self.args.outdir}/{treat}_upregulated.xlsx', sep='\t', index=False)
        filtered_df.to_csv(f'{self.args.outdir}/{treat}_upregulated.csv', sep='\t', index=False)
        filtered_df.to_csv(f'{self.args.outdir}/{treat}_upregulated.tsv', sep='\t', index=False)


        up_control = df
        for group in self.args.groups[:-1]:

            up_control = up_control[(up_control[group] < up_control[self.args.groups[-1]])]
        up_control.to_csv(f'{self.args.outdir}/{treat}_up_in_{self.args.groups[-1]}.csv', sep='\t', index=False)


        # unique = df[(df['UI'] == 0) & (df['102'] == 0) & (df['1640'] > 0)]
        # unique = df[(df[self.args.groups[0]] == 0) & (df[self.args.groups[-1]] > 0)]
        unique = df
        for group in self.args.groups[:-1]:
            print(group)
            print(self.args.groups[-1])
            unique = unique[(unique[group] == 0) & (unique[self.args.groups[-1]] > 0)]

        unique.to_csv(f'{self.args.outdir}/{treat}_unique.xlsx', sep='\t', index=False)
        unique.to_csv(f'{self.args.outdir}/{treat}_unique.csv', sep='\t', index=False)
        unique.to_csv(f'{self.args.outdir}/{treat}_unique.tsv', sep='\t', index=False)

        # unique_control = df[(df[self.args.groups[-1]] == 0) & (df[self.args.groups[0]] > 0)]
        # unique_control.to_csv(f'{self.args.outdir}/{self.args.groups[0]}_unique.csv', sep='\t', index=False)




    def __check_filters(self, prot, proteins):
        add = True
        if self.args.predictedOnly:
            if '_F:' not in prot:
                add = False
        # if prot not in self.filteredMicroproteins:
        #     add = False
        if 'OS=' in proteins and self.args.predictedOnly:
            add = False
        if 'ALFA' in prot and 'ANNO' not in prot:
            add = True
        if '-ALFA' in prot:
            add = True
        if '_F:' in prot and prot not in self.filteredMicroproteins:
            add = False
        return add
