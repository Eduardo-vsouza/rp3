import os

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure


class SpecComparison:
    def __init__(self, args):
        self.args = args

        self.filteredMicroproteins = {}
        self.__read_fasta_results()

        self.proteins = {}
        self.peptides = {}

        self.proteinsUnique = {}

        self.mergedProteins = {'protein': [], 'sequence': []}

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
            proteins = df["proteinIds"].tolist()
            peptides = df["peptide"].tolist()
            for protein, peptide in zip(proteins, peptides):
                protlist = protein.split(",")
                for prot in protlist:
                    # if '_F:' in prot and 'ANNO' not in protein:
                    #     add = True
                    #     print(prot)
                    # else:
                    #     add = False
                    # print(prot)
                    add = self.__check_filters(prot, protein)
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
                        if ',' not in protein:  # check if they are unique
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


    def create_data_frame(self):
        for group in self.args.groups:
            if group not in self.mergedProteins:
                self.mergedProteins[group] = []
        for group in self.args.groups:
            self.mergedProteins[f'{group}_peptides'] = []
        for protein in self.proteins:
            self.mergedProteins['protein'].append(protein)
            self.mergedProteins['sequence'].append(self.filteredMicroproteins[protein])
            for group in self.args.groups:
                self.mergedProteins[f'{group}_peptides'].append(','.join(self.peptides[protein][group]))
            if len(self.filteredMicroproteins[protein]) > 150:
                print(protein)
            for group in self.proteins[protein]:
                self.mergedProteins[group].append(self.proteins[protein][group])
        PipelineStructure.check_dirs([self.args.outdir])
        df = pd.DataFrame(data=self.mergedProteins)
        df.to_csv(f'{self.args.outdir}/group_comparison.xlsx', sep='\t', index=False)
        df.to_csv(f'{self.args.outdir}/group_comparison.csv', sep='\t', index=False)


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
        self.filteredMicroproteins['E4orf6/7-ALFA'] = 'MAANSTSDLHTPGTQLSVADIIVITVYFALNVAVGIWSSCRASRNTVNGYFLAGRDMTWWPIGASLFASSEGSGLFIGLAGSGAAGGLAVAGFEWNATYVLLALAWVFVPIYISSEIVTLPEYIQKRYGGQRIRMYLSVLSLLLSVFTKISLDLYAGALFVHICLGWNFYLSTILTLGITALYTIAGGLAAVIYTDALQTLIMVVGAVILTIKAFDQIGGYGQLEAAYAQAIPSRTIANTTCHLPRTDAMHMFRDPHTGDLPWTGMTFGLTIMATWYWCTDQVIVQRSLSARDLNHAKAGSILASYLKMLPMGLIIMPGMISRALFPDDVGCVVPSECLRACGAEVGCSNIAYPKLVMELMPIGLRGLMIAVMLAALMSSLTSIFNSSSTLFTMDIWRRLRPRSGERELLLVGRLVIVALIGVSVAWIPVLQDSNSGQLFIYMQSVTSSLAPPVTAVFVLGVFWRRANEQGAFWGLIAGLVVGATRLVLEFLNPAPPCGEPDTRPAVLGSIHYLHFAVALFALSGAVVVAGSLLTPPPQSVQIENLTWWTLAQDVPLGTKAGDGQTPQKHAFWARVCGFNAILLMCVNIFFYAYFA'
    def separate_up_regulated(self):
        df = pd.read_csv(f'{self.args.outdir}/group_comparison.xlsx', sep='\t')
        # df = pd.read_csv(f'{self.args.outdir}/group_comparison.csv', sep='\t')

        # for group in self.args.groups:
        #     if group != self.highlightGroup:

        filtered_df = df[(df["IP"] > df['UI']) & (df['IP'] > df['102'])]
        filtered_df.to_csv(f'{self.args.outdir}/ip_upregulated.xlsx', sep='\t', index=False)
        filtered_df.to_csv(f'{self.args.outdir}/ip_upregulated.csv', sep='\t', index=False)

        unique = df[(df['UI'] == 0) & (df['102'] == 0) & (df['IP'] > 0)]
        unique.to_csv(f'{self.args.outdir}/ip_unique.xlsx', sep='\t', index=False)
        unique.to_csv(f'{self.args.outdir}/ip_unique.csv', sep='\t', index=False)


    def __check_filters(self, prot, proteins):
        add = True
        if self.args.predictedOnly:
            if '_F:' not in prot:
                add = False
        if prot not in self.filteredMicroproteins:
            add = False
        if 'OS=' in proteins and self.args.predictedOnly:
            add = False
        if 'ALFA' in prot and 'ANNO' not in prot:
            add = True

        return add
