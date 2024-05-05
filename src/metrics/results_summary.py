import os
import sys

from Bio import SeqIO
import pandas as pd

from ..pipeline_config import PipelineStructure


class MicroproteinCombiner(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.microproteins = {}

    def gather_microprotein_info(self):
        self.show_progress(step=1)
        fasta = self.select_fasta()
        records = SeqIO.parse(f'{fasta}', 'fasta')
        for record in records:
            self.microproteins[str(record.description)] = {'sequence': str(record.seq)}
        # self.__add_signalp()
        self.__add_conservation()
        self.__add_riboseq_coverage()
        self.__add_peptide_data()
        self.__add_smorf_classes()

    def __add_signalp(self):
        self.show_progress(step=2)

        signalp = f'{self.signalPDir}/processed_entries.fasta'
        for mp in self.microproteins:
            self.microproteins[mp]['signalP'] = False
        if os.path.exists(signalp):
            records = SeqIO.parse(signalp, 'fasta')
            for record in records:
                entry = str(record.description)
                if entry in self.microproteins:
                    self.microproteins[entry]['signalP'] = True

    def __add_conservation(self):
        self.show_progress(step=3)

        file = f'{self.phyloDir}/smorfs_entries_per_species.xls'
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            species, smorfs = df["species"].tolist(), df["smorf"].tolist()
            for mp in self.microproteins:
                if 'conservation' not in self.microproteins[mp]:
                    self.microproteins[mp]['conservation'] = ''
            for sp, smorf in zip(species, smorfs):
                if smorf in self.microproteins:

                    self.microproteins[smorf]['conservation'] += f',{sp}'
            for mp in self.microproteins:
                self.microproteins[mp]['conservation'] = self.microproteins[mp]['conservation'][1:]

                self.microproteins[mp]['num_species_conserved'] = len(self.microproteins[mp]['conservation'].split(","))

    def __add_riboseq_coverage(self):
        self.show_progress(step=4)

        if os.path.exists(self.mappingGroups):
            df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
            smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
            for smorf, group in zip(smorfs, groups):
                self.microproteins[smorf]['Ribo-seq mapping group'] = group
        for smorf in self.microproteins:
            if 'Ribo-seq mapping group' not in self.microproteins[smorf]:
                self.microproteins[smorf]['Ribo-seq mapping group'] = 'No coverage'

            # rpkms = pd.read_csv(self.mappingGroupsRPKMs, sep='\t')
            # print(self.microproteins)
            # entries = rpkms["Gene"].tolist()
            # for i, entry in enumerate(entries):
            #     if self.microproteins[entry]['Ribo-seq mapping group'] != 'No coverage':
            #         rpkm = rpkms.at[i, self.microproteins[entry]['Ribo-seq mapping group']]
            #         self.microproteins[entry]['Ribo-seq RPKM'] = rpkm
            #     else:
            #         self.microproteins[entry]['Ribo-seq RPKM'] = 0

    def __add_peptide_data(self):
        self.show_progress(step=5)
        folder = self.postProcessDir
        if os.path.exists(self.rescoredMicroproteinsFasta):
            folder = self.rescorePostProcessDir
        groups = os.listdir(folder)

        for group in groups:
            #
            # dbs = os.listdir(f'{folder}/{group}')
            # for db in dbs:
            db = group

            dbdir = f'{folder}/{group}/{db}'
            groupdir = f'{folder}/{group}/{group}_target_decoy_database.fasta'
            # groupdir = dbdir
            groupdir = f'{folder}/{group}'
            if folder != self.rescorePostProcessDir:
                groupdir = f'{folder}/{group}/db'
            col = f'MS Spectral counts'
            utp_col = f'true UTPs'
            tp_col = 'MS peptides'
            for mp in self.microproteins:
                if 'total spec count' not in self.microproteins[mp]:
                     self.microproteins[mp]['total spec count'] = 0
                if col not in self.microproteins[mp]:
                    self.microproteins[mp][col] = 0
                    self.microproteins[mp][utp_col] = []
                    self.microproteins[mp][tp_col] = []
                    self.microproteins[mp]["num true_UTPs"] = []
                    # self.microproteins[mp][group] = False

            df = pd.read_csv(f'{groupdir}/peptides_fixed.txt', sep='\t')
            df = df[df["proteinIds"].str.contains("_ANNO") == False]
            df = df[df["q-value"] <= 0.01]
            peptides, proteins = df["peptide"].tolist(), df["proteinIds"].tolist()
            for prot, pep in zip(proteins, peptides):
                prot_list = prot.split(",")
                for smorf in prot_list:
                    if smorf in self.microproteins:
                        # self.microproteins[smorf][group] = True
                        self.microproteins[smorf][col] += 1
                        self.microproteins[smorf]['total spec count'] += 1
                        if pep not in self.microproteins[smorf][tp_col]:
                            self.microproteins[smorf][tp_col].append(pep)
                        if pep not in self.microproteins[smorf][utp_col] and ',' not in prot:
                            self.microproteins[smorf][utp_col].append(pep)
            for mp in self.microproteins:
                self.microproteins[mp][f"num true_UTPs"] = len(self.microproteins[mp][utp_col])
        for mp in self.microproteins:
            for col in self.microproteins[mp]:
                if ' UTPs' in col:
                    self.microproteins[mp][col] = ','.join(self.microproteins[mp][col])

    def __add_smorf_classes(self):
        if os.path.exists(f'{self.orfClassDir}/predicted_nonhomolog_smorfs_annotation'):
            df = pd.read_csv(f'{self.orfClassDir}/predicted_nonhomolog_smorfs_annotation', sep='\t',
                             usecols=[0,1,2], names=['smorf', 'class', 'overlapped_gene'])
            smorfs, classes, genes = df["smorf"].tolist(), df["class"].tolist(), df["overlapped_gene"].tolist()
            for i, smorf in enumerate(smorfs):
                if smorf in self.microproteins:
                    # if 'smorf_class' not in self.microproteins[smorf]:
                    self.microproteins[smorf]['smorf_class'] = classes[i]
                    self.microproteins[smorf]['overlapped_gene'] = genes[i]
            for smorf in self.microproteins:

                if 'smorf_class' not in self.microproteins[smorf]:
                    self.microproteins[smorf]['smorf_class'] = 'Intergenic'
                    self.microproteins[smorf]['overlapped_gene'] = None


    def show_progress(self, step):
        step_n = {1: 'Gathering data', 2: 'signalP', 3: 'Conservation', 4: 'Ribo-Seq', 5: 'Peptide data'}
        steps = {'Gathering data': '', 'signalP': '', 'Conservation': '', 'Ribo-Seq': '', 'Peptide data': ''}
        steps[step_n[step]] = '↓'
        prog = f''
        prog_lower = f''
        for s in steps:
            # print(s)
            prog += f'{"."*14}'
            ind = (len(prog)-len(s)/2)-1
            # print(ind)
            # print(prog[5])
            # ind = len(prog)-len(s)
            prog_list = list(prog)
            prog_list[int(ind)] = steps[s]
            # print(prog_list)
            prog = ''.join(prog_list)
            prog_lower += f'{" "*(14-len(s))}{s}'
            # merged = f'{prog}\n{prog_lower}\n'
        print(prog, end='\n')
        print(prog_lower, end='\n')
        print('\n')

        # print(merged, end='\n'

    def save(self):
        data = {'microprotein': []}
        for mp in self.microproteins:
            data['microprotein'].append(mp)
            for col in self.microproteins[mp]:
                # print(self.microproteins[mp][col])
                # print(len(self.microproteins[mp][col]))
                if col not in data:
                    data[col] = []
                data[col].append(self.microproteins[mp][col])
                # print(len(self.microproteins[mp][col]))
                # print(col)
                # data[col].append(self.microproteins[mp][col])
        # for i in data:
        #     print(i, len(data[i]))
        # for col in data:
        #     print(len(data[col]))
        #     print(col, '\n')
        df = pd.DataFrame(data=data)
        df.to_csv(f'{self.mergedResults}/microproteins_summary.txt', sep='\t', index=False)
        df.to_csv(f'{self.mergedResults}/microproteins_summary.xls', sep='\t', index=False)
        df.to_csv(f'{self.mergedResults}/microproteins_summary.xlsx', sep='\t', index=False)









