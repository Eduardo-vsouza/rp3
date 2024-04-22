import sys

import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SeqIO


class BlastParser:
    def __init__(self, xml):
        self.xml = xml
        self.annotatedPeptides = []

        self.conservation = {}
        self.paralogs = {}

    def __add_paralog(self, query, hit, hsp):
        name = f'{hit}:{hsp.sbjct_start}-{hsp.sbjct_end}'
        if query not in self.paralogs:
                self.paralogs[query] = [name]
        else:
            if name not in self.paralogs[query]:
                self.paralogs[query].append(name)

    def save_paralogs(self, outfile):
        data = {'smorfs': [], 'homologs': [], 'number_of_homologs': []}
        for smorf in self.paralogs:
            data['smorfs'].append(smorf)
            data['homologs'].append(','.join(self.paralogs[smorf]))
            data['number_of_homologs'].append(len(self.paralogs[smorf]))
        df = pd.DataFrame(data=data)
        df.to_csv(outfile, sep='\t', index=False)


    def parse(self, evalue=0.001, score=50, pc_id=100, qcov=100, conservation=False, smorfs=None, paralogs=False):
        if smorfs is not None:
            entries = self.__filter_smorfs(fasta=smorfs)
        for record in NCBIXML.parse(open(self.xml)):
            if record.alignments:  # skip queries with no matches
                for align in record.alignments:

                    for hsp in align.hsps:
                        if hsp.expect <= evalue and hsp.score >= score:
                            # print(vars(align))
                            # print(hsp.align_length - hsp.identities)
                            # print(vars(hsp))
                            pc_identity = (hsp.identities / hsp.align_length) * 100
                            # print(hsp.identities)
                            # print()
                            query_cover = (hsp.align_length / record.query_length) * 100
                            # print(hsp.query)
                            # print(hsp.align_length)
                            # print(record.query_length)
                            # print(pc_identity)
                            # print(query_cover)
                            # if (sbjct_end - sbjct_start) >
                            # print([var for var in vars(hsp)])
                            if pc_identity >= pc_id and query_cover >= qcov:
                                if paralogs:
                                    self.__add_paralog(query=record.query, hit=align.hit_def, hsp=hsp)
                                # self.annotatedPeptides.append(hsp.query)
                                # print(record.query)
                                self.annotatedPeptides.append(record.query)
                                if conservation:
                                    add = True
                                    if smorfs is not None:
                                        if record.query in entries:
                                            add = True
                                        else:
                                            add = False
                                    if add:
                                        species = ' '.join(align.hit_def.split(" ")[:2])
                                        if species not in self.conservation:
                                            self.conservation[species] = []
                                        if record.query not in self.conservation[species]:
                                            self.conservation[species].append(record.query)
                            # else:
                            #     if query_cover >= 100:
                            #         print(hsp.query_start)
                            #         print(hsp.align_length)
                            #         print(hsp.identities)
                            #         print(hsp.align_length - hsp.identities)
                            #         if hsp.query_start == 1:
                            #             if hsp.align_length - hsp.identities <= 1:
                            #                 print(hsp.query)
                            #                 self.annotatedPeptides.append(hsp.query)
        return self.annotatedPeptides

    @staticmethod
    def __filter_smorfs(fasta):
        entries = []
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            entries.append(str(record.description))
        return entries

    def filter_results(self, results, output):
        df = pd.read_csv(results, sep='\t')
        df = df[df["Fixed Peptides"].isin(self.annotatedPeptides) == False]
        df.to_csv(output, sep='\t', index=False)

    def filter_peptide_fasta(self, fasta, output):
        print("\n")
        print("Removing peptides based on the Blast results.\n")
        to_write = []
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            entry = str(record.description)
            seq = str(record.seq)
            if entry not in self.annotatedPeptides and seq.startswith("M"):
            # if entry not in self.annotatedPeptides:
                to_write.append(f'>{str(record.description)}\n{seq}\n')
        with open(output, 'w') as outfile:
            outfile.writelines(to_write)
        print(f"Blasted peptides removed from input file and saved to {output}\n")

    def save_conserved(self, output):
        data = {'Species': [], 'Homologs': []}
        for sp in self.conservation:
            data['Species'].append(sp)
            data['Homologs'].append(len(self.conservation[sp]))
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)

    def create_spreadsheet(self, output):
        data = {'species': [], 'smorf': []}
        for sp in self.conservation:
            for smorf in self.conservation[sp]:
                data['species'].append(sp)
                data['smorf'].append(smorf)
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("This script takes as input a Blast XML file and a fasta file. It then removes from the fasta file "
              "any sequence that has a hit in the Blast results with 100% identity and Query coverage. \n"
              "usage: .py <blast_xml> <fasta_to_filter> <output_fasta>")
    else:
        data = BlastParser(xml=sys.argv[1])
        data.parse(evalue=0.001, score=50)
        data.filter_peptide_fasta(fasta=sys.argv[2], output=sys.argv[3])

