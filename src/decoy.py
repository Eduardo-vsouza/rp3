#!/usr/bin/python3

import os
import sys

from Bio import SeqIO


class Decoy(object):
    def __init__(self, db):
        self.df = db
        self.target = []

        self.seqs, self.entries = self.__get_seqs()
        self.reversed = []
        self.path = sys.path[0]


    def __get_seqs(self):
        seqs = []
        entries = []
        records = SeqIO.parse(self.df, 'fasta')
        for record in records:
            self.target.append(f'>{str(record.description)}\n{str(record.seq)}\n')
            seqs.append(record.seq)
            entries.append(str(record.description))
        return seqs, entries

    def reverse_sequences(self):
        """ Reverses the amino acid sequence of a protein, except for the aa in the c-terminal. """
        for seq in self.seqs:
            # cterminus = seq[len(seq)-1]
            # to_reverse = seq[:-1]
            to_reverse = seq
            reversed = to_reverse[::-1]
            # reversed += cterminus
            self.reversed.append(reversed)
        return self

    def add_contaminants(self):
        seqs = []
        records = SeqIO.parse(f'{self.path}/data/contaminants.txt', 'fasta')
        for record in records:
            seqs.append(f'>{record.description.replace(",", "_").replace(" ", "_")}\n{record.seq}\n')
        return seqs

    def to_fasta(self, output, pattern='rev', merge=True):
        out = []
        for i in range(len(self.reversed)):
            string = f">{pattern}_{self.entries[i]}\n{self.reversed[i]}\n"
            out.append(string)
        seqs = self.add_contaminants()
        for seq in seqs:
            to_add = seq.replace("B", "")
            to_add = to_add.replace("X", "")
            to_add = to_add.replace("Z", "")
            out.append(to_add)
        if merge:
            with open(output, 'w') as fa:
                fa.writelines(self.target)
                fa.writelines(out)
        else:
            with open(output, 'w') as fa:
                fa.writelines(out)
        return self


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: generate_decoy.py <target_database> <output_decoy_database> <decoy_pattern>\n'
              'decoy_pattern default: rev')
    else:
        data = Decoy(db=sys.argv[1])
        data.reverse_sequences().to_fasta(output=sys.argv[2], pattern=sys.argv[3])



