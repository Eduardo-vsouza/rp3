import os
import sys

from Bio import SeqIO
#
#
# class ORFCounter:
#     def __init__(self, args):
#         self.args = args
#
#         self.outdir = args.outdir
#         self.resultsDir = f'{self.outdir}/results'
#
#     def count_smorfs(self):
#         proteins = []
#         groups = os.listdir(self.resultsDir)
#         for group in groups:
#             group_dir = f'{self.resultsDir}/{group}'
#             if os.path.isdir(group_dir):
#                 dbs = os.listdir(group_dir)
#                 for db in dbs:
#                     db_dir = f'{group_dir}/{db}'
#                     if os.path.exists(f'{db_dir}/microproteins.fasta'):
#                         records = SeqIO.parse(f'{db_dir}/microproteins.fasta', 'fasta')
#                         for record in records:
#                             if str(record.seq) not in proteins:
#                                 proteins.append(str(record.seq))
#         print(len(set(proteins)))
#
