import os
import sys

from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd

from ..pipeline_config import PipelineStructure


class HomologyFinder(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.homologyDir = f'{self.outdir}/homology'
        self.fastaSeqsToAlignNucDir = f'{self.homologyDir}/fasta_to_align_nuc'
        self.fastaSeqsToAlignProteinDir = f'{self.homologyDir}/fasta_to_align_prot'

        self.nucMSADir = f'{self.homologyDir}/MSA_nuc'
        self.proteinMSADir = f'{self.homologyDir}/MSA_prot'
        self.check_dirs([self.homologyDir, self.fastaSeqsToAlignNucDir, self.fastaSeqsToAlignProteinDir,
                         self.nucMSADir, self.proteinMSADir])

        self.microproteinsMMFastaNucleotide = f'{self.countsDir}/mm_microproteins_nucleotide.fasta'
        self.blastedMMMicroproteinsXMLNucleotide = f'{self.homologyDir}/blasted_mm_microproteins_nucleotide.xml'
        self.blastedMMMicroproteinsXMLProtein = f'{self.homologyDir}/blasted_mm_microproteins_protein.xmls'
        self.blastNMMMicroproteinsResults = f'{self.homologyDir}/blasted_mm_microproteins_results_blastn.txt'
        self.tBlastnMMMicroproteinsResults = f'{self.homologyDir}/blasted_mm_microproteins_results_tblastn.txt'

    def get_mm_nucleotide_sequences(self):
        fasta = []
        # df = pd.read_csv(self.microproteinMappingGroups, sep='\t')
        # df = df[(df["MM_Amb"] >= self.args.minRawCounts) | (df["MM"] >= self.args.minRawCounts)]
        mps = []
        records = SeqIO.parse(self.microproteinsMM, 'fasta')
        for record in records:
            mps.append(str(record.description))
        # mps = df["microprotein"].tolist()

        dbs = os.listdir(self.translationDir)
        for db in dbs:
            files = os.listdir(f'{self.translationDir}/{db}')
            for file in files:
                if file.endswith(".split_nuc"):
                    records = SeqIO.parse(f'{self.translationDir}/{db}/{file}', 'fasta')
                    for record in records:
                        entry = str(record.description)
                        if entry in mps:
                            fasta.append(f'>{entry}\n{str(record.seq)}\n')

        with open(self.microproteinsMMFastaNucleotide, 'w') as outfile:
            outfile.writelines(fasta)

    def blastn_mm_microproteins(self):
        print(f"--Blasting smORFs nucleotide sequences against the genome.")
        cmd = (f'blastn -query {self.microproteinsMMFastaNucleotide} -subject {self.args.genome} -outfmt 5 -out'
               f' {self.blastedMMMicroproteinsXMLNucleotide}')
        os.system(cmd)

    def tblastn_mm_microproteins(self):
        print(f"--Blasting microproteins amino acid sequences against the genome.")
        cmd = (f'tblastn -query {self.select_fasta()} -subject {self.args.genome} -outfmt 5 -out '
               f'{self.blastedMMMicroproteinsXMLProtein}')
        os.system(cmd)

    def __check_blast_type(self, blast):
        if blast == 'blastn':
            blasted = self.blastedMMMicroproteinsXMLNucleotide
            out = self.blastNMMMicroproteinsResults
            aln = self.fastaSeqsToAlignNucDir
            msa = self.nucMSADir
        elif blast == 'tblastn':
            blasted = self.blastedMMMicroproteinsXMLProtein
            out = self.tBlastnMMMicroproteinsResults
            aln = self.fastaSeqsToAlignProteinDir
            msa = self.proteinMSADir
        return blasted, out, aln, msa

    def extract_aligned_sequences(self, blast='tblastn'):
        added = []
        print(f"--Extracting aligned sequences for {blast} results.")
        blasted, out, aln, msa = self.__check_blast_type(blast)
        result_handle = open(blasted)

        blast_records = NCBIXML.parse(result_handle)

        data = {'microprotein': [], 'subject': [], 'subject_sequence': [], 'evalue': [], 'bitscore': [],
                'alignment_length': [], 'coordinates': []}

        # Iterate over each BLAST record
        for blast_record in blast_records:
            # Iterate over each alignment
            for alignment in blast_record.alignments:
                # Iterate over each hsp (High-scoring Segment Pair)
                hit_id = alignment.hit_id
                # Get the hit description (additional information about the hit)
                hit_description = alignment.hit_def
                for hsp in alignment.hsps:
                    # Check if the alignment meets the filtering criteria
                    if hsp.expect < self.args.eValue and hsp.bits > self.args.bitScore:
                        # Print the subject sequence ID and description
                        query_definition = blast_record.query

                        # print("Subject ID:", alignment.hit_id)
                        # print("Description:", alignment.hit_def)
                        # # Optionally, you can print other information like e-value, bitscore, etc.
                        # print("E-value:", hsp.expect)
                        # print("Bit score:", hsp.bits)
                        # print("Alignment length:", hsp.align_length)
                        # print("Alignment:", hsp.align_seq)
                        # print()  # Add an empty line for readability
                        # print("Hit ID:", hit_id)
                        # print("Hit Description:", hit_description)
                        # hit_strand = hsp.sbjct_strand
                        hit_strand = alignment.hsps[0].frame[1]  # [1] index refers to the strand information
                        # Convert the strand information to the appropriate value (1 or -1)
                        hit_strand = '+' if hit_strand > 0 else '-'
                        hit_start = hsp.sbjct_start
                        hit_end = hsp.sbjct_end
                        name = f'{hit_strand}{alignment.hit_id}'
                        if name not in added:
                            data['coordinates'].append(f':{hit_start}-{hit_end}')
                            data['microprotein'].append(query_definition)
                            data['subject'].append(f'{hit_strand}{alignment.hit_id}')
                            data['subject_sequence'].append(hsp.sbjct)
                            data['evalue'].append(hsp.expect)
                            data['bitscore'].append(hsp.bits)
                            data['alignment_length'].append(hsp.align_length)
                            added.append(name)

                        # print("Hit Start:", hit_start)
                        # print("Hit End:", hit_end)
        df = pd.DataFrame(data=data)
        df.to_csv(out, sep='\t', index=False)

    def prepare_alignment_files(self, blast='tblastn'):
        print(f"--Preparing alignments for {blast} results.")
        blasted, out, alns, msa = self.__check_blast_type(blast)
        df = pd.read_csv(out, sep='\t')
        mps = df["microprotein"].tolist()
        chroms = df["subject"].tolist()
        subject_seqs = df["subject_sequence"].tolist()
        coords = df["coordinates"].tolist()

        mp_sequences = self.__get_mp_sequences(mp_entries=mps, blast=blast)

        fastas = {}

        for i, mp in enumerate(mps):
            if mp not in fastas:
                fastas[mp] = [f'>{mp}\n{mp_sequences[mp]}\n']
            entry = f'{chroms[i]}{coords[i]}'
            record = f'>{entry}\n{subject_seqs[i].replace("-", "")}\n'
            fastas[mp].append(record)

        for mp in fastas:
            with open(f'{alns}/{mp}_{len(fastas[mp])-1}_homologs.fasta', 'w') as handler:
                handler.writelines(fastas[mp])

    def __get_mp_sequences(self, mp_entries, blast):
        if blast == 'tblastn':
            fasta = self.microproteinsMM
        elif blast == 'blastn':
            fasta = self.microproteinsMMFastaNucleotide
        mp_sequences = {}
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            entry = str(record.description)
            if entry in mp_entries:
                mp_sequences[entry] = str(record.seq)
        return mp_sequences

    def perform_msa(self, blast='tblastn'):
        print(f"--Performing multiple sequence alignment for {blast} results.")
        blasted, out, aln, msa = self.__check_blast_type(blast)
        files = os.listdir(aln)
        for file in files:
            cmd = (f'msa4u -fa {aln}/{file} -o-aln {msa}/{file}_aligned.fasta -o '
                   f'{msa}/{file}_alignment.pdf')
            os.system(cmd)

