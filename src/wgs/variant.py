import os

import pandas as pd

from ..pipeline_config import PipelineStructure


class Variant(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.print_row(word="Variant detection")
        self.check_wgs_dirs()

        self.readPairs = self.__get_read_pairs()

    def __get_read_pairs(self):
        print("--Obtaining read pairs")
        read_pairs = {}
        df = pd.read_csv(self.args.metadata, sep='\t')
        reads, pairs = df["read"].tolist(), df["pair"].tolist()
        for read, pair in zip(reads, pairs):
            if read not in read_pairs:
                read_pairs[read] = pair
        return read_pairs

    def pre_process_reads(self):
        print(f"--Preprocessing WGS reads")
        reads = os.listdir(self.args.fastq)
        for read in reads:
            if read in self.readPairs:
                print(f"--Trimming reads {read} and {self.readPairs[read]}")
                cmd = (f'{self.toolPaths["trim_galore"]} --paired --fastqc --fastqc_args "-t {len(reads)}" --cores 4 -q 20'
                       f' --length 35 --stringency 3 -o {self.wgsTrimmeddir} '
                       f'{self.args.fastq}/{read} {self.args.fastq}/{self.readPairs[read]}')
                os.system(cmd)
                print(f"--Successfully trimmed reads {read} and {self.readPairs[read]}")
        print(f"--Finished preprocessing of WGS reads")

    def prepare_annotation_files(self):
        print(f"--Generating sequence dictionary")
        cmd = (f'{self.toolPaths["gatk"]} CreateSequenceDictionary -R {self.args.genome} '
               f'-O {self.args.genome.replace(".fasta", "")}.dict')
        os.system(cmd)
        print(f"--Indexing genome")
        cmd = f'{self.toolPaths["samtools"]} faidx {self.args.genome}'
        os.system(cmd)

    def align_reads(self):
        print(f"--Aligning reads to the genome")
        reads = os.listdir(self.wgsTrimmeddir)
        i = 0
        for read in reads:
            if read in self.readPairs:
                i += 1
                sample = f'read{i}'
                print(f"--Aligning {read} and {self.readPairs[read]} to the genome")
                cmd = (f"{self.toolPaths['bwa']} mem -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA' "
                       f"-t {self.args.threads} "
                       f"{self.args.genome} {self.wgsTrimmeddir}/{read} {self.wgsTrimmeddir}/{self.readPairs[read]} "
                       f"> {self.wgsSamDir}/{read}_alignedToGenome.sam")
                os.system(cmd)
        print(f"--Finished aligning reads to the genome")

    def convert_to_sorted_bam(self):
        print(f"--Converting sam files to sorted bam format")
        files = os.listdir(self.wgsSamDir)
        for file in files:
            if file.endswith(".sam"):
                print(f"--Converting {file} to sorted bam format")
                cmd = (f'{self.toolPaths["samtools"]} view -Sb {self.wgsSamDir}/{file} --threads {self.args.threads}'
                       f' | samtools sort --threads {self.args.threads}'
                       f'-o {self.wgsBamDir}/{file.replace(".sam", "_sorted.bam")}')
                os.system(cmd)
        print(f"--Finished generating bam files")

    def mark_duplicates(self):
        print(f"--Marking read duplicates")
        files = os.listdir(self.wgsBamDir)
        for file in files:
            if file.endswith("bam"):
                print(f"--Marking duplicates for {file}")
                cmd = (f'{self.toolPaths["gatk"]} MarkDuplicates -I {self.wgsBamDir}/{file} -o'
                       f' {self.wgsDeduplicatedBamDir}/{file.replace(".bam", "deduplicated.bam")} '
                       f'-M {self.wgsDeduplicatedBamDir}/{file.replace(".bam", "deduplicated.metrics")}')
                os.system(cmd)
        print(f"--Finished marking read duplicates")

    def recalibrate_base_score(self):


