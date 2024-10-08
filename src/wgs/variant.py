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
        genome_index = f'{self.args.genome}.fai'
        run = self.verify_checkpoint(outfile=genome_index, step="genome indexing")
        if run:
            cmd = f'{self.toolPaths["samtools"]} faidx {self.args.genome}'
            os.system(cmd)
        cmd = f'bwa index {self.args.genome}'
        os.system(cmd)

    def align_reads(self):
        print(f"--Aligning reads to the genome")
        print(self.args.genome)
        reads = os.listdir(self.wgsTrimmeddir)
        i = 0
        for read in reads:
            if read.endswith("fq.gz"):
                prefix = f'{read.split(".")[0]}'
                name1 = f'{prefix}.fq.gz'
                original1 = f'{prefix.replace("_val_1", "")}.fastq.gz'
                print(name1, original1)
                print(self.readPairs)
                if original1 in self.readPairs:
                    i += 1
                    pair = f'{self.readPairs[original1].split(".")[0]}_val_2.fq.gz'
                    sample = f'read{i}'
                    print(f"--Aligning {name1} and {pair} to the genome")
                    print(self.args.genome)
                    cmd = (f"{self.toolPaths['bwa']} mem -M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' "
                           f"-t {self.args.threads} "
                           f"{self.args.genome} {self.wgsTrimmeddir}/{name1} {self.wgsTrimmeddir}/{pair} "
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
                       f' -o {self.wgsBamDir}/{file.replace(".sam", "_sorted.bam")}')
                os.system(cmd)
        print(f"--Finished generating bam files")

    def mark_duplicates(self):
        print(f"--Marking read duplicates")
        files = os.listdir(self.wgsBamDir)
        for file in files:
            if file.endswith("bam"):
                print(f"--Marking duplicates for {file}")
                cmd = (f'{self.toolPaths["gatk"]} MarkDuplicates -I {self.wgsBamDir}/{file} -O'
                       f' {self.wgsDeduplicatedBamDir}/{file.replace(".bam", "deduplicated.bam")} '
                       f'-M {self.wgsDeduplicatedBamDir}/{file.replace(".bam", "deduplicated.metrics")}')
                os.system(cmd)
        print(f"--Finished marking read duplicates")

    def recalibrate_base_score(self):
        """
        --known-sites will help differentiating between common germline variants and somatic mutations
        and known indels from sequencing errors
        """
        files = os.listdir(self.wgsDeduplicatedBamDir)
        known_sites = ''
        for site in self.args.knownSites:
            known_sites += f' --known-sites {site}'

        for file in files:
            if file.endswith(".bam"):
                recal_file = f'{self.wgsRecalDir}/{file}_recalData.table'
                run = self.verify_checkpoint(outfile=recal_file,
                                             step="recalibration table generation")
                if run:
                    print(f"--Generating recalibration table for {file}")
                    cmd = (f'{self.toolPaths["gatk"]} BaseRecalibrator -I {self.wgsDeduplicatedBamDir}/{file} '
                           f'-R {self.args.genome}{known_sites} -O {recal_file}')
                    os.system(cmd)

                out_recal_bam = f'{self.wgsRecalDir}/{file.replace(".bam", "recalibrated.bam")}'
                run = self.verify_checkpoint(outfile=out_recal_bam,
                                             step="base score recalibration")
                if run:
                    print(f"--Recalibrating base scores")
                    cmd = (f'{self.toolPaths["gatk"]} ApplyBQSR -R {self.args.genome} '
                           f'-I {self.wgsDeduplicatedBamDir}/{file}'
                           f'--bqsr-recal-file {recal_file} -O {out_recal_bam}')
                    os.system(cmd)

                    print(f"--Indexing recalibrated bam file")
                    cmd = f'{self.toolPaths["samtools"]} index {out_recal_bam}'
                    os.system(cmd)

    def run_mutect2(self):
        print(f"--Running Mutect2 to detect variants")
        files = os.listdir(self.wgsRecalDir)
        for file in files:
            if file.endswith(".bam"):
                print(f"--Running Mutect2 on {file}")
                cmd = (f'{self.toolPaths["gatk"]} Mutect2 -R {self.args.genome} '
                       f'-I {self.wgsRecalDir}/{file}'
                       f'-germline-resource {self.args.germlineResource} '
                       f'-pon {self.args.panelOfNormals} '
                       f'--f1r2-tar-gz {self.wgsMutectDir}/{file}_f1r2.tar.gz '
                       f'-O {self.wgsMutectDir}/{file}_unfiltered.vcf')
                os.system(cmd)

    def filter_vcf(self):
        ...

    def __learn_read_orientation(self):
        files = os.listdir(self.wgsMutectDir)
        for file in files:
            if file.endswith("f1r2.tar.gz"):
                cmd = (f'{self.toolPaths["gatk"]} LearnReadOrientationModel '
                       f'-I {self.wgsMutectDir}/{file} '
                       f'-O {self.wgsMutectDir}/{file}_read_orientation_model.tar.gz')
                os.system(cmd)

    def __get_pileup_summaries(self):
        ...



