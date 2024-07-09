import os
import sys

from ..pipeline_config import PipelineStructure


class StringTieAssembly(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.print_row(word="Assembly mode")
        self.mode = 'rna'
        self.strandness = args.strandness
        self.gtf = args.gtf
        self.threads = args.threads
        # self.intermediateAssembliesDir = f'{self.assemblyDir}/intermediate_assemblies'
        # self.check_dirs([self.intermediateAssembliesDir])

        self.check_dirs([self.rnaDir, self.rnaTrimmedDir, self.assemblyDir, self.rnaAlnDir, self.rnaQCDir,
                         self.interAssemblyDir])
        self.check_dirs([self.genomeIndexDir])


        # self.bamDir = f''
        self.readPairs = {}
        self.__get_read_pairs()

    def __get_read_pairs(self):
        reads = os.listdir(self.args.fastq)

        for read in reads:
            path = f'{self.args.fastq}/{read}'
            if self.args.libType == 'paired':
                if f"_{self.args.pairedPattern}2_" not in read:
                    self.readPairs[path] = path.replace(f"_{self.args.pairedPattern}1_", f"_{self.args.pairedPattern}2_")
            else:
                self.readPairs[path] = None

    def clean_reads(self):
        for read in self.readPairs:
            if self.readPairs[read] is not None:
                out = f'{self.rnaTrimmedDir}/{read.split("/")[-1].split(".")[0]}_val_1.fq.gz'
                run = self.verify_checkpoint(outfile=out, step=f"read trimming for {read}")
                if run:

                    print(f"--Trimming read pairs {read} and {self.readPairs[read]}.")
                    cmd_trim = f'{self.toolPaths["trim_galore"]} -q 20 --fastqc_args "--outdir {self.rnaQCDir} --threads 8" --gzip ' \
                               f'--cores 8 -o {self.rnaTrimmedDir} --paired {read} {self.readPairs[read]}'
                    os.system(cmd_trim)
                    print(f"--Finished trimming read pair.")
        print(f"--Completed read trimming.")

    def __get_trimmed_pairs(self):
        pairs = {}
        files = os.listdir(self.rnaTrimmedDir)
        for file in files:
            if f"_{self.args.pairedPattern}2_" not in file and not file.endswith(".txt"):
                pairs[file] = file.replace(f"_{self.args.pairedPattern}1_", f"_{self.args.pairedPattern}2_").replace("_val_1", "_val_2")
        return pairs

    def remove_contaminants(self):
        self.print_row(word="Contaminant removal")
        pairs = self.__get_trimmed_pairs()
        files = os.listdir(self.rnaTrimmedDir)
        for file in files:
            if not file.endswith(".txt") and file in pairs and not os.path.isdir(f'{self.rnaTrimmedDir}/{file}'):
                run = self.verify_checkpoint(outfile=f'{self.rnaNoContDir}/no_contaminant_{file}Unmapped.out.mate1',
                                             step="contaminant removal")
                if file.endswith(".gz"):
                    zcat = f' --readFilesCommand zcat'
                else:
                    zcat = ''
                if run:
                    print(f"--Removing contaminants from {file}.")
                    cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --outReadsUnmapped Fastx --genomeDir '
                           f'{self.args.contIndex} --runThreadN {self.args.threads} --readFilesIn'
                           f' {self.rnaTrimmedDir}/{file} {self.rnaTrimmedDir}/{pairs[file]} '
                           f'--outFileNamePrefix {self.rnaNoContDir}/no_contaminant_{file}{zcat}')
                    os.system(cmd)
                    print(f"--Finished removing contaminants from {file}.")
        print(f"--Finished removing contaminants.")

    def __get_decont_pairs(self):
        pairs = {}
        files = os.listdir(self.rnaNoContDir)
        for file in files:
            if file.endswith("mate1"):
                # if f"_{self.args.pairedPattern}2_" not in file and not file.endswith(".txt"):
                pairs[file] = file.replace("mate1", "mate2")
        return pairs

    def align_to_genome(self):
        self.print_row(word="Genome alignment")
        pairs = self.__get_decont_pairs()
        files = os.listdir(self.rnaNoContDir)
        for file in files:
            if 'Unmapped.out.mate1' in file or file.endswith(".fastq"):
                out = f'{self.rnaAlnDir}/aligned_to_genome_{file}Aligned.sortedByCoord.out.bam'
                run = self.verify_checkpoint(outfile=out, step=f"alignment of clean RNA-Seq reads to the genome")
                if run:
                    if file.endswith(".gz"):
                        gz = ' --readFilesCommand zcat'
                    else:
                        gz = ''
                    filepath = f'{self.rnaNoContDir}/{file}'
                    print(f"Aligning RPF reads from {file} to {self.args.genomeIndex}")
                    cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --genomeDir {self.args.genomeIndex}'
                           f' --runThreadN {self.args.threads} --readFilesIn {filepath} {self.rnaNoContDir}/{pairs[file]} --outFileNamePrefix '
                           f'{self.rnaAlnDir}/aligned_to_genome_{file} --outSAMtype BAM SortedByCoordinate'
                           f' --outSAMattributes All{gz}')
                    # os.system(cmd)
                    print(f'Finished aligning {file}.')

    def assemble_transcriptome(self):
        self.print_row(word="Transcriptome assembly")
        files = os.listdir(self.rnaAlnDir)
        for file in files:
            if file.endswith(".bam"):
                print(f'--Performing the assembly for {file}\n')
                cmd = f'stringtie --{self.strandness} -G {self.gtf} -o {self.interAssemblyDir}/{file}.gtf -m 50 -p' \
                      f' {self.args.threads} -A {self.interAssemblyDir}/{file}_gene_abundances.txt {self.rnaAlnDir}/{file}'
                os.system(cmd)
                print(f'--Transcriptome assembled at {self.interAssemblyDir}/{file}\n')

    def merge_transcriptomes(self):
        file_list = ''
        files = os.listdir(f'{self.interAssemblyDir}')
        for file in files:
            if file.endswith(".gtf") and os.path.getsize(f'{self.interAssemblyDir}/{file}') > 1000:
                file_list += f' {self.interAssemblyDir}/{file}'
        print(f"--Merging StringTie assemblies")
        cmd = (f'{self.toolPaths["stringtie"]} --merge -G {self.args.gtf} -o {self.assemblyDir}/merged_assembly.gtf '
               f'-m 50 -p {self.args.threads}{file_list}')
        os.system(cmd)
        print(f"--Finished merging StringTie assemblies")

    #
    # def check_index(self):
    #     if self.args.index is None:
    #         print("Index not provided. Generating index for HISAT2 from genome fasta file.")
    #         cmd_splice = f'{self.toolPaths["hisat2_extract_splice_sites"]} {self.gtf} > {self.spliceSites}'
    #         print("--Extracting splice sites from GTF file.")
    #         os.system(cmd_splice)
    #
    #         cmd_exon = f'{self.toolPaths["hisat2_extract_exons"]} {self.gtf} > {self.exons}'
    #         print("--Extracting exons from GTF file.")
    #         os.system(cmd_exon)
    #
    #         hisat = f'{self.toolPaths["hisat2-build"]}'
    #         cmd = (f'{hisat} -p {self.args.threads} {self.args.genome} {self.genomeIndex} --exon {self.exons} '
    #                f'--ss {self.spliceSites}')
    #         print("Generating genome index.")
    #         os.system(cmd)
    #     else:
    #         self.genomeIndex = self.args.index
    #
    # def align(self):
    #     groups = os.listdir(self.args.reads_folder)
    #     for group in groups:
    #         outdir = f'{self.rnaAlnDir}/{group}'
    #         self.check_dirs([outdir])
    #
    #         groupdir = f'{self.args.reads_folder}/{group}'
    #         files = os.listdir(groupdir)
    #
    #         for file in files:
    #             if '_2.fastq' not in file:
    #                 sorted_dir = f'{outdir}/sorted_alignments'
    #                 self.check_dirs([sorted_dir])
    #             if self.args.lib_type == 'single':
    #                 cmd = (f'{self.toolPaths["hisat2"]} -x {self.genomeIndex} --threads {self.args.threads} --dta '
    #                        f'--rna-strandness {self.args.strandness} -U {self.args.reads_folder}/{group}/{file} '
    #                        f'-S {outdir}/{file}.sam')
    #                 os.system(cmd)
    #             else:
    #                 if '_2.fastq' not in file:
    #                     cmd = (f'{self.toolPaths["hisat2"]} -x {self.genomeIndex} --threads {self.args.threads} --dta '
    #                            f'--rna-strandness {self.args.strandness.upper()} -1 {self.args.reads_folder}/{group}/{file} '
    #                            f'-2 {self.args.reads_folder}/{group}/{file.replace("_1.fastq", "_2.fastq")} -S '
    #                            f'{outdir}/{file}.sam')
    #                     os.system(cmd)
    #             if '_2.fastq' not in file:
    #
    #                 cmd_sam = (f'{self.toolPaths["samtools"]} view --threads {self.args.threads} -bS {outdir}/{file}.sam '
    #                            f'| {self.toolPaths["samtools"]} sort -o {sorted_dir}/{file}_sorted.bam')
    #                 os.system(cmd_sam)
    #                 print(f"Finished aligning {file}. Results at {sorted_dir}/{file}_sorted.bam")
    #
    #
    # def assemble_transcriptomes(self):
    #     groups = os.listdir(self.rnaAlnDir)
    #     print(groups)
    #     for group in groups:
    #         if os.path.isdir(f'{self.rnaAlnDir}/{group}'):
    #             aln_dir = f'{self.rnaAlnDir}/{group}/sorted_alignments'
    #
    #             files = os.listdir(aln_dir)
    #             print(files)
    #             groupdir = f'{self.assemblyDir}/{group}'
    #             inter_groupdir = f'{self.assemblyDir}/{group}/intermediate_assemblies'
    #             self.check_dirs([groupdir, inter_groupdir])
    #
    #
    #             for file in files:
    #                 print(f'Performing the assembly for {file}\n')
    #                 cmd = f'stringtie --{self.strandness} -G {self.gtf} -o {inter_groupdir}/{file}.gtf -m 50 -p' \
    #                       f' {self.threads} -A {inter_groupdir}/{file}_gene_abundances.txt {aln_dir}/{file}'
    #                 os.system(cmd)
    #                 print(f'Transcriptome assembled at {inter_groupdir}/{file}\n')
    #
    # def merge_assemblies(self):
    #     groups = os.listdir(self.assemblyDir)
    #     for group in groups:
    #         if os.path.isdir(f'{self.assemblyDir}/{group}'):
    #             inter_dir = f'{self.assemblyDir}/{group}/intermediate_assemblies'
    #             files = os.listdir(inter_dir)
    #
    #             intermeds = []
    #             for file in files:
    #                     if file.endswith('.gtf'):
    #                         intermeds.append(f'{inter_dir}/{file}')
    #             file_list = ''
    #             for file in intermeds:
    #                 file_list += f' {file}'
    #             print(f'Merging the intermediate transcriptome assemblies into a single GTF file.\n')
    #             cmd = f'stringtie --merge -G {self.gtf} -o {self.assemblyDir}/{group}/merged_assembly.gtf -m 50 -p ' \
    #                   f'{self.threads}{file_list}'
    #             os.system(cmd)
    #             print(f'Transcriptome assembly done. Check results at {self.assemblyDir}/{group}/merged_assembly.gtf\n')
    #

