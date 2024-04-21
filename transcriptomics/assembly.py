import os
import sys

from ..pipeline_config import PipelineStructure


class StringTieAssembly(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.mode = 'rna'
        self.strandness = args.strandness
        self.gtf = args.gtf
        self.threads = args.threads
        # self.intermediateAssembliesDir = f'{self.assemblyDir}/intermediate_assemblies'
        # self.check_dirs([self.intermediateAssembliesDir])

        self.check_dirs([self.rnaDir, self.assemblyDir, self.rnaAlnDir])
        self.check_dirs([self.genomeIndexDir])


        # self.bamDir = f''


    def check_index(self):
        if self.args.index is None:
            print("Index not provided. Generating index for HISAT2 from genome fasta file.")
            cmd_splice = f'{self.toolPaths["hisat2_extract_splice_sites"]} {self.gtf} > {self.spliceSites}'
            print("--Extracting splice sites from GTF file.")
            os.system(cmd_splice)

            cmd_exon = f'{self.toolPaths["hisat2_extract_exons"]} {self.gtf} > {self.exons}'
            print("--Extracting exons from GTF file.")
            os.system(cmd_exon)

            hisat = f'{self.toolPaths["hisat2-build"]}'
            cmd = (f'{hisat} -p {self.args.threads} {self.args.genome} {self.genomeIndex} --exon {self.exons} '
                   f'--ss {self.spliceSites}')
            print("Generating genome index.")
            os.system(cmd)
        else:
            self.genomeIndex = self.args.index

    def align(self):
        groups = os.listdir(self.args.reads_folder)
        for group in groups:
            outdir = f'{self.rnaAlnDir}/{group}'
            self.check_dirs([outdir])

            groupdir = f'{self.args.reads_folder}/{group}'
            files = os.listdir(groupdir)

            for file in files:
                if '_2.fastq' not in file:
                    sorted_dir = f'{outdir}/sorted_alignments'
                    self.check_dirs([sorted_dir])
                if self.args.lib_type == 'single':
                    cmd = (f'{self.toolPaths["hisat2"]} -x {self.genomeIndex} --threads {self.args.threads} --dta '
                           f'--rna-strandness {self.args.strandness} -U {self.args.reads_folder}/{group}/{file} '
                           f'-S {outdir}/{file}.sam')
                    os.system(cmd)
                else:
                    if '_2.fastq' not in file:
                        cmd = (f'{self.toolPaths["hisat2"]} -x {self.genomeIndex} --threads {self.args.threads} --dta '
                               f'--rna-strandness {self.args.strandness.upper()} -1 {self.args.reads_folder}/{group}/{file} '
                               f'-2 {self.args.reads_folder}/{group}/{file.replace("_1.fastq", "_2.fastq")} -S '
                               f'{outdir}/{file}.sam')
                        os.system(cmd)
                if '_2.fastq' not in file:

                    cmd_sam = (f'{self.toolPaths["samtools"]} view --threads {self.args.threads} -bS {outdir}/{file}.sam '
                               f'| {self.toolPaths["samtools"]} sort -o {sorted_dir}/{file}_sorted.bam')
                    os.system(cmd_sam)
                    print(f"Finished aligning {file}. Results at {sorted_dir}/{file}_sorted.bam")


    def assemble_transcriptomes(self):
        groups = os.listdir(self.rnaAlnDir)
        print(groups)
        for group in groups:
            if os.path.isdir(f'{self.rnaAlnDir}/{group}'):
                aln_dir = f'{self.rnaAlnDir}/{group}/sorted_alignments'

                files = os.listdir(aln_dir)
                print(files)
                groupdir = f'{self.assemblyDir}/{group}'
                inter_groupdir = f'{self.assemblyDir}/{group}/intermediate_assemblies'
                self.check_dirs([groupdir, inter_groupdir])


                for file in files:
                    print(f'Performing the assembly for {file}\n')
                    cmd = f'stringtie --{self.strandness} -G {self.gtf} -o {inter_groupdir}/{file}.gtf -m 50 -p' \
                          f' {self.threads} -A {inter_groupdir}/{file}_gene_abundances.txt {aln_dir}/{file}'
                    os.system(cmd)
                    print(f'Transcriptome assembled at {inter_groupdir}/{file}\n')

    def merge_assemblies(self):
        groups = os.listdir(self.assemblyDir)
        for group in groups:
            if os.path.isdir(f'{self.assemblyDir}/{group}'):
                inter_dir = f'{self.assemblyDir}/{group}/intermediate_assemblies'
                files = os.listdir(inter_dir)

                intermeds = []
                for file in files:
                        if file.endswith('.gtf'):
                            intermeds.append(f'{inter_dir}/{file}')
                file_list = ''
                for file in intermeds:
                    file_list += f' {file}'
                print(f'Merging the intermediate transcriptome assemblies into a single GTF file.\n')
                cmd = f'stringtie --merge -G {self.gtf} -o {self.assemblyDir}/{group}/merged_assembly.gtf -m 50 -p ' \
                      f'{self.threads}{file_list}'
                os.system(cmd)
                print(f'Transcriptome assembly done. Check results at {self.assemblyDir}/{group}/merged_assembly.gtf\n')


