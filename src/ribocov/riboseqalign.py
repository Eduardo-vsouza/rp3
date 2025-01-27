import os
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

from ..pipeline_config import PipelineStructure



class RiboSeqAlign(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.set_ribocov_params()
        self.adapterSequence = args.adapter

        self.index = f'{sys.path[0]}/data/STAR_indexes/hg19.star'
        self.indexCont = f'{sys.path[0]}/data/STAR_indexes/hg19_cont.star'
        self.__check_indexes()

    def trim_reads(self):
        if not self.args.skip_trimming:
            files = os.listdir(self.args.fastq)
            for file in files:
                # run = self.verify_checkpoint(outfile=f'{self.riboSeqTrimmedDir}/trimmed_{file}',
                #                              step=f'Ribocov trimming for {file}.')
                # if run:
                if file.endswith("gz"):
                    os.system(f'gzip -d {self.args.fastq}/{file}')
                # if file.endswith("fastq.gz")
                if file.endswith(".fastq") or file.endswith(".fq"):
                    print(f'Performing read trimming for {file}.')

                    full = f'{self.args.fastq}/{file.replace(".gz", "")}'
                    if self.args.trimmer == 'fastX':

                        cmd = (f'{self.args.fastx_clipper_path} -Q33 -l 20 -n -v -c -a {self.adapterSequence} -i {full} | '
                               f'{self.args.fastx_trimmer_path} -Q33 -f 1 2> {self.riboSeqTrimmedDir}/trim.log > '
                               f'{self.riboSeqTrimmedDir}/trimmed_{file}')
                        os.system(cmd)
                    elif self.args.trimmer == 'cutadapt':
                        threads = self.args.threads
                        if threads > 8:
                            threads = 8
                        cmd = (f'{self.toolPaths["cutadapt"]} -a {self.adapterSequence} -o {self.riboSeqTrimmedDir}/trimmed_{file} {full} '
                               f'-j {threads}')
                        os.system(cmd)
                    print(f"Finished trimming {file}.")

    def run_single_cutadapt(self, input_file, output_file):
        """
        Helper function to run cutadapt for a single file.
        This is what we'll submit to ProcessPoolExecutor.
        """
        # Cap cutadapt at 8 threads
        threads_for_cutadapt = min(self.args.threads, 8)

        cmd = [
            self.toolPaths["cutadapt"],
            "-a", self.adapterSequence,
            "-o", output_file,
            input_file,
            "-j", str(threads_for_cutadapt)
        ]
        print(f"[Cutadapt] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"[Cutadapt] Finished trimming {os.path.basename(input_file)}")

    def run_trimming_parallel(self):
        """
        Main trimming function that:
         1) Finds FASTQ files in self.args.fastq
         2) Decompresses .gz files if needed
         3) Submits each file to ProcessPoolExecutor for parallel cutadapt runs
        """

        if self.args.skip_trimming:
            print("Skipping trimming step, as requested.")
            return

        # Gather input files (FASTQ or FQ, possibly gzipped)
        files = os.listdir(self.args.fastq)
        fastq_files = [
            f for f in files
            if (f.endswith(".fastq") or f.endswith(".fq") or f.endswith(".fastq.gz") or f.endswith(".fq.gz"))
        ]

        # Figure out how many parallel tasks to launch at once
        # e.g. if self.args.threads=32 and cutadapt uses 8 threads each => 4 parallel jobs
        max_cutadapt_threads = 8
        if self.args.threads <= max_cutadapt_threads:
            n_processes = 1
        else:
            n_processes = self.args.threads // max_cutadapt_threads
            if n_processes < 1:
                n_processes = 1

        print(f"Trimming with up to {n_processes} parallel processes "
              f"(each cutadapt job uses up to {max_cutadapt_threads} threads).")

        # Use a ProcessPoolExecutor for parallel jobs
        futures = []
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            for filename in fastq_files:
                full_path = os.path.join(self.args.fastq, filename)

                # If file is gzipped, decompress it first
                if filename.endswith(".gz"):
                    print(f"Decompressing {filename}...")
                    cmd_unzip = f"gzip -d {full_path}"
                    subprocess.run(cmd_unzip, shell=True, check=True)

                    # e.g. "sample.fastq.gz" => "sample.fastq"
                    full_path = full_path.rsplit(".gz", 1)[0]
                    filename = filename.rsplit(".gz", 1)[0]

                trimmed_output = os.path.join(self.riboSeqTrimmedDir, f"trimmed_{filename}")

                # Submit the job to the executor
                futures.append(
                    executor.submit(
                        self.run_single_cutadapt,
                        full_path,
                        trimmed_output
                    )
                )

            # Optionally wait for each future to complete
            for future in as_completed(futures):
                try:
                    future.result()  # Raises any exception from run_single_cutadapt
                except Exception as e:
                    print(f"Error while trimming: {e}")

        print("All trimming jobs completed.")



    def __check_indexes(self):
        if self.args.aln is None:
            if self.args.genome_index is not None:
                if self.args.genome_index == 'hg38':
                    self.index = f'{self.indexesDir}/hg38.star'
                    self.indexCont = f'{self.indexesDir}/hg38cont.star'
                elif self.args.genome_index == 'hg19':
                    self.index = self.index
                    self.indexCont = self.indexCont
                else:
                    self.index = self.args.genome_index
                    self.indexCont = self.args.cont_index
            else:
                print(f"Genome index not provided. Using default index at {self.index}.")



    def remove_ribosome(self):
        files = os.listdir(self.riboSeqTrimmedDir)
        if self.args.skip_trimming:
            files = os.listdir(self.args.fastq)
            self.riboSeqTrimmedDir = self.args.fastq
        for file in files:
            if file.endswith(".fastq") or file.endswith("fastq.gz") or self.args.skip_trimming:
                out = f'{self.riboSeqContaminantAlnDir}/no_contaminant_{file}Unmapped.out.mate1'
                run = self.verify_checkpoint(outfile=out, step=f"removal of contaminants for {file}.")
                if run:
                    print(f"Removing contaminants from {file}")
                    v = f"{self.toolPaths['STAR']} --version"
                    os.system(v)
                    # print(v)
                    if file.endswith(".gz"):
                        zcat = f' --readFilesCommand zcat'
                    else:
                        zcat = ''
                    cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --outReadsUnmapped Fastx '
                           f'--alignEndsType EndToEnd --genomeDir '
                           f'{self.indexCont} --runThreadN {self.args.threads} --readFilesIn {self.riboSeqTrimmedDir}/{file} '
                           f'--outFileNamePrefix {self.riboSeqContaminantAlnDir}/no_contaminant_{file}{zcat}')
                    os.system(cmd)
                    print(f"Finished removing contaminants from {file}")

    def align_rpfs(self):
        files = os.listdir(self.riboSeqContaminantAlnDir)
        for file in files:
            if 'Unmapped.out.mate1' in file or file.endswith(".fastq"):
                out = f'{self.riboSeqAlnDir}/aligned_to_genome_{file}Aligned.out.sam'
                run = self.verify_checkpoint(outfile=out, step=f"alignment of clean Ribo-Seq reads to the genome")
                if run:
                    if file.endswith(".gz"):
                        gz = ' --readFilesCommand zcat'
                    else:
                        gz = ''
                    filepath = f'{self.riboSeqContaminantAlnDir}/{file}'
                    # filepath = f'{self.riboSeqTrimmedDir}/{file}'
                    print(f"Aligning RPF reads from {file} to {self.index}")
                    # cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --genomeDir {self.index} --runThreadN '
                    #        f'{self.args.threads} --readFilesIn {filepath} --outFileNamePrefix '
                    #        f'{self.riboSeqAlnDir}/aligned_to_genome_{file} --outFilterMismatchNmax 2 '
                    #        f'--outFilterMultimapNmax {self.args.multimappings} --chimScoreSeparation 10 --chimScoreMin '
                    #        f'20 --chimSegmentMin 15 --outSAMattributes All{gz}')
                    clipbases = ''
                    if self.args.clip5pNbases is not None:
                        clipbases = f' {self.args.clip5pNbases}'
                    cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --genomeDir {self.index} --runThreadN '
                           f'{self.args.threads} --readFilesIn {filepath} --outFileNamePrefix '
                           f'{self.riboSeqAlnDir}/aligned_to_genome_{file} '
                           f'--outFilterMultimapNmax {self.args.multimappings} '
                           f'--outSAMattributes All{gz}{clipbases}')
                    os.system(cmd)
                    print(f'Finished aligning {file}.')



