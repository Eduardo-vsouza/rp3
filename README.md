# rp3
Ribosome Profiling and Proteogenomics Pipeline (RP3) for the identification of novel microproteins encoded by smORFs

Table of contents
1. [Introduction[(## 1. Introduction)
2. [Database generation[](## 5. Database)

## 1. Introduction
RP3 (Ribosome Profiling and Proteogenomics pipeline) was developed to integrate the analyses for three different multi-omics techniques: RNA-Seq, Ribo-Seq and Proteogenomics. Its overarching goal is to identify novel microproteins (shorter than 100/150 aa) encoded by small Open Reading Frames (smORFs). Then, it will check for translational evidence in the Ribo-Seq data for these novel smORFs.

  
#### Ribo-Seq 
Ribosome Profiling, or Ribo-Seq, is a widely used technique for the identification of novel coding genes, as it provides direct evidence of translation. It consists of the sequencing of fragments of mRNA that are being actively read by the ribosome (RPFs). Pipelines usually map these RPFs to transcripts from a de novo or reference-guided transcriptome assembly to check for the existence of novel transcripts that might carry smORFs (or any gene) actively being translated by the ribosome. Here, the Ribo-Seq reads will be used to check for translational evidence of novel microprotein-encoding smORFs identified by the Proteogenomics part of the pipeline (Fig. 1 steps VI and VII). 
  
![image](https://github.com/Eduardo-vsouza/rp3/assets/60533781/8c689efa-e7d4-4501-92ce-f604572e82ac)


  
*Disclaimer*: Note that this figure differs from the one in the paper published at *Journal*. In the paper, the workflow in this figure was used to explain the sequence of analyses that were performed. The RP3 pipeline does not cover the usage of RibORF for smORF identification. Instead, it covers the steps necessary to identify microproteins with Proteogenomics. Additionally, it checks for translational evidence for PG smORFs by making use of the Ribo-Seq reads. It does so by mapping the Ribo-Seq reads back to the PG smORFs coordinates in the genome.

## 2. Installation
Typical installation time will vary depending on how many dependencies requirements are already met. Downloading the release from the GitHub page should not take more than 5 min. at average connection speed, and installing all dependencies should take no more than 20 min.

1. Download the latest version of the pipeline from its GitHub releases page.
2. ``cd`` to the directory containing the .tar.gz file and run ``tar -xzvf`` on the terminal.
3. Install required Python packages
	1. ``$ cd`` to the RP3 folder
	2. ``$ pip install -r requirements.txt
4. Install dependencies
	To avoid incompatibility, make sure to have installed the versions for each tool that was used during the development of RP3. You can install alternative versions at your own risk. Compatibility is not guaranteed.
	-  Obtain a license and install MSFragger v3.5 (https://msfragger.nesvilab.org/)
	-  Install Percolator v3.06.1 (percolator.ms)
	-  Install NCBI Blast v2.12.0 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-src.tar.gz)
	-  Install STAR v2.7.4a (https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz)
	-  Install StringTie v2.2.1
	-  Install the subread package v1.6.3. It contains the tool featureCounts, used for the ``ribocov`` mode.
	-  Install samtools
	-  Install MSBooster (only necessary if using the flag --msBooster to predict retention times.)
	
5.  **Important**. Configure the Paths to each of the dependencies in the file config.txt located inside the RP3 folder. Replace the $PATH to each tool in its respective column. By default, the pipeline will look for the tools in your $PATH. 
#### Testing the installation
We provide a *demo* mode with reduced datasets so the user can check if the installation is working properly. This mode will check the 6 main modes (*rna, translation, database, search, postms,* and *ribocov*). To check all modes at once, simply run 
`$ rp3.py demo --threads 8`
This will use 8 threads to test all 6 main modes of the RP3 pipeline. Typical run time for this is 20-30 min, but can vary depending on available computational resources. 
The output files will be generated at `demo_outdir` inside the pipeline's installation directory.
##### Test dataset 
The test data is composed of single files (to enable fast testing of the software's core functionalities) from studies used in the published manuscript. For each mode:
``rna``: GEO Accession SRR8449601_1.fastq and SRR8449601_2.fastq paired end fastq files from RNA-Seq experiments from GEO Series (GSE198109)
``search``: the mzML file 20130328_EXQ1_MiBa_SA_HCC1937.mzML from MassIVE (accession MSV000089022).
``ribocov``: SRR8449580.fastq file containing Ribo-Seq reads from GEO Series (GSE198109). 
Additionally, reference annotation files are included in the testing datasets. These are used for nearly every mode of the pipeline:
Reference GTF, rRNA and tRNA fasta, and genome Fasta files from hg19 assembly from UCSC. Human RefSeq from latest assembly from NCBI. Human reference proteome from Uniprot.

Every file is located inside the `demo_data` directory, located inside the Rp3 directory.
## 3. Quick start

1. ##### I already have an assembled transcriptome
If you already have a transcriptome, start the pipeline on the ``database`` mode and provide the path to the folder containing your GTF files with the ``--gtf_folder`` parameter (skip to step 3)

2.  ##### I have RNA-Seq reads, but still need to assemble a transcriptome
Run the RP3 pipeline on ``rna`` mode.
$ ``rp3.py rna --outdir /path/to/output/directory --threads 8 --strandness rf --gtf /path/to/gtf/file --lib_type single --reads_folder /path/to/folder/containing/ribo-seq/fastq/files``

3. ##### Generate a custom proteomics database
Run RP3 on the ``database`` mode
``$ rp3.py database --outdir <path/to/output/directory> --threads 8 --genome <path/to/genome.fasta> --gtf_folder <path/to/gtf/folder> --proteome <path/to/reference_proteome.fasta``

- The specified ``--gtf_folder`` should be the folder containing the GTF file with the transcriptome assembly. If you have assembled it with the ``rna`` mode, it is located inside the output directory at ``/transcriptomics/assembly``. If you already have a GTF file, provide the folder containing it and other GTF files to be included in the analysis. **Do not provide the path directly to the GTF file. The pipeline expects it to be inside a directory, optionally with other GTF files.**

4. ##### Perform the peptide search
Run RP3 on the ``search`` mode.
``$ rp3.py search --mzml /path/to/mass/spec/folder --outdir path/to/output/directory/ --threads 8`` 
	- The ``--mzml`` expects a folder containing subdirectories corresponding to each group (see Peptide search section)

5. ##### Re-score your results
Run the pipeline on ``rescore`` mode.
``$ rp3.py rescore --outdir /path/to/output/directory --threads 8 --mzml /path/to/mzmz/files --proteome /path/to/reference/proteome --msPattern mzML``


6. ##### Check if the proteogenomics smORFs have Ribo-Seq coverage
Run the pipeline on ``ribocov`` mode.
``$ rp3.py ribocov --outdir /path/to/output/directory --threads 8 --fastq /path/to/fastq/folder --gtf /path/to/gtf/file --genome_index /path/to/genome/index --cont_index /path/to/contaminants/index --plots``

##### Notes 
- Always use the same output directory when running different modes for the same analysis.

## 4. Transcriptome assembly

1. Run the transcriptome assembly mode of the pipeline

	This mode focuses on the generation of a transcriptome assembly from RNA-Seq reads. This will be used during the next step to generate a custom mass spectrometry database containing all the possible (or most likely) smORFs in the assembled transcriptome.
	
	**Note**: this step can be skipped if you already have a transcriptome. All you need to provide is the GTF files in the next step. Alternatively, you can use a reference GTF file to identify novel smORFs in already annotated transcripts.
	
2. Parameters
	Run ``$ rp3.py rna -h`` at the terminal to print this message:
```
General Parameters:
  rna
  --outdir OUTDIR, -o OUTDIR
                        Inform the output directory (default: foopipe)
  --threads THREADS, -p THREADS
                        Number of threads to be used. (default: 1)

rna options:
  --strandness STRANDNESS
  --gtf GTF
  --lib_type LIB_TYPE   Single or Paired. (default: None)
  --reads_folder READS_FOLDER
                        Folder containing sequencing reads in fastq format. (default: None)
  genome_index          The genome index generated for hisat2 with the --ss and --exon
                        arguments. If not provided, a new index will be generated. This might
                        take a while.

```

**Note**. The ``--reads_folder`` flag expects that reads are organized by groups, such as:
	
```
reads_folder/
	├── group_1/
	│   └── reads_1.fastq
	├── group_2/
	│   └── reads_2.fastq
	└── group_3/
	    └── reads_3.fastq
 ```
In case you have a single group and thus only one transcriptome to assemble, put all the reads in the same group folder inside the reads_folder directory, like:
```
reads_folder/
	├── group_1/
	   └── reads_1.fastq
	   └── reads_2.fastq
	   └── reads_3.fastq
```

## 5. Database

1. Run the RP3 pipeline in its ``database`` mode to generate a custom database for mass spectrometry-based proteomics. Run the pipeline with ``$ rp3.py database -h`` to print all available commands for this mode.
 
2. Run the database mode
	 - Alternatively, you can run the ``translation`` mode separately, in case you do not need the reference proteome, decoy sequences and contaminants added to the fasta file. You can always run the ``database`` mode with the ``--skip_translation`` flag later on to generate a database suited for peptide search for mass spectrometry data.
	 - By default, ``database`` automatically executes the ``translation`` mode.
	 
3. Parameters
	Run ``$ rp3.py database -h`` on the command line to check which commands are available. This will print out this message:
		
```
General Parameters:
  database
  --outdir OUTDIR, -o OUTDIR
						Inform the output directory (default: foopipe)
  --threads THREADS, -p THREADS
						Number of threads to be used. (default: 1)

database options:
  --proteome PROTEOME   Reference proteome (default: None)
  --genome GENOME
  --gtf_folder GTF_FOLDER
  --external_database EXTERNAL_DATABASE
  --skip_translation
  ```
		    
		  
4. Example code
	4.1. Generating a custom database for proteomics from a transcriptome assembly or custom GTF files
	
	``$ rp3.py database --outdir <path/to/output/directory> --threads 8 --genome <path/to/genome.fasta> --gtf_folder <path/to/gtf/folder> --proteome <path/to/reference_proteome.fasta``
			
	This will tell the pipeline to generate databases in the ``databases`` directory inside the output directory provided with ``--outdir``. 
	
	**Remember to always specify the same output directory in the other modes, as each step is dependent on the previous one. **
	
	4.2. In case you already have your own database with putative sequences, you can provide it to the pipeline with the ``--external_database`` flag along with the reference proteome. In that case, the pipeline will append the reference proteome to your fasta file along with the decoys and contaminants. It will then tag each sequence for the subsequent steps.
			
	2.3. Notes
		- ``--gtf_folder`` expects a folder with one or more GTF files.

5. Output files
	Database files will be generated inside the output directory in the ``databases`` folder. You can check for database metrics at ``metrics/database_metrics.txt``.
	
			
## 6. Peptide search

1. During this step, the pipeline will use the database generated in the ``database`` to look for evidence of microproteins in the mass spectrometry data. 
2. Make sure to have the data stored in the centroided ``.mzML`` format. In case it isn't, use the tool msconvert from the ProteoWizard suite (https://proteowizard.sourceforge.io/).
3. By default, the ``search`` mode will run the ``postms`` mode as well, which employs Percolator to assess the FDR and performs the necessary cutoffs. You can run the ``search`` mode by itself by specifying the flag ``--no_post_process``. In that case, only MSFragger will be run and you will have ``.pin`` files generated from the search. You can always run the ``postms`` mode later on to process the search files.
4. Requirements
	1. Raw data from label-free LC-MS/MS experiments (``mzML``).
	2. Database generated by the ``database`` mode OR provided with the flag ``--external_database``. 
5. Parameters
	Run ``$ rp3.py search -h`` at the terminal to print this message:
```
General Parameters:
  search
  --outdir OUTDIR, -o OUTDIR
                        Inform the output directory (default: foopipe)
  --threads THREADS, -p THREADS
                        Number of threads to be used. (default: 1)

search options:
  --mzml MZML
  --digest_max_length DIGEST_MAX_LENGTH
  --digest_min_length DIGEST_MIN_LENGTH
  --std_proteomics
  --quantify
  --mod MOD
  --create_gtf
  --cat                 Perform the search using a concatenated target and decoy database.
                        Requires the databases to be generated using the 'cat' flag.
                        (default: False)
  --tmt_mod TMT_MOD
  --fragment_mass_tolerance FRAGMENT_MASS_TOLERANCE
  --refseq REFSEQ
  --groups GROUPS       Tab-delimited file associating a database with a raw file. Should
                        contain two columns: files, groups. Groups should have the same name
                        as the generated databases. If not specified, the pipeline will
                        search every raw file using every GTF file provided. (default: None)

```

Just like in the ``rna`` mode, the ``-mzml`` flag expects the mzml folder to contain groups as well, such as:
```
mass_spec_folder/
	├── group_1/
	│   └── lc_ms-ms_1.mzML
	├── group_2/
	│   └── lc_ms-ms_2.mzML
	└── group_3/
	    └── lc_ms-ms_3.mzML
	    └── lc_ms-ms_4.mzML
```

In case you have a single group/condition, put all the mzML files inside a folder in the ``-mzml`` directory, such as:
```
mass_spec_folder/
	├── group_1/
	    └── lc_ms-ms_1.mzML
	    └── lc_ms-ms_2.mzML
	    └── lc_ms-ms_3.mzML
	    └── lc_ms-ms_4.mzML
```

In that case, specify the ''mass_spec_folder'' for the ``--mzml``  parameter, and **not** the group folder.


###### Notes
- The ``--groups`` parameter allows you to specify which GTF files should be used for each mzML group. This is useful in case you have a transcriptome for condition X, Y and mass spec groups for the conditions X and Y. In that case, you need to specify, for each mass spec file, the GTF group that should be used with it. Note that it should be the same name as the groups provided in the ``-gtf_folder`` parameter for the ``database`` and/or ``translation`` modes. **In case the ``--groups`` file is not provided, all databases generated from the provided GTF files will be used to search each mzML file.**
	The groups.txt file should be organized as a tab-separated table, like:
	files    groups
	file_1.mzML    Y
	file_2.mzML    X
- Always specify the same ``--outdir`` previously used for the other modes. 
- The ``--refseq`` parameter will accept a fasta file containing a reference annotation, such as the one from NCBI RefSeq. This is used as an extra sanity check to make sure we remove all annotated microproteins, even those that might be missed by the pipeline in case the provided reference proteome (from Uniprot, for instance) is not comprehensive enough. This parameter is optional, but recommended.
- 


6. Example code
	``$ rp3.py search --mzml /path/to/mass/spec/folder --outdir path/to/output/directory/ --threads 8 

7. Output files
	RP3 will produce output files in fasta format for each of the provided groups. Look for them inside the output directory at ``summarized_results/group_name``. Merged files from all the groups are located inside ``summarized_results/merged``.

## 7. Post-processing

1. The RP3 pipeline contains a re-scoring mode called ``rescore``. This is intended to perform a second round of searches, now using as a proteomics database the results from the first proteogenomics search (the fasta file generated by the ``search`` mode) appended to the reference proteome. This is useful because the FDR assessment from the first search is not very accurate, as the database generated from the three-frame translation of the transcriptome contains millions of predicted sequences. This bloated database results in false positives and false negatives during FDR assessment. To correct for this, we select the hits at an FDR < 0.01 from the first search and look for them again, now with a smaller database to obtain more accurate hits. This mode will reduce the final number of novel microproteins.
2. After running the ``search`` mode, run the ``rescore`` in the same output directory:
	``$ rp3.py rescore --outdir /path/to/output/directory --threads 8 --mzml /path/to/mzmz/files --proteome /path/to/reference/proteome --msPattern mzML
###### Notes
- The ``--msPattern`` specifies the format of the files (usually mzML or bruker (.d) format).

3. Output files
	Look for output files in fasta and gtf format in the ``rescore/summarized_results`` directory inside the output directory. 

## 8. Ribocov
This mode will check for Ribo-Seq coverage for the microproteins identified with proteogenomics. To do so, it will run featureCounts on a custom GTF file automatically generated by the pipeline. The available parameters are:
```
General Parameters:
  ribocov
  --outdir OUTDIR, -o OUTDIR
                        Inform the output directory (default: None)
  --threads THREADS, -p THREADS
                        Number of threads to be used. (default: 1)

ribocov options:
  --fastq FASTQ         Provide the path to the folder containing fastq files
                        to be aligned to the genome. If the --aln argument is
                        provided, this is not necessary. (default: None)
  --gtf GTF             Reference gtf file containing coordinates for
                        annotated genes. The novel smORFs sequences from the
                        proteogenomics analysis will be appended to it.
                        (default: None)
  --genome_index GENOME_INDEX
                        Path to the genome STAR index. If not provided, it
                        will use the human hg19 index available at /data/
                        (default: None)
  --cont_index CONT_INDEX
                        STAR index containing the contaminants (tRNA/rRNA
                        sequences). Reads mapped to these will be excluded
                        from the analysis. (default: None)
  --aln ALN             Folder containing bam or sam files with Ribo-Seq reads
                        aligned to the genome. In case this is provided,
                        indexes are not required and the alignment step will
                        be skipped. (default: None)
  --rpkm RPKM           RPKM cutoff to consider whether a smORF is
                        sufficiently covered by RPFs or not. (default: 1)
  --multimappings MULTIMAPPINGS
                        max number of multimappings to be allowed. (default:
                        99)
  --adapter ADAPTER     Provide the adapter sequence to be removed. (default:
                        AGATCGGAAGAGCACACGTCT)
  --plots
  --fastx_clipper_path FASTX_CLIPPER_PATH
  --fastx_trimmer_path FASTX_TRIMMER_PATH

```

To run the RP3 pipeline on ribocov mode, run:
``rp3.py --outdir path/to/output/directory --threads 8 --gtf path/to/gtf/file --fastq path/to/fastq/folder``
This will use the provided genome indexes for the human hg19 assembly located inside the STAR_indexes directory, located inside the rp3 main directory. The user can also generate new indexes if they require to do so. In that case, provide the path to them using the parameters ``--genome_index`` and ``cont_index``. Make sure to change the ``--adapter`` parameter to suit the adapter sequence used for your Ribo-Seq experiment.
The output files will be located inside the ``counts`` directory. They will include a heatmap showing the overall Ribo-Seq coverage for the proteogenomics smORFs, as well as a table containing information about the mapping groups. If the ``--plots`` argument was specified, a plot showing the number of Ribo-Seq-covered smORFs in each mapping group will be generated at ``counts/plots``. 
	

  
  
  
