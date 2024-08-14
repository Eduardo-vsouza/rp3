Rp3
=============================
<img src="https://github.com/Eduardo-vsouza/rp3/assets/60533781/ee5208ad-5a1a-42a0-830e-13fbfc8cef20" width="230" height="150" /> \
``Design by Manami Nishikawa``

Rp3 (Ribosome Profiling and Proteogenomics pipeline) was developed to integrate the analyses for three different multi-omics techniques: RNA-Seq, Ribo-Seq and Proteogenomics. Its overarching goal is to identify novel microproteins (shorter than 100/150 aa) encoded by small Open Reading Frames (smORFs). Then, it will check for translational evidence in the Ribo-Seq data for these novel smORFs. 
  
### Ribo-Seq 
Ribosome Profiling, or Ribo-Seq, is a widely used technique for the identification of novel coding genes, as it provides direct evidence of translation. It consists of the sequencing of fragments of mRNA that are being actively read by the ribosome (RPFs). Pipelines usually map these RPFs to transcripts from a de novo or reference-guided transcriptome assembly to check for the existence of novel transcripts that might carry smORFs (or any gene) actively being translated by the ribosome. Here, the Ribo-Seq reads will be used to check for translational evidence of novel microprotein-encoding smORFs identified by the Proteogenomics part of the pipeline (Fig. 1 steps VI and VII). 
  
![image](https://github.com/Eduardo-vsouza/rp3/assets/60533781/8c689efa-e7d4-4501-92ce-f604572e82ac)

 
*Disclaimer*: Note that this figure differs from the one in the paper published in *Nature Communications*. In the paper, the workflow in this figure was used to explain the sequence of analyses that were performed. The RP3 pipeline does not cover the usage of RibORF for smORF identification. Instead, it covers the steps necessary to identify microproteins with Proteogenomics. Additionally, it checks for translational evidence for PG smORFs by making use of the Ribo-Seq reads. It does so by mapping the Ribo-Seq reads back to the PG smORFs coordinates in the genome.

## System requirements
The pipeline was tested on Ubuntu 22.04 LTS. 
A typical desktop with 8 cores and 16 GB RAM should be enough to run most of the analyses, although if running the workflow with these specs, it might take a long time to perform every step.


## Installation
Typical installation time will vary depending on how many dependencies requirements are already met. Downloading the release from the GitHub page should not take more than 5 min. at average connection speed, and installing all dependencies should take no more than 20 min.

1. Download the latest version of the pipeline from its GitHub releases page.
2. ``cd`` to the directory containing the .tar.gz file and run ``tar -xzvf`` on the terminal.
3. Install required Python packages
  	1. Install Python3.10. Although some backward compatibility is expected, other versions of Python haven't been tested.
	2. ``$ cd`` to the RP3 folder
	3. ``$ pip install -r requirements.txt
5. Install dependencies
	To avoid incompatibility, make sure to have installed the versions for each tool that was used during the development of RP3. You can install alternative versions at your own risk. Compatibility is not guaranteed.

| Tool          | Mode                                     | Version          | Purpose   | Link                                                                                                   |
|--------------|------------------------------------------|------------------|-----------|--------------------------------------------------------------------------------------------------------|
| MSFragger    | search, rescore                          | 3.5              | Match MS spectra against protein database   | [Download](https://msfragger.nesvilab.org/)                                                            |
| Percolator    | search, postms, rescore                  | 3.06.1           | Post-process MSFragger results and infer FDR | [Download](percolator.ms)                                                                              |
| NCBI Blast	| search (--refseq), anno (--conservation) | 2.12.0	      | Local alignment of microproteins | [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-src.tar.gz) |
| STAR			| 	        ribocov, rna                    | 2.7.4a	      | Alignment of short-reads to the genome in a splice-aware manner | [Download](https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz)                                    |
| StringTie		| 	             rna                        | 2.2.1 	      | Reference-guided transcriptome assembly | [Download](https://github.com/gpertea/stringtie/releases/tag/v2.2.1)                                   |
| Subread		| 	                  ribocov               | 1.6.3	      | Contains featureCounts, used to perform read counting for Ribo-Seq data in the ``ribocov`` mode | [Download](https://sourceforge.net/projects/subread/files/subread-1.6.3/)                              |
| samtools		| 	                    ribocov, rna        | 1.18	      | Performs a variety of data processing for ``.sam`` and ``.bam`` files | [Download](https://github.com/samtools/samtools/releases/tag/1.18)                                     |
| MSBooster		| 	                  search, rescore       | 1.2.1	      | Predicts Retention Times and add delta_RT_loess to ``.pin`` files coming from MSFragger | [Download](https://github.com/Nesvilab/MSBooster/releases/tag/v1.2.1)                                  |

	
6.  **Important**. Configure the Paths to each of the dependencies in the file config.txt located inside the RP3 folder. Replace the $PATH to each tool **executable** in its respective column. By default, the pipeline will look for the tools in your $PATH.

### Testing the installation
We provide a *demo* mode with reduced datasets so the user can check if the installation is working properly. This mode will check the 5 main modes (*translation, database, search, postms,* and *ribocov*). 
- First, you need to download the demo_data from the releases page and put it in the same directory as the Rp3.py script.
- Then, to check all modes at once, simply run 
`$ rp3.py demo --threads 8 --outdir demo_outdir`
This will use 8 threads to test all 5 main modes of the RP3 pipeline. Typical run time to test every mode is ~30-50 min, but can vary depending on available computational resources. The whole workflow of the pipeline is time-consuming because it has to deal with multiple types of omics datasets. 
The output files will be generated at `demo_outdir`, or another specified directory.
- If you want to skip a mode during testing, pass the argument --skip_database, for instance. All parameters available for the ``demo`` mode can be checked with ``rp3.py demo -h``.

- If testing also the Ribocov mode, you will need to provide STAR indexes for the hg19 genome. These are very large files and we make them available for download in this box folder: https://salkinstitute.box.com/s/5uyf0tdfm7w1zgx3kqr1u367zlg9ocoy. Alternatively, you can use your own indexes. For demo mode, make sure to put these indexes inside the sofware directory under the folder ``STAR_indexes``.
###### Test dataset 
The test data is composed of single files (to enable fast testing of the software's core functionalities) from studies used in the published manuscript. For each mode:
``search``: the mzML file 20130328_EXQ1_MiBa_SA_HCC1937.mzML from MassIVE (accession MSV000089022).
``ribocov``: SRR8449580.fastq file containing Ribo-Seq reads from GEO Series (GSE198109). 
Additionally, reference annotation files are included in the testing datasets. These are used for nearly every mode of the pipeline:
Reference GTF, rRNA and tRNA fasta, and genome Fasta files from hg19 assembly from UCSC. Human RefSeq from latest assembly from NCBI. Human reference proteome from Uniprot.

Every file is located inside the `demo_data` directory, located inside the Rp3 directory. If you want to download the demo_data, get it from the release page, separately from the source code.

## Quick start

1. ###### I already have an assembled transcriptome
If you already have a transcriptome, start the pipeline on the ``database`` mode and provide the path to the folder containing your GTF files with the ``--gtf_folder`` parameter (skip to step 3)

2.  ###### I have RNA-Seq reads, but still need to assemble a transcriptome
You should assemble your transcriptome using a transcriptome assembler such as StringTie or Cufflinks after aligning the reads to the genome with a splice-aware aligner such as STAR or hisat2. Rp3 needs the GTF file produced by these assemblers.

###### Using the reference GTF
Alternatively, you can just use a reference GTF file from Ensembl or NCBI, for instance. The GTF files has to contain 'transcript' and 'exon' features. The limitation when using a reference GTF is that you cannot include novel transcripts for your samples. The pipeline will then perform the 3-frame translation on already annotated transcripts. Annotated ORFs will be removed afterwards.

3. ###### Generate a custom proteomics database
Run RP3 on the ``database`` mode
``$ rp3.py database --outdir <path/to/output/directory> --threads 8 --genome <path/to/genome.fasta> --gtf_folder <path/to/gtf/folder> --proteome <path/to/reference_proteome.fasta``

- The specified ``--gtf_folder`` should be the folder containing the GTF file with the transcriptome assembly. If you already have a GTF file, provide the folder containing it and other GTF files to be included in the analysis. **Provide the path to a GTF file.**

4. ###### Perform the peptide search
Run RP3 on the ``search`` mode.
``$ rp3.py search --mzml /path/to/mass/spec/folder --outdir path/to/output/directory/ --threads 8`` 
	- The ``--mzml`` expects a folder containing subdirectories corresponding to each group (see Peptide search section)

5. ###### Re-score your results
Run the pipeline on ``rescore`` mode.
``$ rp3.py rescore --outdir /path/to/output/directory --threads 8 --mzml /path/to/mzmz/files --proteome /path/to/reference/proteome --msPattern mzML``


6. ###### Check if the proteogenomics smORFs have Ribo-Seq coverage
Run the pipeline on ``ribocov`` mode.
``$ rp3.py ribocov --outdir /path/to/output/directory --threads 8 --fastq /path/to/fastq/folder --gtf /path/to/gtf/file --genome_index /path/to/genome/index --cont_index /path/to/contaminants/index --plots``

###### Notes 
- Always use the **same output directory** when running different modes for the same analysis. 


## Database

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
	
	4.2. In case you already have your own database with putative sequences, you can provide it to the pipeline with the ``--external_database`` flag along with the reference proteome. In that case, the pipeline will append the reference proteome to your fasta file along with the decoys and contaminants. It will then tag each sequence for the subsequent steps. This is useful when not translating a transcriptome assembly to the 3-reading frames, like when searching for mass spectrometry evidence for microproteins identified with Ribo-Seq. 
			
	2.3. Notes
		- ``--gtf_folder`` expects a folder with one or more GTF files or the path to a GTF file.

5. Output files
	Database files will be generated inside the output directory in the ``databases`` folder. You can check for database metrics at ``metrics/database_metrics.txt``.
### Notes
Unless specified in the next modes, databases will have annotated proteins excluded based on the ``-proteome`` provided. This behavior is intended to remove, by default, any already annotated microprotein and peptides that match to them. If you are interested in annotated proteins as well, there are two ways to prevent this from happening:
- Specify the parameter ``--uniprotAnnotation`` during the ``database`` mode. This will use a data frame obtainable from Uniprot containing the annotation level for each protein in the ``--proteome``. Rp3 will then tag differently each protein in the provided proteome based on their annotation level. As proteins in Uniprot vary in their annotation, it might be wise to include those with very low annotation levels when identifying unannotated proteins, as they are poorly characterized. The annotation level ranges from 1 to 5. To control which should be kept, provide the argument ``--annotationLevel`` followed by a number from 1 to 5. ``--annotationLevel 4`` will keep only those proteins with an annotation level equal to or lower than 4. This is the default behavior if ``--uniprotAnnotation`` is provided. **Requires ``--includeLowAnnotation`` to be specified during the ``search`` mode.**
- Specify ``--keepAnnotated`` in ``search`` mode, the next step. This will treat every protein, annotated or unannotated, the same.
			
## Peptide search

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

The ``--mzml`` flag expects the mzml folder to contain groups, such as:
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

In case you have a single group/condition, put all the mzML files inside a folder in the ``--mzml`` directory, such as:
```
mass_spec_folder/
	├── group_1/
	    └── lc_ms-ms_1.mzML
	    └── lc_ms-ms_2.mzML
	    └── lc_ms-ms_3.mzML
	    └── lc_ms-ms_4.mzML
```

In that case, specify the ''mass_spec_folder'' for the ``--mzml``  parameter, and **not** the group folder.


### Notes
- The ``--groups`` parameter allows you to specify which GTF files should be used for each mzML group. This is useful in case you have a transcriptome for condition X, Y and mass spec groups for the conditions X and Y. In that case, you need to specify, for each mass spec file, the GTF group that should be used with it. Note that it should be the same name as the groups provided in the ``-gtf_folder`` parameter for the ``database`` and/or ``translation`` modes. **In case the ``--groups`` file is not provided, all databases generated from the provided GTF files will be used to search each mzML file.**
	The groups.txt file should be organized as a tab-separated table, like:
	files    groups
	file_1.mzML    Y
	file_2.mzML    X
- Always specify the same ``--outdir`` previously used for the other modes.
  ##### Extra parameters
- The ``--refseq`` parameter will accept a fasta file containing a reference annotation, such as the one from NCBI RefSeq. This is used as an extra sanity check to make sure we remove all annotated microproteins, even those that might be missed by the pipeline in case the provided reference proteome (from Uniprot, for instance) is not comprehensive enough. This parameter is optional, but recommended.
- You can also include the execution of  ``rescore`` mode following the search and post-processing with Percolator. To do so, specify the flag ``--rescore`` and provide the path to the proteome fasta file with ``--proteome``. The proteome should be the same one used to generate the database in the first step. The proteome is required if rescoring the results, but not for the first round of searches. For details, see [Post-processing](#post-processing)
- The flag ``--MSBooster`` will generate a spectral library with predicted retention times (RT) and delta RT loess compared to the experimental data. These values will be incorporated in the ``.pin`` file used as input for post-processing with Percolator. This can either increase or reduce the number of identifications depending on the analysis, but should improve confidence.


6. Example code
	``$ rp3.py search --mzml /path/to/mass/spec/folder --outdir path/to/output/directory/ --threads 8 --MSBooster --rescore --proteome path/to/proteome.fasta``

7. Output files
	RP3 will produce output files in fasta format for each of the provided groups. Look for them inside the output directory at ``summarized_results/group_name``. Merged files from all the groups are located inside ``summarized_results/merged``.
 
## Post-processing

1. The RP3 pipeline contains a re-scoring mode called ``rescore``. This is intended to perform a second round of searches, now using as a proteomics database the results from the first proteogenomics search (the fasta file generated by the ``search`` mode) appended to the reference proteome. This is useful because the FDR assessment from the first search is not very accurate, as the database generated from the three-frame translation of the transcriptome contains millions of predicted sequences. This bloated database results in false positives and false negatives during FDR assessment. To correct for this, we select the hits at an FDR < 0.01 from the first search and look for them again, now with a smaller database to obtain more accurate hits. This mode will reduce the final number of novel microproteins.
2. After running the ``search`` mode, run the ``rescore`` in the same output directory:
	``$ rp3.py rescore --outdir /path/to/output/directory --threads 8 --mzml /path/to/mzmz/files --proteome /path/to/reference/proteome --msPattern mzML
####### Notes
- The ``--msPattern`` specifies the format of the files (usually mzML or bruker (.d) format).

3. Output files
	Look for output files in fasta and gtf format in the ``rescore/summarized_results`` directory inside the output directory. 

## Ribocov
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

## Annotation mode
Rp3 provides an additional mode, ``anno``, to provide additional information on the identified microproteins. Running 
Rp3 with ``rp3.py anno --help`` will return:
```commandline
 ____       _____
|  _ \ _ __|___ /
| |_) | '_ \ |_ \
|  _ <| |_) |__) |
|_| \_\ .__/____/
      |_|  
RP3 v1.1.0
usage: /home/microway/PycharmProjects/rp3/rp3.py anno [-h] [--outdir OUTDIR]
                                                      [--threads THREADS]
                                                      [--overwrite]
                                                      [--signalP]
                                                      [--organism ORGANISM]
                                                      [--conservation]
                                                      [--blast_db BLAST_DB]
                                                      [--rescored]
                                                      [--uniprotTable UNIPROTTABLE]
                                                      [--orfClass]
                                                      [--paralogy] [--mhc]
                                                      [--repeats] [--isoforms]
                                                      [--exclusiveMappingGroups]
                                                      [--affinity AFFINITY]
                                                      [--affinityPercentile AFFINITYPERCENTILE]
                                                      [--filterPipeResults]
                                                      [--genome GENOME]
                                                      [--alignToTranscriptome]
                                                      [--maxMismatches MAXMISMATCHES]
                                                      [--gtf GTF]
                                                      [--repeatsFile REPEATSFILE]
                                                      [--refGTF REFGTF]
                                                      anno

Run pipeline_config in anno mode

options:
  -h, --help            show this help message and exit

General Parameters:
  anno
  --outdir OUTDIR, -o OUTDIR
                        Inform the output directory (default: None)
  --threads THREADS, -p THREADS
                        Number of threads to be used. (default: 1)
  --overwrite

anno options:
  --signalP
  --organism ORGANISM
  --conservation
  --blast_db BLAST_DB
  --rescored            Use this flag if the 'rescore' mode was used to
                        perform a second round of search using the results
                        from the first search. Only the rescored microproteins
                        will be analyzed for conservation in this case.
                        (default: False)
  --uniprotTable UNIPROTTABLE
  --orfClass
  --paralogy
  --mhc
  --repeats
  --isoforms
  --exclusiveMappingGroups

MHC detection parameters.:
  --affinity AFFINITY
  --affinityPercentile AFFINITYPERCENTILE
  --filterPipeResults

Paralogy parameters.:
  --genome GENOME
  --alignToTranscriptome
  --maxMismatches MAXMISMATCHES

ORF Classification parameters.:
  --gtf GTF             reference GTF file. For better accuracy in annotation,
                        this should be a GTF file from Ensembl. They contain
                        more terms that help better classifying the smORF.
                        (default: None)

Repeats parameters.:
  --repeatsFile REPEATSFILE

Isoforms parameters.:
  --refGTF REFGTF
```

With this mode, it's possible to identify signal peptides (running signalP6.0), conservation, orf classes, presence of 
MHC epitopes, and presence of paralogs in the genome. To identify signal peptides and annotate smORF classes, run Rp3 
as:

```$ rp3.py anno --outdir /path/to/outdir/from/previous/modes --signalP --orfClass --gtf /path/to/ensembl/gtf```

To define smORF classes in the manuscript, we used the annotation from Ensembl, which we believe to be very 
comprehensive and allows us to get better insight into our data. To obtain a GTF file from the human genome 
assembly hg19, for instance, go to: https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/ and download the 
appropriate GTF file.


## Ways to run Rp3
 Rp3 allows you the perform proteogenomics analysis in a few different ways. 
 Every workflow will require:
 - Mass Spectrometry data (LC-MS/MS, DDA) in ``.mzML``
 - Reference proteome in ``.fasta``
 - Output directory ``--outdir``. This has to be the same for every mode for a given analysis.
### Identify unannotated microproteins with proteomics and translational evidence from a transcriptome
This will use the proteogenomics approach to identify any possible microprotein in the provided transcriptome. Then, Rp3 will check for translational evidence for these microproteins using the Ribo-Seq reads.
#### Additional requirements: 
- Ribo-Seq reads in ``.fastq`` 
- GTF file (either a *de novo* assembled or a reference transcriptome) OR - Fasta file containing a database of interest
- Genome ``.fasta`` file. Must be the same genome assembly as the GTF file, as the three-frame translation of the GTF will be based on this ``.fasta``.
#### Workflow
- Run the ``database`` mode providing a ``--gtf_folder``, ``--proteome``, ``--genome``, and ``--cat``.
- Run ``search`` mode providing ``--mzml``, ``--postms_mode cat``.
	- Additionally, you may provide ``-refseq`` with a ``.fasta`` file containing additional protein sequences you want to exclude from the results.
 	- You may also rescore the proteogenomics results either by providing the argument ``--rescore`` or by running the pipeline on ``rescore`` mode after finishing the ``search``. This is **highly recommended** for proteogenomics searches starting from the three-frame translation of the transcriptome, as the search databases become highly bloated in those cases. The ``--rescore`` argument also requires ``--proteome`` to be specified, and it should be the same proteome in ``.fasta`` format provided in the ``database`` mode.
- Run the ``ribocov`` mode providing ``--fastq``, ``--gtf``, ``--genome_index``, ``--cont_index``, ``--adapter adapter_sequence``, ``--rescored``

## Controlling the FDR
Rp3 provides a few ways to control the way we infer the FDR. 
- You can rescore your proteogenomics results providing the ``--rescore``. This will append the positive results from the first peptide search to the reference proteome and perform a second round of searches with a smaller database, to better control the FDR and reduce false positives.
- If you want to assess the FDR at the protein level, specify ``--proteinFDR`` during the ``search`` mode. By default, the pipeline will assess the FDR at the peptide level.
- Rp3 can also assess the FDR for unannotated microproteins and annotated proteins separately by providing the ``--groupedFDR`` flag. This will likely increase the number of identifications in each group, but will also increase false positives. This flag can be provided during the ``search`` and ``rescore`` mode. If provided in ``search``, it also requires the ``--rescore`` flag.
- MSBooster may also be run after the peptide search and before the post-processing with Percolator by specifying the flag ``--MSBooster``. This will include predicted retention times to the ``.pin`` file from the search with MSFragger used as input for Percolator. This can either increase or reduce identifications, heavily dependant on the dataset. It should make the analysis more robust and we recommend it, however. Requires ``--rescore``. 
- If you have already run the ``search`` or ``postms`` mode, and wish to recalculate the FDR by providing a different set of arguments, you can run either mode by providing ``--recalculateFDR`` along with any desired arguments.
  
## Pipeline output
```
outdir
├── counts	-- contains visualizations and counts table for ribocov mode
│   ├── plots	-- bar plots summarizing ribocov results
│   ├── raw
│   └── rpkm
├── databases
├── homology
│   ├── fasta_to_align
│   └── MSA
├── logs
├── metrics
├── multimappers
├── peptide_search
│   └── group
│       ├── gtf_name_target_database.fasta   
│       └── gtf_name_target_decoy_database.fasta	-- .pin files generated by MSFragger.
├── phylogenetics	-- results from the annotate mode run with the flag --conservation. Formatted for visualization in Evolview.
│   └── blast	-- tblastn results for the identified smORFs
├── post_processing
│   └── group
│       └── db
├── rescore	-- folder containing rescored results
│   ├── databases	-- rescored databases in .fasta format. These contain the microproteins from the first search + the provided proteome
│   ├── peptide_search	-- divided into subdirectories for each group. They contain the .pin files from MSFragger for the second round of searches
│   │   └── group
│   ├── post-processing		-- results from Percolator, divided into subdirectories for each group
│   │   └── group
│   └── summarized_results	-- folder containing the most important results
├── results
│   └── group
│       └── db
├── ribo-seq	-- folder containing results from the ribocov mode
│   ├── alignments	-- alignments in .sam format
│   ├── contaminant_alns	-- reads aligned to rRNAs or tRNAs (unused)
│   ├── sorted_bam_alignments	-- sorted alignments in .bam format
│   └── trimmed_reads	-- trimmed ribo-seq reads without adapters
├── spectra	-- contains annotated spectra for each identified microprotein if the spectrum mode was run
│   └── annotated_spectra	
├── summarized_results		-- main results before rescoring
│   ├── db
│   └── merged		-- merged results from every database (each gtf file generates a db)
└── translation		-- contains the 3 frame translation of the provided GTF files
    └── db
```

## Contact
The Rp3 pipeline is in continuous development. For any questions or suggestions, please contact me at esouza@salk.edu or open an issue on this GitHub page.

## Citation
If you used the Rp3 pipeline in your work, please cite:/
Vieira de Souza E, L. Bookout A, Barnes CA, Miller B, Machado P, Basso LA, Bizarro CV, Saghatelian A. Rp3: Ribosome profiling-assisted proteogenomics improves coverage and confidence during microprotein discovery. Nature Communications. 2024 Aug 9;15(1):6839.

Additionally, depending on which mode you ran and which parameters you used, cite the appropriate software used during the workflow:\
### ``database`` mode
- Cite the paper that presented GTFtoFasta
  - Martinez, T.F., Chu, Q., Donaldson, C., Tan, D., Shokhirev, M.N. and Saghatelian, A., 2020. Accurate annotation of human protein-coding small open reading frames. Nature chemical biology, 16(4), pp.458-468.\
### ``search`` and ``rescore`` mode
- Cite MSFragger:
  - Kong, A.T., Leprevost, F.V., Avtonomov, D.M., Mellacheruvu, D. and Nesvizhskii, A.I., 2017. MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics. Nature methods, 14(5), pp.513-520.\
- Cite Percolator 
  - Käll, L., Canterbury, J.D., Weston, J., Noble, W.S. and MacCoss, M.J., 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nature methods, 4(11), pp.923-925.
- If you used the flag ``--MSBooster``, cite MSBooster:
  - Yang, K.L., Yu, F., Teo, G.C., Li, K., Demichev, V., Ralser, M. and Nesvizhskii, A.I., 2023. MSBooster: improving peptide identification rates using deep learning-based features. Nature Communications, 14(1), p.4539.\
### ``ribocov`` mode
- Cite STAR:
  - Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), pp.15-21.
- Cite featureCounts:
  - Liao, Y., Smyth, G.K. and Shi, W., 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), pp.923-930.
- Cite samtools:
  - Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M. and Li, H., 2021. Twelve years of SAMtools and BCFtools. Gigascience, 10(2), p.giab008.


