import os
import sys

from tqdm import tqdm


class Iterator:
    def __init__(self, mzml_folders, outdir, proteome, rp3_folder_to_copy=None):
        self.mzmlFolders = mzml_folders
        self.outdir = outdir
        self.proteome = proteome
        # self.rp3Folder = rp3_folder_to_copy
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    def iterate_tmt_searches(self):
        batches = os.listdir(self.mzmlFolders)
        for mzml in tqdm(batches):
            if mzml.startswith("b"):
                # if not os.path.exists(f'{self.outdir}/{mzml}'):
                #     os.mkdir(f'{self.outdir}/{mzml}')

                # outdir = self.copy_db(mzml)
                outdir = f'{self.outdir}/{mzml}'

                checker = f'{outdir}/summarized_results/merged/filtered_rescored_microproteins_150.fasta'

                if not os.path.exists(checker):
                    cmd = (f'python3 ~/PycharmProjects/rp3/rp3.py database --outdir {outdir} --proteome /home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta '
                           f'--external_database /ceph/pbla/brendan/assemblies/cpm05/2024_12_21_cpm05_brain_espresso_shortstop_no_duplicates.fa '
                           f'--skip_translation --genomeAssembly hg38')
                    os.system(cmd)
                    cmd = (f'python3 ~/PycharmProjects/rp3/rp3.py search --outdir {outdir} --mzml {self.mzmlFolders}/{mzml} --proteome {self.proteome} '
                           f'--postms_mode sep --threads 128 --mod "229.16293_S_1" --msBooster --fragment_mass_tolerance 20')
                    os.system(cmd)
                else:
                    print(f"--BLING Skipping {mzml}. ")

    def iterate_tmt_searches_proteogenomics(self):
        batches = os.listdir(self.mzmlFolders)
        for mzml in tqdm(batches):
            if mzml.startswith("b"):
                # if not os.path.exists(f'{self.outdir}/{mzml}'):
                #     os.mkdir(f'{self.outdir}/{mzml}')

                # outdir = self.copy_db(mzml)
                outdir = f'{self.outdir}/{mzml}'

                checker = f'{outdir}/summarized_results/merged/filtered_rescored_microproteins_150.fasta'

                if not os.path.exists(checker):
                    cmd = (f'python3 ~/PycharmProjects/rp3/rp3.py database --outdir {outdir} --proteome /home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta '
                           f'--gtf_folder /ceph/pbla/brendan/assemblies/cpm05/ESPRESSO_syn52047893_hg38NCBI_N2_R0_updated.sorted_medianCPM05_GENE_CLEANED.gtf'
                           f' --genomeAssembly hg38')
                    os.system(cmd)
                    cmd = (f'python3 ~/PycharmProjects/rp3/rp3.py search --outdir {outdir} --mzml {self.mzmlFolders}/{mzml} --proteome {self.proteome} '
                           f'--postms_mode sep --threads 128 --mod "229.16293_S_1" --msBooster --fragment_mass_tolerance 20')
                    os.system(cmd)
                else:
                    print(f"--BLING Skipping {mzml}. ")

    # def copy_db(self, batch):
    #     outdir = f'{self.outdir}/{batch}'
    #
    #     print(f"--Copying DBs to {outdir}")
    #     cmd = f'cp -r {self.rp3Folder}/.* {outdir}/.'
    #     os.system(cmd)
    #     return outdir

if __name__ == '__main__':
    """
    Brendan's brain smORFome
    rp3 proteogenomics
    CPM cutoff 2
    """
    # rp3 = Iterator(mzml_folders='/ceph/pbla/brendan/rosmap_tmt_ms/round1',
    #                outdir='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/proteogenomics_DB',
    #                rp3_folder_to_copy='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/proteogenomics_DB/241219_brainSmorfome',
    #                proteome='/home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta')
    # rp3.iterate_tmt_searches()

    """
    Brendan's brain smORFome
    rp3 shortStop
    CPM cutoff 0.5
    """
    # rp3 = Iterator(mzml_folders='/ceph/pbla/brendan/rosmap_tmt_ms/round1',
    #                outdir='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/cpmCutoff_005/shortStop_Rp3/25-01-02_shortStop_DB',
    #                rp3_folder_to_copy='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/cpmCutoff_005/shortStop_Rp3/241224_shortStopDB',
    #                proteome='/home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta')
    # rp3.iterate_tmt_searches()
    """
    proteogenomics
    """
    # rp3 = Iterator(mzml_folders='/ceph/pbla/brendan/rosmap_tmt_ms/round1',
    #                outdir='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/cpmCutoff_005/shortStop_Rp3/25-01-03_proteogenomics',
    #                proteome='/home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta')
    # rp3.iterate_tmt_searches_proteogenomics()

    """ 
    proteogenomics 
    round2 
    """
    rp3 = Iterator(mzml_folders='/ceph/pbla/brendan/rosmap_tmt_ms/round2',
                   outdir='/ceph/pbla/eduardo/brendan_brainSmorfome/rp3_runs/cpmCutoff_005/shortStop_Rp3/25-01-07_proteogenomics_round2',
                   proteome='/home/microway/proteomes/homoSapiens_uniprotkb_proteome_UP000005640_AND_revi_2024_07_30.fasta')
    rp3.iterate_tmt_searches_proteogenomics()

