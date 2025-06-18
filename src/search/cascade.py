import os
import pandas as pd
from pyopenms import MzMLFile, MSExperiment

from ..pipeline_config import PipelineStructure

class Cascade(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
        self.firstPassScans = {}
        self.zeroPassScans = {}

    def get_zero_pass_scans(self):
        """
        This will remove scans coming from contaminant sequences.
        """
        self.check_dirs([self.cascadeZeroPassDir, self.cascadeZeroPassMzmlDir])
        self.__retrieve_scans(scan_dict=self.zeroPassScans, pin_dir=self.cascadeZeroPassDir)        
        self.filter_mzml(scan_dict=self.zeroPassScans,
                         mzml_dir=self.args.mzml,  # we start by removing scans from the original mzml files
                         outdir=self.cascadeZeroPassMzmlDir)  # first pass should start here

    def get_first_pass_scans(self):
        print("--Gathering first pass scans")
        self.__retrieve_scans(scan_dict=self.firstPassScans, pin_dir=self.cascadeFirstPassDir)
        self.filter_mzml(scan_dict=self.firstPassScans,
                         mzml_dir=self.cascadeZeroPassMzmlDir,
                         outdir=self.cascadeMzmlDir)

    def __retrieve_scans(self, scan_dict, pin_dir):
        """
        This will feed the provided dictionary with the identified scans to be removed later 
        by calling filter_mzml().
        The pin_dir should contain the pin files with scans to be removed before the next search.
        """
        files = os.listdir(pin_dir)
        total_scans = 0
        for file in files:
            if not file.endswith(".pin"):
                continue
            df = pd.read_csv(
                os.path.join(pin_dir, file),
                sep='\t',
                usecols=range(10),
                comment='#',
                index_col=False
            )
            df["deltCn"] = pd.to_numeric(df["deltCn"], errors="coerce")
            df = df[df["deltCn"] >= 0.1]
            for mzml_id, scan in zip(df["SpecId"], df["ScanNr"]):
                if not isinstance(mzml_id, str):
                    continue
                key = '_'.join(mzml_id.split("_")[:-3])
                scan_dict.setdefault(key, []).append(int(scan))
                total_scans += 1
        print(f"--Found {total_scans} scans from this pass.")

    def filter_mzml(self, scan_dict, mzml_dir=None, outdir=None):
        os.makedirs(outdir, exist_ok=True)
        for fname in os.listdir(mzml_dir):
            if not fname.endswith(self.args.fileFormat):
                continue
            mzml_path = os.path.join(mzml_dir, fname)
            self._filter_single_mzml(mzml_path, outdir, scan_dict=scan_dict)

    def _filter_single_mzml(self, mzml_file, outdir, scan_dict):
        base = os.path.splitext(os.path.basename(mzml_file))[0]
        key = base.rsplit('_', 3)[0]

        # Load full experiment
        exp = MSExperiment()
        MzMLFile().load(mzml_file, exp)

        # Build set of scans to remove
        remove_scans = set(scan_dict.get(key, []))

        # Filter spectra
        out_exp = MSExperiment()
        for spec in exp.getSpectra():
            if spec.getMSLevel() == 2:
                native = spec.getNativeID()
                try:
                    scan = int(native.split('scan=')[-1].split()[0])
                except Exception:
                    # fallback to index-based positioning
                    continue
                if scan in remove_scans:
                    continue
            out_exp.addSpectrum(spec)

        # Write filtered experiment
        filtered_file = os.path.join(outdir, f"{base}_filtered.mzML")
        MzMLFile().store(filtered_file, out_exp)
        print(f"--Wrote filtered mzML via pyOpenMS to {filtered_file}")


    def concatenate_pin_files(self):
        print(f"--Concatenating pin files from first and second pass of cascade search...")
        cmd = f'cat {self.cascadeFirstPassDir}/*pin {self.cascadeSecondPassDir}/*pin > {self.searchDir}/cascade_search_unfixed.pin'
        os.system(cmd)

        cmd = f"grep -v 'SpecId' {self.searchDir}/cascade_search_unfixed.pin > {self.searchDir}/cascade_search_filtered.pin"
        os.system(cmd)

        cmd = f"awk 'FNR<2' {self.searchDir}/cascade_search_unfixed.pin | cat - {self.searchDir}/cascade_search_filtered.pin > {self.searchDir}/cascade_search.pin"
        os.system(cmd)