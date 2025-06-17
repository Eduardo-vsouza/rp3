import os
import pandas as pd
from pyopenms import MzMLFile, MSExperiment

from ..pipeline_config import PipelineStructure

class Cascade(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
        self.firstPassScans = {}

    def get_first_pass_scans(self):
        print("--Gathering first pass scans")
        files = os.listdir(self.cascadeFirstPassDir)
        total_scans = 0
        for file in files:
            if not file.endswith(".pin"):
                continue
            df = pd.read_csv(
                os.path.join(self.cascadeFirstPassDir, file),
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
                self.firstPassScans.setdefault(key, []).append(int(scan))
                total_scans += 1
        print(f"--Found {total_scans} scans from the first pass.")

    def filter_mzml(self, mzml_dir=None, outdir=None):
        os.makedirs(outdir, exist_ok=True)
        for fname in os.listdir(mzml_dir):
            if not fname.endswith(self.args.fileFormat):
                continue
            mzml_path = os.path.join(mzml_dir, fname)
            self._filter_single_mzml(mzml_path, outdir)

    def _filter_single_mzml(self, mzml_file, outdir):
        base = os.path.splitext(os.path.basename(mzml_file))[0]
        key = base.rsplit('_', 3)[0]

        # Load full experiment
        exp = MSExperiment()
        MzMLFile().load(mzml_file, exp)

        # Build set of scans to remove
        remove_scans = set(self.firstPassScans.get(key, []))

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
        cmd = f'cat {self.cascadeFirstPassDir}/*pin {self.cascadeSecondPassDir}/*pin > {self.searchOutdir}/cascade_search_unfixed.pin'
        os.system(cmd)

        cmd = f"grep -v 'SpecId' {self.searchOutdir}/cascade_search_unfixed.pin > {self.searchOutdir}/cascade_search_filtered.pin"
        os.system(cmd)

        cmd = f"awk 'FNR<2' {self.searchOutdir}/cascade_search_unfixed.pin | cat - {self.searchOutdir}/cascade_search_filtered.pin > {self.searchOutdir}/cascade_search.pin"
        os.system(cmd)