import os
import sys

import pandas as pd
from pymzml import run

from ..pipeline_config import PipelineStructure

class Cascade(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)

        self.firstPassScans = {}

    def get_first_pass_scans(self):
        files = os.listdir(self.cascadeFirstPassDir)
        total_scans = 0
        for file in files:
            df = pd.read_csv(os.path.join(self.cascadeFirstPassDir, file), sep='\t')
            df = df[df["deltCn"] >= 0.1]
            mzml_files = df["SpecId"].tolist()
            scans = df["ScanNr"].tolist()
            for mzml, scan in zip(mzml_files, scans):
                renamed = '_'.join(mzml.split("_")[:-3])
                if renamed not in self.firstPassScans:
                    self.firstPassScans[renamed] = []
                self.firstPassScans[renamed].append(scan)
                total_scans += 1
        print(f"--Found {total_scans} from the reference proteome in the first-pass Comet search.")

    def filter_mzml(self, mzml_dir=None, outdir=None):
        files = os.listdir(mzml_dir)
        for file in files:
            if file.endswith(self.args.fileFormat):
                self._filter_single_mzml(mzml_file=os.path.join(mzml_dir, file),
                                         outdir=outdir)

    def _filter_single_mzml(self, mzml_file, outdir):
        mzml = run.Reader(mzml_file)
        filtered = []   
        for spec in mzml:
            if spec.ms_level == 2 and spec.scan_number not in self.firstPassScans.get(mzml_file, []):
                filtered.append(spec)
        if filtered:
            filtered_file = os.path.join(outdir, f"{mzml_file.replace('.mzML', '_filtered.mzML')}".split("/")[-1])
            with run.Writer(filtered_file) as writer:
                for spec in filtered:
                    writer.write(spec)
