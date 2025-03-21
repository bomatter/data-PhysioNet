import os
import pandas as pd
import shutil
import warnings

from tqdm import tqdm
from collections import Counter
from pathlib import Path

from mne_bids import read_raw_bids, write_raw_bids, get_bids_path_from_fname


# Suppress known warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Converting data files to BrainVision.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*There are channels without locations*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Not setting positions of 7 ecg/emg/eog/misc channels*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*Encountered unsupported non-voltage units*")

arousal_labels = [
    "arousal_bruxism",
    "arousal_noise",
    "arousal_plm",
    "arousal_rera",
    "arousal_snore",
    "arousal_spontaneous",
    "resp_centralapnea",
    "resp_cheynestokesbreath",
    "resp_hypopnea",
    "resp_hypoventilation",
    "resp_mixedapnea",
    "resp_obstructiveapnea",
    "resp_partialobstructive",
]

current_dir = Path(__file__).parents[1]
os.chdir(current_dir)

rawdata_dir = Path("rawdata")
derivatives_dir = Path("derivatives")
deriv_root = derivatives_dir / "arousal_labels"

# Delete derivative root if it already exists to ensure a clean slate
if deriv_root.exists():
    print(f"Deleting existing derivative directory: {deriv_root}")
    shutil.rmtree(deriv_root)

# Training data (annotations are only available for training data)
files = [f for f in rawdata_dir.glob("**/sub-tr*.vhdr")]

errors = []
print("Creating derivative files...")
for file in tqdm(files):
    try:
        bp = get_bids_path_from_fname(file)
        raw = read_raw_bids(bp, verbose=False)
        annotations_df = raw.annotations.to_data_frame(time_format=None)
        
        run_counter = Counter()
        for _, annotation in annotations_df.iterrows():
            if annotation["description"] in arousal_labels:
                # Extract the segment corresponding to the annotation
                segment = raw.copy().crop(
                    tmin=annotation["onset"],
                    tmax=annotation["onset"] + annotation["duration"]
                )

                run_counter[annotation["description"]] += 1
                new_bp = bp.copy().update(
                    root = deriv_root,
                    description = annotation["description"].replace("_", ""),
                    run = run_counter[annotation["description"]],
                )

                write_raw_bids(
                    segment,
                    new_bp,
                    overwrite=True,
                    allow_preload=True,
                    format="BrainVision",
                    verbose=False,
                )

                # Update scans.tsv with sleep stage information
                scans_tsv_path = new_bp.fpath.parents[1] / ("sub-" + new_bp.subject + "_scans.tsv")
                scans_tsv = pd.read_csv(scans_tsv_path, sep="\t")
                scans_tsv.loc[scans_tsv.filename == "eeg/" + new_bp.basename, "arousal_label"] = annotation["description"]
                scans_tsv.to_csv(scans_tsv_path, sep="\t", index=False, na_rep="n/a")

    except Exception as e:
        error_msg = f"Error processing file {str(file)}: {e}"
        print(error_msg)
        errors.append(error_msg)
        continue

# Update participants.tsv file: we use the participants.tsv file from the rawdata
# directory but only keep the rows for which we have derivative files
participants_tsv_raw = pd.read_csv(rawdata_dir / "participants.tsv", sep="\t")
participants_tsv_deriv = pd.read_csv(deriv_root / "participants.tsv", sep="\t")
participants_tsv_filtered = participants_tsv_raw[participants_tsv_raw["participant_id"].isin(participants_tsv_deriv["participant_id"])]
participants_tsv_filtered.to_csv(deriv_root / "participants.tsv", sep="\t", index=False, na_rep="n/a")


# Print final summary
print(f"Arousal labels derivative creation completed. {len(files) - len(errors)} / {len(files)} files were successfully processed.")
if errors:
    print("Errors occurred for the following files:")
    print("\n".join(errors))