import os
import re
import shutil
import numpy as np
import pandas as pd

import warnings
from warnings import warn

from pathlib import Path
from tqdm import tqdm
from wfdb import rdrecord, rdann
from pymatreader import read_mat
from sklearn.model_selection import train_test_split

import mne
from mne import create_info
from mne.io import RawArray
from mne_bids import write_raw_bids, BIDSPath

# Suppress known warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Converting data files to BrainVision.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message='.*Encountered data in "double" format.*')
warnings.filterwarnings("ignore", category=UserWarning, message=".*Encountered unsupported non-voltage units*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*No events found or provided.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*annotation\(s\) that were expanding outside the data range.*")

# Rename channel names for better interoperability with standard montages in MNE
channel_name_mapping = {
    'F3-M2': 'F3',
    'F4-M1': 'F4',
    'C3-M2': 'C3',
    'C4-M1': 'C4',
    'O1-M2': 'O1',
    'O2-M1': 'O2',
    'E1-M2': 'EOG',
    'Chin1-Chin2': 'EMG',
    'ABD': 'ABD',
    'CHEST': 'CHEST',
    'AIRFLOW': 'AIRFLOW',
    'SaO2': 'SaO2',
    'ECG': 'ECG'
}

# Map channel names to channel types
channel_type_mapping = {
    'F3-M2': 'eeg',
    'F4-M1': 'eeg',
    'C3-M2': 'eeg',
    'C4-M1': 'eeg',
    'O1-M2': 'eeg',
    'O2-M1': 'eeg',
    'E1-M2': 'eog',  # left EOG
    'Chin1-Chin2': 'emg',
    'ABD': 'misc',  # abdominal
    'CHEST': 'misc',
    'AIRFLOW': 'misc',
    'SaO2': 'misc',  # oxygen saturation
    'ECG': 'ecg'
}

# Arousal labels used in WFDB (.arousal files) files
# Note: in the files, these labels appear with parentheses to indicate
# on and offsets, e.g. "(arousal_rera" and "arousal_rera)"
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

def parse_samplewise_annotations(samplewise_annotations, category, sampling_frequency):
    """
    Creates lists of onsets, durations, and descriptions from numpy
    arrays of sample-wise annotations and the provided category.
    
    sample_wise_annotations: array of shape (n_samples,) with 1s
    indicating samples of the provided category and 0s elsewhere.

    sample_frequency is used to convert sample indices to seconds
    as required by MNE.
    """
    changes = np.diff(samplewise_annotations, prepend=0, append=0)
    onsets = (np.where(changes == 1)[0] + 1) / sampling_frequency
    offsets = (np.where(changes == -1)[0] + 1) / sampling_frequency
    durations = offsets - onsets
    descriptions = [category] * len(onsets)
    
    return onsets, durations, descriptions

def parse_arousal_annotations(annotations, label, sampling_frequency):
    """
    Parses arousal annotations from WFDB annotations object and returns
    onsets, durations, and descriptions for the provided label.
    """

    # Filter out onsets for the specific label, e.g. "(arousal_rera"
    onsets = sorted([
        t / sampling_frequency for desc, t
        in zip(annotations.aux_note, annotations.sample)
        if desc == "(" + label
    ])

    # Filter out offsets for the specific label, e.g. "arousal_rera)"
    offsets = sorted([
        t / sampling_frequency for desc, t 
        in zip(annotations.aux_note, annotations.sample)
        if desc == label + ")"
    ])

    durations = [off - on for on, off in zip(onsets, offsets)]
    descriptions = [label] * len(onsets)

    return onsets, durations, descriptions


current_dir = Path(__file__).parents[1]
os.chdir(current_dir)

sourcedata_dir = Path("sourcedata/physionet.org/files/challenge-2018/1.0.0/")
rawdata_dir = Path("rawdata")  # Where the BIDSified data should be written to

# Delete rawdata directory if it already exists to ensure a clean slate
if rawdata_dir.exists():
    print(f"Deleting existing rawdata directory: {rawdata_dir}")
    shutil.rmtree(rawdata_dir)

# Collect file names
pattern = re.compile(r'(tr|te)\d{2}-\d{4}\.mat$')
files = [f for f in sourcedata_dir.glob("**/*.mat") if pattern.match(f.name)]
assert len(files) == 1983, f"Expected 1982 files (994 training + 989 test), found {len(files)}"

errors = []
print("Converting PhysioNet files to BIDS...")
for file in tqdm(files):
    try:
        subject = file.stem.replace("-", "")  # MNE BIDS does not allow dashes in subject IDs

        if subject.startswith("tr"):
            split = "train"
        elif subject.startswith("te"):
            split = "test"
        else:
            raise ValueError(f"Unknown split for subject {subject}")

        # Create BIDSPath
        bids_path = BIDSPath(
            root=rawdata_dir,
            subject=subject,
            task="sleep",
        )

        # Load data from WFDB file
        data = rdrecord(file.with_suffix(''))
        
        # Extract information
        sampling_frequency = data.fs
        channel_names = data.sig_name
        channel_units = data.units

        # Set channel types and standardize names
        channel_types = [channel_type_mapping[ch] for ch in channel_names]
        channel_names = [channel_name_mapping[ch] for ch in channel_names]

        # Convert from microvolts to volts (MNE expects volts)
        scaling_factors = np.ones(len(channel_names)).reshape(-1, 1)
        for i, (type, unit) in enumerate(zip(channel_types, channel_units)):
            if type in ['eeg', 'eog', 'emg', 'ecg']:
                if unit == 'uV':
                    scaling_factors[i] = 1e-6
                elif unit == 'mV':
                    scaling_factors[i] = 1e-3
                elif unit != 'V':
                    raise ValueError(f"Unknown unit {unit} for channel {channel_names[i]}")

        # Create MNE raw object
        raw = RawArray(
            data.p_signal.T * scaling_factors,
            create_info(
                ch_names=channel_names,
                sfreq=sampling_frequency,
                ch_types=channel_types,
            ),
            verbose=False
        )

        # Add annotations from sample-wise matlab annotations
        if split == "train":
            annotations_mat = read_mat(file.with_stem(file.stem + "-arousal"))

            onsets = []
            durations = []
            descriptions = []

            # Sleep stage annotations
            for label, values in annotations_mat["data"]["sleep_stages"].items():
                if label != "undefined":  # no need to create annotations for undefined class
                    on, dur, desc = parse_samplewise_annotations(values, label, sampling_frequency)
                    onsets.extend(on)
                    durations.extend(dur)
                    descriptions.extend(desc)

            # Binary arousal annotations from sample-wise matlab annotations
            # Note: annotations["data"]["arousals"] contains values 0 (non-arousal),
            # 1 (arousal), and -1 (not scored in competition)

            # Arousal
            on, dur, desc = parse_samplewise_annotations(
                annotations_mat["data"]["arousals"] == 1,
                "arousal",
                sampling_frequency
            )
            onsets.extend(on)
            durations.extend(dur)
            descriptions.extend(desc)

            # Non-arousal
            on, dur, desc = parse_samplewise_annotations(
                annotations_mat["data"]["arousals"] == 0,
                "non-arousal",
                sampling_frequency
            )
            onsets.extend(on)
            durations.extend(dur)
            descriptions.extend(desc)

            # Add detailed arousal annotations from WFDB (.arousal) files
            annotations_rdann = rdann(str(file.with_suffix('')), 'arousal')

            # Loop over all arousal labels in the annotation file
            for label in {l.strip("()") for l in annotations_rdann.aux_note if l.strip("()") in arousal_labels}:
                on, dur, desc = parse_arousal_annotations(annotations_rdann, label, sampling_frequency)
                onsets.extend(on)
                durations.extend(dur)
                descriptions.extend(desc)

            # Create MNE annotations object and add to raw object
            annotations = mne.Annotations(
                onset=onsets,
                duration=durations,
                description=descriptions
            )

            # Remove annotations with onset outside of the recording
            remove_indices = np.where(annotations.onset > raw.times[-1])[0]
            if len(remove_indices):
                warn(f"Removing {len(remove_indices)} annotations with onsets outside of the recording.")
                annotations.delete(remove_indices)

            # Fix annotations with durations that extend beyond the recording
            remove_indices = np.where(annotations.onset + annotations.duration > raw.times[-1])[0]
            for idx in remove_indices:
                a = annotations[idx]
                # Check if the duration is off by less than 20 samples (0.1 s at a sampling frequency of 200 Hz).
                # Discrepancies up to this tolerance are fixed, otherwise we raise an error.
                dt = 1 / raw.info["sfreq"]
                tolerance = 20 * dt
                if (a["onset"] + a["duration"]) - raw.times[-1] < tolerance:
                    a["duration"] = raw.times[-1] - a["onset"]  # Adjust duration to end of recording
                    annotations.delete(idx)
                    annotations.append(onset=a["onset"], duration=a["duration"], description=a["description"])
                else:
                    raise RuntimeError(f"Annotation {a} extends beyond the recording and exceeds the tolerance of 10 samples.")

            raw.set_annotations(annotations)

        # Add channel locations
        montage = mne.channels.make_standard_montage('standard_1020')
        raw.set_montage(montage)

        # Write BIDS
        write_raw_bids(
            raw,
            bids_path=bids_path,
            overwrite=True,
            verbose=False,
            allow_preload=True,
            format="BrainVision"
        )

    except Exception as e:
        print(f"Error processing file {file}: {e}")
        errors.append(str(file))
        continue

# Curate participants.tsv file
participants_tsv = pd.read_csv(rawdata_dir / "participants.tsv", sep="\t")

# Read age and sex info csv and format so we can merge with participants.tsv
age_sex_info = pd.read_csv(sourcedata_dir / "age-sex.csv")
age_sex_info["participant_id"] = "sub-" + age_sex_info["Record"].str.replace("-", "")
age_sex_info = age_sex_info.drop(columns=["Record"])
age_sex_info = age_sex_info.rename(columns={"Age": "age", "Sex": "sex"})

# Merge age and sex info into participants.tsv
participants_tsv = participants_tsv.drop(columns=["age", "sex", "hand", "weight", "height"])  # these columns only contain NaNs
participants_tsv = participants_tsv.merge(age_sex_info, on="participant_id", how="left")
participants_tsv.reindex(columns=["participant_id", "age", "sex"])  # reorder columns

# Standardise sex column
participants_tsv['sex'] = participants_tsv['sex'].str.upper()

# Add official split information
participants_tsv["split"] = np.where(
    participants_tsv["participant_id"].str.startswith("sub-tr"), "train", "test"
)

# Split train into train and validation, stratifying by age and sex
participants_train = participants_tsv[participants_tsv['split'] == 'train']
age_bins = pd.cut(participants_train['age'], bins=[0, 20, 30, 40, 50, 60, 70, 80, 100])
stratification_vars = age_bins.astype(str) + '_' + participants_train['sex']
train, val = train_test_split(
    participants_train.participant_id,
    stratify=stratification_vars,
    test_size=0.1, random_state=9292, shuffle=True
)

participants_tsv['train_val_test_split'] = participants_tsv['participant_id'].apply(
    lambda x: 'train' if x in train.values else 'val' if x in val.values else 'test'
)

# Check that all participants in the official test split are also test in the train_val_test_split
assert (participants_tsv.loc[participants_tsv['split'] == 'test', 'train_val_test_split'] == 'test').all()

# Since labels are only available for the training set, we also create
# a labled_train_val_test_split column that only contains participants
# with labels, i.e. those in the official train split.
participants_train = participants_tsv[participants_tsv['split'] == 'train']
age_bins = pd.cut(participants_train['age'], bins=[0, 25, 50, 75, 100])
stratification_vars = age_bins.astype(str) + '_' + participants_train['sex']
train, val_and_test = train_test_split(
    participants_train.participant_id,
    stratify=stratification_vars,
    test_size=0.15, random_state=9292, shuffle=True
)

val, test = train_test_split(
    val_and_test,
    stratify=val_and_test.map(lambda pid: stratification_vars[participants_train.participant_id == pid].values[0]),
    test_size=0.5, random_state=9292, shuffle=True
)

participants_tsv['labeled_train_val_test_split'] = 'NaN'
participants_tsv.loc[participants_tsv['participant_id'].isin(train), 'labeled_train_val_test_split'] = 'train'
participants_tsv.loc[participants_tsv['participant_id'].isin(val), 'labeled_train_val_test_split'] = 'val'
participants_tsv.loc[participants_tsv['participant_id'].isin(test), 'labeled_train_val_test_split'] = 'test'


participants_tsv.to_csv(rawdata_dir / "participants.tsv", sep="\t", index=False, na_rep="n/a")


print(f"BIDS conversion completed. {len(files)- len(errors)}/{len(files)} files were successfully processed.")
if errors:
    print("Errors occurred for the following files:")
    print("\n".join(errors))
