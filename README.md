# PhysioNet Challenge 2018 Dataset

This is a curated version of the [PhysioNet Challenge 2018 Dataset (v1.0.0)](https://doi.org/10.13026/6phb-r450) in BIDS format.



- Participant IDs correspond to the subject names but without hyphen (e.g. tr030005 instead of tr03-0005) for BIDS compatibility
- Channel names and types are standardised
- Channel locations are set from the standard 10-20 montage
- Age and sex information is added to the participants.tsv file from age-sex.csv
- The official split is added to the participants.tsv file
- Annotations (only available for the train split) are extracted from the matlab and WFDB files and added as events:
  - Sleep stage annotations with labels from the matlab files ("wake", "rem", "nonrem1", ...)
  - Binary arousal regions used in the PhysioNet Challenge 2018 ("arousal" and "non-arousal")
  - Detailed arousal labels extracted from the WFDB files ("arousal_rera", "resp_centralapnea", ...)
- Derivatives with segmented data for the different annotations (`derivatives/sleep_stages`, `derivatives/arousal_binary`, `derivatives/arousal_labels`): segments corresponding to individual annotations (e.g. a continuous segment of a particular sleep stage) are saved as separate files. File names contain the `desc` entity to indicate the corresponding label and the `run` entitiy to enumerate multiple segments of the same label. For example: `sub-tr030146_task-sleep_run-01_desc-rem_eeg.vhdr`. Furthermore, the label for each file is also noted in the associated `scans.tsv` files.




## Reproduce from the Source Data

1. Clone repository

   ```
   git clone https://github.com/bomatter/data-PhysioNet.git PhysioNet
   cd PhysioNet
   ```

2. Install dependencies or use an existing environment with mne, mne-bids, and wfdb installed.
   Example using mamba:

   ```
   mamba create -n bidsify python=3.10 mne>=1.8 mne-bids wfdb
   mamba activate bidsify
   ```

3. (optional) Modify or delete the `.gitignore` file to start tracking the data folders with [DataLad](https://www.datalad.org/).

   ```
   mamba install datalad
   rm .gitignore
   datalad save -m"start tracking data folders"
   ```

4. Download the PhysioNet data to `sourcedata/`

   ```
   mkdir sourcedata
   wget -r -N -c -np -P sourcedata/ https://physionet.org/files/challenge-2018/1.0.0/
   ```

   Save if you are using DataLad to track the data folders:

   ```
   datalad save -m"downloaded sourcedata"
   ```

5. Run the BIDS conversion script.

   ```
   python code/convert_PhysioNet_to_BIDS.py
   ```

   or using DataLad:

   ```
   datalad run \
      -m "run data curation" \
      -i "sourcedata/*" \
      -o "rawdata/*" \
      -o "derivatives/*" \
      "python code/convert_PhysioNet_to_BIDS.py"
   ```

6. (optional) Create derivatives with data segmented according to annotations.

   ```
   python code/create_sleep_stages_derivative.py
   python code/create_arousal_binary_derivative.py
   python code/create_arousal_labels_derivative.py
   ```

   (or use the same format as above for DataLad)