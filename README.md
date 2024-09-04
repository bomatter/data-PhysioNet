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




## Reproduce from the Source Data

1. Clone repository

   ```
   git clone https://github.com/bomatter/data-PhysioNet.git PhysioNet
   cd PhysioNet
   ```

2. Install dependencies or use an existing environment with mne, mne-bids, and wfdb installed.
   Example using mamba:

   ```
   mamba create -n bidsify python=3.10 mne mne-bids wfdb
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