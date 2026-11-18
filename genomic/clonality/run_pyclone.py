import os
import pandas as pd
import subprocess

# Define input and output directories
input_dir = 'pyclone_input/'
output_dir = 'pyclone_output/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop over all TSV files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.tsv'):
        patient_id = filename.replace('.tsv', '')
        input_file_path = os.path.join(input_dir, filename)
        
        # Read in the TSV and extract tumor_purity
        df = pd.read_csv(input_file_path, sep = '\t')
        tumor_purity = df['tumor_purity'].iloc[0]

        # Prepare output subdirectory for this patient
        patient_output_dir = os.path.join(output_dir, patient_id)
        if not os.path.exists(patient_output_dir):
            os.makedirs(patient_output_dir)

        # Build and run PyClone command
        command = [
            'pyclone', 'run_analysis_pipeline',
            '--in_files', input_file_path,
            '--working_dir', patient_output_dir,
            '--tumour_contents', str(tumor_purity)
        ]
        print(patient_id)
        subprocess.call(command)