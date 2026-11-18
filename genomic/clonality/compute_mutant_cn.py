import os
import pandas as pd

input_dir = 'pyclone_input/'
output_dir = 'pyclone_output/'
results = []

# Loop through each input file
for file in os.listdir(input_dir):
    patient_id = file.replace('.tsv', '')
    input_path = os.path.join(input_dir, file)

    # Load input data and extract NF1 mutation row
    input_df = pd.read_csv(input_path, sep = '\t')
    nf1_row = input_df[input_df['mutation_id'].str.contains('NF1_chr17')]

    # Extract values
    mutation_id = nf1_row['mutation_id'].values[0]
    ref_counts = nf1_row['ref_counts'].values[0]
    var_counts = nf1_row['var_counts'].values[0]
    normal_cn = nf1_row['normal_cn'].values[0]
    minor_cn = nf1_row['minor_cn'].values[0]
    major_cn = nf1_row['major_cn'].values[0]
    tumor_purity = nf1_row['tumor_purity'].values[0]

    # Find matching PyClone output
    pyclone_subdir = os.path.join(output_dir, patient_id, 'tables', 'loci.tsv')
    if not os.path.exists(pyclone_subdir):
        print(f'No PyClone output for {patient_id}')
        continue

    loci_df = pd.read_csv(pyclone_subdir, sep = '\t')
    loci_row = loci_df[loci_df['mutation_id'] == mutation_id]
    if loci_row.empty:
        print(f'No matching mutation_id in PyClone output for {patient_id}')
        continue

    cellular_prevalence = loci_row['cellular_prevalence'].values[0]

    # Calculate variables
    VAF_obs = var_counts / (ref_counts + var_counts)
    total_cn = major_cn + minor_cn
    mutant_cn = (VAF_obs / tumor_purity) * ((tumor_purity * total_cn) + ((1 - tumor_purity) * normal_cn)) / cellular_prevalence

    # Store result
    results.append({
        'patient_id': patient_id,
        'mutation_id': mutation_id,
        'cellular_prevalence': cellular_prevalence,
        'tumor_purity': tumor_purity,
        'VAF_obs': VAF_obs,
        'major_cn': major_cn,
        'minor_cn': minor_cn,
        'mutant_cn': mutant_cn
    })

# Write output file
summary_df = pd.DataFrame(results)
summary_df.to_csv('nf1_clonality_zygosity.tsv', sep = '\t', index = False)