import pandas as pd
import requests
import os
import json

def ASCN_analysis():
    # Read in purity and ploidy data
    purity_ploidy_data = pd.read_csv('purity_ploidy_data.txt', header = 0, sep = '\t')

    # Read in metadata
    with open('ascat3_metadata.json', 'r') as metadata_file:
        metadata = json.load(metadata_file)

    # Assign each ASCAT3 file to the corresponding submitter_id
    submitter_ids = {}
    for file_dict in metadata:
        submitter_ids[file_dict['file_name']] = file_dict['submitter_id']

    tumor_data = {}
    # Loop through the folders, extracting the ASCAT3 file (starts with TCGA)
    for subdir, dirs, files in os.walk('ascat3_allele_specific_data/'):
        for file in files:
            if file.startswith('TCGA'):
                file_path = os.path.join(subdir, file)
                file_data = pd.read_csv(file_path, header = 0, sep = '\t')
                submitter_id = submitter_ids[file]

                # Examine the 11 homozygous deletion cases
                hom_del_ids = ['TCGA-61-2000', 'TCGA-23-2072', 'TCGA-25-2399', 'TCGA-57-1584', 'TCGA-24-1464', 'TCGA-24-1423',
                               'TCGA-13-1817', 'TCGA-24-1603', 'TCGA-04-1336', 'TCGA-09-1669', 'TCGA-24-1549']
                
                for id in hom_del_ids:
                    if id in submitter_id:
                        file_data.to_csv(f'{id}.csv', sep = '\t')

                # Find overlapping segment with NF1 gene
                NF1_start = 31094977
                NF1_end = 31377675
                # 31094977
                # 31377675
                # filtered_file_data = file_data[(file_data['Chromosome'] == 'chr17') & (file_data['Start'] <= NF1_start) & (file_data['Start'] > NF1_start - 10000000) & (file_data['End'] >= NF1_end) & (file_data['End'] < NF1_end + 10000000)]
                filtered_file_data = file_data[(file_data['Chromosome'] == 'chr17') & (file_data['Start'] <= NF1_start) & (file_data['End'] >= NF1_end)]
                if not filtered_file_data.empty:
                    total_copy_number = filtered_file_data.iloc[0]['Copy_Number']
                    major_copy_number = filtered_file_data.iloc[0]['Major_Copy_Number']
                    minor_copy_number = filtered_file_data.iloc[0]['Minor_Copy_Number']
                    region_start = filtered_file_data.iloc[0]['Start']
                    region_end = filtered_file_data.iloc[0]['End']
                    sample_purity_ploidy = purity_ploidy_data[purity_ploidy_data['array'].apply(lambda x: x in submitter_id)]
                    if not sample_purity_ploidy.empty:
                        tumor_purity = sample_purity_ploidy.iloc[0]['purity']
                        tumor_ploidy = sample_purity_ploidy.iloc[0]['ploidy']
                        tumor_data[submitter_id] = {'total_copy_number': total_copy_number, 'major_copy_number': major_copy_number,
                                            'minor_copy_number': minor_copy_number, 'tumor_purity': tumor_purity, 'tumor_ploidy': tumor_ploidy,
                                            'region_start': region_start, 'region_end': region_end}
    
    print(tumor_data)
    print(len(tumor_data))

    non_wildtype = []

    # Get IDs with NF1 short mutation
    NF1_short_mutations = pd.read_csv('NF1_short_mutation_IDs.txt', header = 0, sep = '\t')
    NF1_short_mutation_IDs = NF1_short_mutations['submitter_id'].tolist()

    # Separate NF1 short mutations into separate files based on mutation type
    short_mutations_data = pd.read_csv('TCGA-OV.somaticmutation_wxs.tsv', header = 0, sep = '\t')
    NF1_short_mutations_data = short_mutations_data[short_mutations_data['gene'] == 'NF1']
    for _, row in NF1_short_mutations_data.iterrows():
        mutation_type = row['effect'].split(';')[0]
        sample_id = row['sample'][:-1] # Remove the "A/B/C" at the end of -01 tag
        file_name = f'NF1_{mutation_type}.txt'
        # Don't add sample if duplicate
        if os.path.exists(file_name):
            with open(file_name, 'r') as file:
                existing_samples = [line.strip() for line in file]
        else:
            existing_samples = []
        with open(file_name, 'a') as file:
            if sample_id not in existing_samples:
                non_wildtype.append(sample_id)
                file.write(sample_id + '\n')

    # Filter out those with purity < 30% and with a short mutation in NF1
    filtered_tumor_data = {}
    for tumor in tumor_data:
        if tumor_data[tumor]['tumor_purity'] >= 0.3 and not any(id in tumor for id in NF1_short_mutation_IDs):
            filtered_tumor_data[tumor] = tumor_data[tumor]

    # Save files of IDs for each NF1 genotype group
    # Homozygous deletion
    hom_dels = ['TCGA-04-1336', 'TCGA-09-1669', 'TCGA-13-1817', 'TCGA-23-2072', 'TCGA-24-1423', 'TCGA-24-1464', 'TCGA-24-1549', 'TCGA-24-1603', 'TCGA-25-2399', 'TCGA-57-1584',
                'TCGA-61-2000']
    with open('NF1_hom_del.txt', 'w') as hom_del_file:
        for id in hom_dels:
            non_wildtype.append(id + '-01')
            hom_del_file.write(id + '-01' + '\n')

    # Heterozygous deletion
    with open('NF1_het_del.txt', 'w') as het_del_file:
        for tumor in filtered_tumor_data:
            if filtered_tumor_data[tumor]['total_copy_number'] == 1:
                patient_id = '-'.join(tumor.split('-')[:3])
                non_wildtype.append(patient_id + '-01')
                het_del_file.write(patient_id + '-01' + '\n')

    # Mutation + LOH
    with open('NF1_short_mutations.txt', 'w') as mutations_file:
        for id in NF1_short_mutation_IDs:
            mutations_file.write(id + '-01' + '\n')

    # ================================================================
    # Wild-type
    # ================================================================
    # Classify tumors by relative CN ratio
    summary_rows = []
    loss_samples = []
    neutral_samples = []
    gain_samples = []

    for tumor, vals in filtered_tumor_data.items():
        total_cn = vals['total_copy_number']
        ploidy = vals['tumor_ploidy']

        cn_ratio = total_cn / ploidy
        patient_id = '-'.join(tumor.split('-')[:3])
        sample_id = patient_id + '-01'

        # Skip non-wildtype samples
        if sample_id in non_wildtype:
            continue

        # Classify based on CN ratio
        if cn_ratio < 0.75:
            CN_call = 'loss'
            loss_samples.append(sample_id)
        elif cn_ratio > 1.25:
            CN_call = 'gain'
            gain_samples.append(sample_id)
        else:
            CN_call = 'neutral'
            neutral_samples.append(sample_id)

        # Record full metadata for summary table
        summary_rows.append({'sample_id': sample_id, 'total_copy_number': total_cn, 'tumor_ploidy': ploidy, 'cn_ratio': round(cn_ratio, 3), 'CN_call': CN_call})

    # Save separate lists
    with open('NF1_CN_loss.txt', 'w') as f:
        f.write('\n'.join(loss_samples) + '\n')
    with open('NF1_CN_neutral_wildtype.txt', 'w') as f:
        f.write('\n'.join(neutral_samples) + '\n')
    with open('NF1_CN_gain.txt', 'w') as f:
        f.write('\n'.join(gain_samples) + '\n')

    # Save summary table
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv('NF1_CN_ratio_summary.tsv', sep = '\t', index = False)

def get_wild_type_NF1():
    # Read in purity and ploidy data
    purity_ploidy_data = pd.read_csv('purity_ploidy_data.txt', header = 0, sep = '\t')

    # Read in metadata
    with open('ascat3_metadata.json', 'r') as metadata_file:
        metadata = json.load(metadata_file)

    # Assign each ASCAT3 file to the corresponding submitter_id
    submitter_ids = {}
    for file_dict in metadata:
        submitter_ids[file_dict['file_name']] = file_dict['submitter_id']

    tumor_data = {}
    # Loop through the folders, extracting the ASCAT3 file (starts with TCGA)
    for subdir, dirs, files in os.walk('ascat3_allele_specific_data/'):
        for file in files:
            if file.startswith('TCGA'):
                file_path = os.path.join(subdir, file)
                file_data = pd.read_csv(file_path, header = 0, sep = '\t')
                submitter_id = submitter_ids[file]
    
                # Find overlapping segment with NF1 gene
                NF1_start = 31094977
                NF1_end = 31377675
                filtered_file_data = file_data[(file_data['Chromosome'] == 'chr17') & (file_data['Start'] <= NF1_start) & (file_data['End'] >= NF1_end)]
                if not filtered_file_data.empty:
                    total_copy_number = filtered_file_data.iloc[0]['Copy_Number']
                    major_copy_number = filtered_file_data.iloc[0]['Major_Copy_Number']
                    minor_copy_number = filtered_file_data.iloc[0]['Minor_Copy_Number']
                    sample_purity_ploidy = purity_ploidy_data[purity_ploidy_data['array'].apply(lambda x: x in submitter_id)]
                    if not sample_purity_ploidy.empty:
                        tumor_purity = sample_purity_ploidy.iloc[0]['purity']
                        tumor_ploidy = sample_purity_ploidy.iloc[0]['ploidy']
                        tumor_data[submitter_id] = {'total_copy_number': total_copy_number, 'major_copy_number': major_copy_number,
                                            'minor_copy_number': minor_copy_number, 'tumor_purity': tumor_purity, 'tumor_ploidy': tumor_ploidy}
                        if tumor_data[submitter_id]['total_copy_number'] == 2 and tumor_data[submitter_id]['major_copy_number'] == 1 and tumor_data[submitter_id]['minor_copy_number'] == 1 and \
                        tumor_data[submitter_id]['tumor_ploidy'] >= 1.7 and tumor_data[submitter_id]['tumor_ploidy'] <= 2.3:
                            print(file_data)
    
    # print(tumor_data)
    # print(len(tumor_data))

def get_matched_normal_samples():
    # Only Solid Tissue Normal samples, no blood derived samples
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "project.project_id", "value": ["TCGA-OV"]}},
            {"op": "=", "content": {"field": "samples.sample_type", "value": ["Solid Tissue Normal"]}}
        ]
    }
    fields = [
        "submitter_id",
        "samples.submitter_id",
        "samples.sample_type"
    ]
    params = {
        "filters": filters,
        "fields": ",".join(fields),
        "format": "JSON",
        "size": "10000"
    }
    response = requests.post("https://api.gdc.cancer.gov/cases", json = params)
    data = response.json()["data"]["hits"]

    # Collect all sample barcodes that are Solid Tissue Normal
    sample_ids = []
    for case in data:
        for sample in case.get("samples", []):
            if sample.get("sample_type") == "Solid Tissue Normal":
                sample_ids.append(sample.get("submitter_id"))
    sample_ids = sorted(set(sample_ids))

    # Convert to DataFrame
    df = pd.DataFrame(sample_ids, columns = ["sample_id"])

    # Save to file
    df.to_csv("solid_tissue_normal_ids.txt", sep = "\t", index = False, header = False)

def get_gtex_samples():
    gtex_annotations = "gtex_annotations.txt"
    df = pd.read_csv(gtex_annotations, sep = '\t')
    ovary_df = df[(df["SMTS"] == "Ovary") & (df["SMAFRZE"] == "RNASEQ")]
    ovary_df[["SAMPID"]].to_csv("gtex_ovary_sample_ids.txt", sep = '\t', index = False, header = False)

def match_clonality_and_functional_group():
    expression_df = pd.read_csv("NF1_expression_functional_groups.tsv", sep = "\t")
    clonality_df = pd.read_csv("nf1_clonality_zygosity.tsv", sep = "\t")

    # Make mapping dictionary
    mapping_dict = expression_df.set_index("sample")[["group", "NF1_functional_status"]].to_dict(orient = "index")

    # Match patient IDs to sample IDs
    clonality_df["group"] = None
    clonality_df["NF1_functional_status"] = None
    for i, row in clonality_df.iterrows():
        patient_id = row["patient_id"]
        matches = [sample_id for sample_id in mapping_dict.keys() if patient_id in sample_id]
        if matches:
            match = matches[0]
            clonality_df.at[i, "group"] = mapping_dict[match]["group"]
            clonality_df.at[i, "NF1_functional_status"] = mapping_dict[match]["NF1_functional_status"]

    # === Save merged output ===
    clonality_df.to_csv("nf1_clonality_functional_group.tsv", sep = "\t", index = False)

if __name__ == '__main__':
    ASCN_analysis()
    # get_wild_type_NF1()
    # df = pd.read_csv('TcgaTargetGtex_gene_expected_count.gz', sep = '\t', compression = 'gzip')
    # df.to_pickle('tcga_gtex_rna_seq_counts.pkl')
    # df = pd.read_pickle('tcga_gtex_rna_seq_counts.pkl')
    # print(df)

    # get_matched_normal_samples()
    # get_gtex_samples()
    # match_clonality_and_functional_group()