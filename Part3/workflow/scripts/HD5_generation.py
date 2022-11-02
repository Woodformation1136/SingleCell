#%%
import os
import h5py
import pandas as pd
import numpy as np
import scipy.sparse as sp_sparse
import argparse


#%%
# Define functions
# Read and merge UMI tables
def read_and_merge_UMI_table(file_list: str):
    output_df = pd.read_csv(
        file_list[0], sep = '\t'
    )
    for f in file_list[1:]:
        output_df = output_df.merge(
            pd.read_csv(f, sep='\t'),
            left_index = True,
            right_index = True
        )
    output_df.columns = [x.replace('_', '') for x in output_df.columns]
    return output_df

# To save UMI table as hd5 file
def save_matrix_as_hd5(UMI_table, output_filename: str):
    # Sparse and reorganize the data
    UMI_CSC_matrix = sp_sparse.csc_matrix(UMI_table.to_numpy())

    matrix_barcodes = UMI_table.columns.values.astype('S')

    matrix_data = UMI_CSC_matrix.data.astype('<i4')
    matrix_indices = UMI_CSC_matrix.indices.astype('<i8')
    matrix_indptr = UMI_CSC_matrix.indptr.astype('<i8')
    matrix_shape = np.array(list(UMI_CSC_matrix.shape), dtype = '<i4')

    features_id = UMI_table.index.values.astype('S')
    features_name = UMI_table.index.values.astype('S')

    features_all_tag_keys = np.array([b'genome'])
    feature_type = np.array([b'Gene Expression'] * len(features_id))
    features_genome = np.array([b'refgenome_'] * len(features_id))

    # Create dict for output
    dict_hf = {
        'chemistry_description': b"Single Cell from Ku",
        'filetype': 'matrix',
        'library_ids': np.array([b'Ku_scseq'], dtype='S'),
        'original_gem_groups': np.array([1], dtype='int64'),
        'software_version': 'cellranger-7.0.0',
        'version': 2
    }
    dict_matrix = {
        'barcodes': matrix_barcodes,
        'data': matrix_data,
        'indices': matrix_indices,
        'indptr': matrix_indptr,
        'shape': matrix_shape
    }
    dict_features = {
        '_all_tag_keys': features_all_tag_keys,
        'feature_type': feature_type,
        'genome': features_genome,
        'id': features_id,
        'name': features_name
    }

    # Save the data as hd5 file
    output_hf = h5py.File(output_filename, 'w')
    output_matrix = output_hf.create_group('matrix')
    output_features = output_matrix.create_group('features')

    for (key, value) in dict_hf.items():
        output_hf.attrs.create(
            name = key, data = value
        )
    for (key, value) in dict_matrix.items():
        output_matrix.create_dataset(
            name = key, data = value, compression = 'gzip'
        )
    for (key, value) in dict_features.items():
        output_features.create_dataset(
            name = key, data = value, compression = 'gzip'
        )

    output_hf.close()

# Extract UMI table using CDS list and save as hd5 file
def read_to_save(
    UMI_folder_path: str,
    CDS_file_path: str,
    hd5_file_path: str,
    UMI_criteria: int = 100
):
    UMI_file_list = [
        UMI_folder_path + '/' +
        f for f in sorted(os.listdir(UMI_folder_path))
    ]
    UMI_table = read_and_merge_UMI_table(UMI_file_list)

    CDS_list = pd.read_csv(
        CDS_file_path, header = None
    ).values.flatten()

    UMI_extracted_table = UMI_table.loc[CDS_list,]

    UMI_counts = UMI_extracted_table.sum(axis = 0)
    UMI_filtered_table = UMI_extracted_table.loc[:, UMI_counts >= UMI_criteria]

    if "Ta_CDS" in CDS_file_path:
        UMI_filtered_table = UMI_filtered_table.rename(
            index=lambda x: "evm.TU." + x
        )

    save_matrix_as_hd5(UMI_filtered_table, hd5_file_path)


#%%
parser = argparse.ArgumentParser()
parser.add_argument(
    '--umi_dir', help='The directory of input UMI matrix from MARSseq'
)
parser.add_argument(
    '--cds', help='The path of input CDS list'
)
parser.add_argument(
    '--output_hd5', help='The path of output hd5 file'
)
parser.add_argument(
    '--umi_criteria', help='Only cells with at least this number of UMI counts would be included', type=int
)
args = parser.parse_args()

#%%
read_to_save(
    UMI_folder_path=args.umi_dir,
    CDS_file_path=args.cds,
    hd5_file_path=args.output_hd5,
    UMI_criteria=args.umi_criteria
)


