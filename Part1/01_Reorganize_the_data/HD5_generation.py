#%%
import os
import h5py
import pandas as pd
import numpy as np
import scipy.sparse as sp_sparse

os.chdir('c:/Users/ftroo/Desktop/SC_LCM/00_reorganize_the_data')

#%%
def read_and_merge_UMI_table(file_list):
    output_df = pd.read_csv(file_list[0],sep='\t')
    for f in file_list[1:]:
        temp_df = pd.read_csv(f,sep='\t')
        output_df = output_df.merge(temp_df,left_index=True,right_index=True)
    return(output_df)

#%%
def save_matrix_as_hd5(UMI_table,output_filename):
    # Sparse and reorganize the data
    UMI_CSC_matrix = sp_sparse.csc_matrix(UMI_table.to_numpy())

    matrix_barcodes = UMI_table.columns.values.astype('S')

    matrix_data = UMI_CSC_matrix.data.astype('<i4')
    matrix_indices = UMI_CSC_matrix.indices.astype('<i8')
    matrix_indptr = UMI_CSC_matrix.indptr.astype('<i8')
    matrix_shape = np.array(list(UMI_CSC_matrix.shape),dtype='<i4')

    features_id = UMI_table.index.values.astype('S')
    features_name = UMI_table.index.values.astype('S')

    features_all_tag_keys = np.array([b'genome'])
    feature_type = np.array([b'Gene Expression']*len(features_id))
    features_genome = np.array([b'refgenome_']*len(features_id))

    # Save the data as hd5 file
    output_hf = h5py.File(output_filename,'w')

    output_hf.attrs.create('chemistry_description',b"Single Cell from Ku")
    output_hf.attrs.create('filetype','matrix')
    output_hf.attrs.create('library_ids',np.array([b'Ku_scseq'], dtype='S'))
    output_hf.attrs.create('original_gem_groups',np.array([1], dtype='int64'))
    output_hf.attrs.create('software_version','cellranger-5.0.1')
    output_hf.attrs.create('version',2)
    output_matrix = output_hf.create_group('matrix')
    output_features = output_matrix.create_group('features')

    output_matrix.create_dataset('barcodes',data=matrix_barcodes,compression='gzip')
    output_matrix.create_dataset('data',data=matrix_data,compression='gzip')
    output_matrix.create_dataset('indices',data=matrix_indices,compression='gzip')
    output_matrix.create_dataset('indptr',data=matrix_indptr,compression='gzip')
    output_matrix.create_dataset('shape',data=matrix_shape,compression='gzip')

    output_features.create_dataset('_all_tag_keys',data=features_all_tag_keys,compression='gzip')
    output_features.create_dataset('feature_type',data=feature_type,compression='gzip')
    output_features.create_dataset('genome',data=features_genome,compression='gzip')
    output_features.create_dataset('id',data=features_id,compression='gzip')
    output_features.create_dataset('name',data=features_name,compression='gzip')

    output_hf.close()

    return(0)




#####Eucalpytus grandis#####
#%%
path = 'MARS-seq_UMItable/MARS-seq_UMItable_Eg/'
Eg_UMI_file_list = [path + f for f in os.listdir(path)]
Eg_UMI_table = read_and_merge_UMI_table(Eg_UMI_file_list)
# Eg_UMI_table.shape == (106468, 7296)

#%%
whole_CDS_list = np.array(Eg_UMI_table.index.values)
Eg_CDS_list = pd.read_csv('MARS-seq_UMItable\Eg_CDS.txt',header=None).values.flatten()
set(Eg_CDS_list).issubset(set(whole_CDS_list)) # True

#%%
Eg_UMI_extracted_table = Eg_UMI_table.loc[Eg_CDS_list,]
# Eg_UMI_extracted_table.shape == (36349, 7296)

Eg_UMI_counts = Eg_UMI_extracted_table.sum(0)
Eg_UMI_filtered_table = Eg_UMI_extracted_table.loc[:,Eg_UMI_counts>=100]
# Eg_UMI_filtered_table.shape == (36349, 5491)

#%%
# save_matrix_as_hd5(Eg_UMI_extracted_table,'Eg_UMI.h5')
save_matrix_as_hd5(Eg_UMI_filtered_table,'Eg_UMI_filtered100.h5')



#####Trochodendron aralioides#####
#%%
path = 'MARS-seq_UMItable/MARS-seq_UMItable_Ta/'
Ta_UMI_file_list = [path + f for f in os.listdir(path)]
Ta_UMI_table = read_and_merge_UMI_table(Ta_UMI_file_list)
# Ta_UMI_table.shape == (106468, 3840)

#%%
whole_CDS_list = np.array(Ta_UMI_table.index.values)
Ta_CDS_list = pd.read_csv('MARS-seq_UMItable\Ta_CDS.txt',header=None).values.flatten()
set(Ta_CDS_list).issubset(set(whole_CDS_list)) # True

#%%
Ta_UMI_extracted_table = Ta_UMI_table.loc[Ta_CDS_list,]
# Ta_UMI_extracted_table.shape == (35328, 3840)

Ta_UMI_counts = Ta_UMI_extracted_table.sum(0)
Ta_UMI_filtered_table = Ta_UMI_extracted_table.loc[:,Ta_UMI_counts>=100]
Ta_UMI_filtered_table = Ta_UMI_filtered_table.rename(index=lambda x: x.replace('Ta','Kun'))
# Ta_UMI_filtered_table = Ta_UMI_filtered_table.rename(index=lambda x: x.split('_')[0])
# Ta_UMI_filtered_table.shape == (35328, 1993)

#%%
# save_matrix_as_hd5(Ta_UMI_extracted_table,'Ta_UMI.h5')
save_matrix_as_hd5(Ta_UMI_filtered_table,'Ta_UMI_filtered100.h5')


# %%
