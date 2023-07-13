origin = open('/data/jiayi/IVT_results/normalU/eventalign_normalU.txt','r')
origin2 = open('','r')
new = open('/data/jiayi/IVT_results/Root','w')

for idx, i in enumerate(origin):
    if idx == 0:
        i = i.replace("\n", '') + '\t' + 'modification_status\n'
        new.write(i)
        print(i)

        continue
    else:
        i = i.replace("\n", '')
        new.write(i + '\t' + "0\n")

for idx, i in enumerate(origin2):
    if idx == 0:
        continue
    else:
        i = i.replace("\n", "")
        new.write(i + "\t" + "1\n")
new.close()
origin.close()
origin2.close()

# import numpy as np
# import pandas as pd
# pd.set_option('display.max_columns',20)
#
#
# def index(eventalign_result,pos_start):
#    eventalign_result = eventalign_result.set_index(['contig','read_index'])
#    # print(eventalign_result)
#    pos_end=pos_start
#    with open('/data/jiayi/IVT_results/check/Check_2*1000/eventalign.index','a') as f_index:
#
#        # print(eventalign_result.index)
#        for index, modification_status in zip(list(dict.fromkeys(eventalign_result.index)),eventalign_result['modification_status']):
#            print(index)
#            print(modification_status)
#            transcript_id,read_index = index
#            pos_end += eventalign_result.loc[index]['line_length'].sum()
#            f_index.write('%s,%d,%d,%d,%d\n' %(transcript_id,read_index,pos_start,pos_end,modification_status))
#            pos_start = pos_end
#        f_index.close()
#
# with open('/data/jiayi/IVT_results/check/Check_2*1000/eventalign.index','w') as f:
# # f.write('transcript_id,read_index,pos_start,pos_end\n') # header
#         f.write('transcript_id,read_index,pos_start,pos_end,modification_status\n')
# eventalign_file = open('/data/jiayi/IVT_results/check/Check_2*1000/eventalign1.txt','r')
# pos_start = len(eventalign_file.readline()) #remove header
# chunk_split = None
# index_features = ['contig','read_index','line_length','modification_status']
# for chunk in pd.read_csv('/data/jiayi/IVT_results/check/Check_2*1000/eventalign1.txt', chunksize=1000000,sep='\t'):
#         chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
#         # print("chunk['read_index']",chunk['read_index'])
#         # print("chunk.iloc[-1]['read_index']",chunk.iloc[-1]['read_index'])
#         # print(chunk_complete)
#         chunk_concat = pd.concat([chunk_split,chunk_complete])
#         # print("chunk_concat.index",chunk_concat.index)
#         chunk_concat_size = len(chunk_concat.index)
#         ## read the file at where it left off because the file is opened once ##
#         lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
#         # print("lines",lines)
#         chunk_concat.loc[:, 'line_length'] = np.array(lines)
#         # print(chunk_concat[index_features])
#         index(chunk_concat[index_features],pos_start)




# eventalign_result = pd.read_csv('/data/jiayi/IVT_results/check/Check_2*1000/eventalign1.txt',delimiter='\t',names=['contig','position','reference_kmer','read_index','strand','event_index','event_level_mean','event_stdv','event_length','model_kmer','model_mean','model_stdv','standardized_level','start_idx','end_idx','modification_status'])
#
# cond_successfully_eventaligned = eventalign_result['reference_kmer'] == eventalign_result['model_kmer']
#
# eventalign_result = eventalign_result[cond_successfully_eventaligned]
#
# keys = ['read_index','contig','position','reference_kmer','modification_status'] # for groupby
# eventalign_result.loc[:, 'length'] = pd.to_numeric(eventalign_result['end_idx'])-pd.to_numeric(eventalign_result['start_idx'])
# eventalign_result.loc[:, 'sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) * eventalign_result['length']
# eventalign_result.loc[:, 'sum_norm_std'] = pd.to_numeric(eventalign_result['event_stdv']) * eventalign_result['length']
# eventalign_result.loc[:, 'sum_dwell_time'] = pd.to_numeric(eventalign_result['event_length']) * eventalign_result['length']
#
# eventalign_result = eventalign_result.groupby(keys)
# # print(eventalign_result.shape)
# sum_norm_mean = eventalign_result['sum_norm_mean'].sum()
# sum_norm_std = eventalign_result["sum_norm_std"].sum()
# sum_dwell_time = eventalign_result["sum_dwell_time"].sum()
#
# start_idx = eventalign_result['start_idx'].min()
# end_idx = eventalign_result['end_idx'].max()
# total_length = eventalign_result['length'].sum()
#
# # print(eventalign_result.head())
#
#
# eventalign_result = pd.concat([start_idx, end_idx], axis=1)
# # print(eventalign_result.shape)
# eventalign_result['norm_mean'] = (sum_norm_mean / total_length).round(1)
# eventalign_result["norm_std"] = sum_norm_std / total_length
# eventalign_result["dwell_time"] = sum_dwell_time / total_length
# # print(eventalign_result.head())
# eventalign_result = eventalign_result.reset_index()
# print(eventalign_result.head())
# eventalign_result['transcript_id'] = eventalign_result['contig']    #### CHANGE MADE ####
# eventalign_result['transcriptomic_position'] = pd.to_numeric(eventalign_result['position']) + 2
# print(eventalign_result.head())# the middle position of 5-mers.