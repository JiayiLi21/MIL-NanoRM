import argparse
import multiprocessing
import os
import warnings
from collections import defaultdict
from io import StringIO
from itertools import groupby
from operator import itemgetter

import numpy as np
import pandas as pd
import ujson

from PseudoNet.scripts import helper
from PseudoNet.scripts.constants import IVT_KMERS, NUM_NEIGHBORING_FEATURES
from PseudoNet.utils import misc

# function: get_args()
# required:
eventalign_filepath = '/data/jiayi/IVT_results/sample_select/eventalign.txt'
out_dir = '/data/jiayi/IVT_results/sample_select'
# optional:
n_processes = 1
chunk_size = 1000000
readcount_min = 20  # only transcripts with at least readcount_min reads will be reserved in processing
readcount_max = 1000
min_segment_count = 20  # only positions with at least min_segment_count reads will be reserved
skip_index = False
n_neighbors = 1
"""
Using head -n1000  and tail -n100 dataset as demo to comprehensively aquire the principle and workflow of dataprep, including detailed transformation and algorithm complexity
eventalign.txt    contig,position,reference_kmer,read_index,stran,event_index,event_level_mean,event_stdv,event_length,model_kmer,model_mean,model_stdv,standardized_level,start_idx,end_idx,modification_status
eventalign.index  transcript_id,read_index,pos_start,pos_end
data.index  transcript_id,transcript_position,start,end
data.json
data.log
data.readcount.labelled   transcript_id,transcript_position,n_reads,modification_status,set_type
"""

def index(eventalign_result,pos_start,out_paths,locks):
    eventalign_result=eventalign_result.set_index(['contig','read_index'])
    print(eventalign_result)
    pos_end=pos_start
    with locks['index'], open(out_paths['index'],'a') as f_index:
        for index in list(dict.fromkeys(eventalign_result.index)):
            transcript_id,read_index =index
            pos_end +=eventalign_result.loc[index]['line_length'].sum()
            f_index.write('%s,%d,%d,%d\n' %(transcript_id,read_index,pos_start,pos_end))
            pos_start=pos_end


def parallel_index(eventalign_filepath, chunk_size, out_dir, n_processes):
    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['index']:
        out_paths[out_filetype] = os.path.join(out_dir, 'eventalign.%s' % out_filetype)
        locks[out_filetype] = multiprocessing.Lock()
        print(locks)
    # TO DO: resume functionality for index creation

    with open(out_paths['index'], 'w') as f:
        # f.write('transcript_id,read_index,pos_start,pos_end\n') # header
        f.write('transcript_id,read_index,pos_start,pos_end\n')

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue, task_function=index, locks=locks) for i in range(n_processes)]
    for p in consumers:
        p.start()

    ## Load tasks into task_queue. A task is eventalign information of one read.
    eventalign_file = open(eventalign_filepath, 'r')
    pos_start = len(eventalign_file.readline())  # remove header
    chunk_split = None
    index_features = ['contig', 'read_index', 'line_length']
    for chunk in pd.read_csv(eventalign_filepath, chunksize=chunk_size, sep='\t'):
        chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
        chunk_concat = pd.concat([chunk_split, chunk_complete])
        chunk_concat_size = len(chunk_concat.index)
        ## read the file at where it left off because the file is opened once ##
        lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
        chunk_concat.loc[:, 'line_length'] = np.array(lines)
        task_queue.put((chunk_concat[index_features], pos_start, out_paths))
        pos_start += sum(lines)
        chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']].copy()
    ## the loop above leaves off w/o adding the last read_index to eventalign.index
    chunk_split_size = len(chunk_split.index)
    lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
    chunk_split.loc[:, 'line_length'] = np.array(lines)
    task_queue.put((chunk_split[index_features], pos_start, out_paths))

    # Put the stop task into task_queue.
    task_queue = helper.end_queue(task_queue, n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()

# parallel_preprocess_tx() --start
def parallel_preprocess_tx(eventaalign_file_path,out_dir,n_processes,readcount_min,read_count_max,n_neighbors,min_segment_count):
    out_paths,locks=dict(),dict()
    for out_filetype in ['json','index','log','readcount']:
        out_paths[out_filetype]=os.path.join(out_dir,'data.%s' %out_filetype)
        locks[out_filetype]=multiprocessing.lock()
    # write headings for files
    open(out_paths['json'],'w').close()
    with open(out_paths['index'],'w') as f:
        f.write('transcript_id,transcript_position,start,end\n')
    with open(out_paths['readcount'],'w') as f:
        f.write('transcript_id,transcript_position,n_reads,modification_status\n')
    open(out_paths['log'],'w').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [helper.Consumer(task_queue=task_queue, task_function=preprocess_tx, locks=locks) for i in
                 range(n_processes)]
    for p in consumers:
        p.start()

    df_eventalign_index=pd.read_csv(os.path.join(out_dir,'eventalign.index'))
    df_eventalign_index['transcript_id']=df_eventalign_index['transcript_id']
    tx_ids=df_eventalign_index['transcript_id'].values.tolist()
    tx_ids=list(dict.fromkeys(tx_ids))

    df_eventalign_index=df_eventalign_index.set_index('transcript_id')
    with open(eventaalign_filepath,'r') as eventalign_result:
        for tx_id in tx_ids:
            data_dict=dict()
            readcount=0
            for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                read_index,pos_start,pos_end=row['read_index'],row['pos_start'],row['pos_end']
                eventalign_result.seek(pos_start,0)
                events_str=eventalign_result.read(pos_end-pos_start)

                ## data = combine(events_str) -- start
                # If a read covers the same transcript position multiple times, they will be combined

                # convert the string to pandas dataframe
                f_string = StringIO(events_str)
                eventalign_df = pd.read_csv(f_string, delimiter='\t',
                                            names=['contig', 'position', 'reference_kmer', 'read_index', 'strand',
                                                   'event_index', 'event_level_mean', 'event_stdv', 'event_length',
                                                   'model_kmer', 'model_mean', 'model_stdv', 'standardized_level',
                                                   'start_idx', 'end_idx'])
                f_string.close()
                # ToDo: check what is the model_kmer
                cond_successfully_eventaligned = eventalign_df['reference_kmer'] == eventalign_df['model_kmer']
                if cond_successfully_eventaligned.sum() != 0:
                    eventalign_df = eventalign_df[cond_successfully_eventaligned]

                    keys = ['read_index', 'contig', 'position', 'reference_kmer']  # for groupby
                    eventalign_df.loc[:, 'length'] = pd.to_numeric(eventalign_df['end_idx']) - \
                                                     pd.to_numeric(eventalign_df['start_idx'])
                    eventalign_df.loc[:, 'sum_norm_mean'] = pd.to_numeric(eventalign_df['event_level_mean']) * \
                                                            eventalign_df['length']
                    eventalign_df.loc[:, 'sum_norm_std'] = pd.to_numeric(eventalign_df['event_stdv']) * \
                                                           eventalign_df['length']
                    eventalign_df.loc[:, 'sum_dwell_time'] = pd.to_numeric(eventalign_df['event_length']) * \
                                                             eventalign_df['length']

                    eventalign_df = eventalign_df.groupby(keys)

                    sum_norm_mean = eventalign_df['sum_norm_mean'].sum()
                    sum_norm_std = eventalign_df["sum_norm_std"].sum()
                    sum_dwell_time = eventalign_df["sum_dwell_time"].sum()

                    start_idx = eventalign_df['start_idx'].min()
                    end_idx = eventalign_df['end_idx'].max()
                    total_length = eventalign_df['length'].sum()

                    # read_index, contig, position and reference_kmer are index
                    eventalign_df = pd.concat([start_idx, end_idx], axis=1)
                    eventalign_df['norm_mean'] = (sum_norm_mean / total_length).round(1)
                    eventalign_df["norm_std"] = sum_norm_std / total_length
                    eventalign_df["dwell_time"] = sum_dwell_time / total_length
                    eventalign_df = eventalign_df.reset_index()

                    eventalign_df['transcript_id'] = eventalign_df['contig']  #### CHANGE MADE ####
                    # the middle position of 5-mers, i.e, the previous position is the start position of 5-mers
                    eventalign_df['transcriptomic_position'] = pd.to_numeric(eventalign_df['position']) + 2
                    features = ['transcript_id', 'read_index', 'transcriptomic_position', 'reference_kmer', 'norm_mean',
                                'norm_std', 'dwell_time']
                    df_events = eventalign_df[features]
                    np_events = np.rec.fromrecords(df_events, names=[*df_events])
                else:
                    np_events = np.array([])
                data = np_events
                """
                e.g., data:
                [('ENST00000361055.8', 82380, 553, 'TGGAC', 125.2, 2.364, 0.00299)
                 ('ENST00000361055.8', 82380, 554, 'GGACT', 122.5, 4.523, 0.00764)
                 ('ENST00000361055.8', 82380, 555, 'GACTC',  92.6, 2.74 , 0.00498)]
                """

                # data = combine(events_str) -- end

                if data.size > 1:
                    data_dict[read_index] = data
                readcount += 1
                # for one transcript, we store at most readcount_max number of reads
                if readcount > readcount_max:
                    break


def main():
    if not skip_index:
        parallel_index(eventalign_file_path,chunk_size,out_dir,n_processes)
    #parallel_preprocess_tx(eventalign_filepath,out_dir,n_process,readcount_min,readcount_max,n_neighbours,min_segment_count)


if __name__ =='___main__':
    main()