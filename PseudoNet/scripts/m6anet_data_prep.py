import os
import ujson
import numpy as np
import pandas as pd
from io import StringIO
from itertools import groupby
from operator import itemgetter
from collections import defaultdict
from m6anet.scripts.constants import M6A_KMERS, NUM_NEIGHBORING_FEATURES


# function: get_args()
# required:
eventalign_filepath = './demo/eventalign.txt'
out_dir = 'dataprep_results/'
# optional:
n_processes = 1
chunk_size = 1000000
readcount_min = 20  # only transcripts with at least readcount_min reads will be reserved in processing
readcount_max = 1000
min_segment_count = 20  # only positions with at least min_segment_count reads will be reserved
skip_index = False
n_neighbors = 1

"""
# parallel_index(...) -- start
# get file with transcript_id, read_index, pos start and end (txt position)

out_paths = dict()
for out_filetype in ['index']:
    # {'index': 'dataprep_results/eventalign.index'}
    out_paths[out_filetype] = os.path.join(out_dir, 'eventalign.%s' %out_filetype)

with open(out_paths['index'], 'w') as f:
    f.write('transcript_id,read_index,pos_start,pos_end\n')  # write header for eventalign.index

eventalign_file = open(eventalign_filepath, 'r')  # read eventalign.txt
pos_start = len(eventalign_file.readline())  # length of first line
chunk_split = None
index_features = ['contig', 'read_index', 'line_length']
for chunk in pd.read_csv(eventalign_filepath, chunksize=chunk_size, sep='\t'):
    chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]  # chunk except the last read
    chunk_concat = pd.concat([chunk_split, chunk_complete])
    chunk_concat_size = len(chunk_concat.index)  # number of lines
    ## read the file at where it left off because the file is opened once ##
    lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
    chunk_concat.loc[:, 'line_length'] = np.array(lines)

    # task_queue.put((chunk_concat[index_features], pos_start, out_paths)) -- start

    # def index
    eventalign_result = chunk_concat[index_features].set_index(['contig', 'read_index'])
    # MultiIndex e.g., ('ENST00000361055.8', 82380)
    pos_end = pos_start
    # pos_start and pos_end: the position in txt file and count by characters not lines
    with open(out_paths['index'], 'a') as f_index:
        for index in list(dict.fromkeys(eventalign_result.index)):
            transcript_id, read_index = index  # e.g., transcript_id = 'ENST00000361055.8', read_index = 82380
            pos_end += eventalign_result.loc[index]['line_length'].sum()  # first 3 lines in demo eventalign.txt
            f_index.write('%s,%d,%d,%d\n' % (transcript_id, read_index, pos_start, pos_end))
            pos_start = pos_end

    # task_queue.put((chunk_concat[index_features], pos_start, out_paths)) -- end

    pos_start += sum(lines)
    # leave the last read to next loop or the ending processing
    chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']].copy()

# the loop above leaves off w/o adding the last read_index to eventalign.index
chunk_split_size = len(chunk_split.index)
lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
chunk_split.loc[:, 'line_length'] = np.array(lines)

# task_queue.put((chunk_split[index_features], pos_start, out_paths)) -- start

# same as above
eventalign_result = chunk_split[index_features].set_index(['contig', 'read_index'])
pos_end = pos_start
with open(out_paths['index'], 'a') as f_index:
    for index in list(dict.fromkeys(eventalign_result.index)):
        transcript_id, read_index = index
        pos_end += eventalign_result.loc[index]['line_length'].sum()
        f_index.write('%s,%d,%d,%d\n' % (transcript_id, read_index, pos_start, pos_end))
        pos_start = pos_end

# task_queue.put((chunk_split[index_features], pos_start, out_paths)) -- end

# parallel_index(...) -- end
"""

# parallel_preprocess_tx(...) -- start

out_paths = dict()
for out_filetype in ['json', 'index', 'log', 'readcount']:
    # data.json, data.index, data.log, data.readcount
    out_paths[out_filetype] = os.path.join(out_dir, 'data.%s' % out_filetype)

# write headers
open(out_paths['json'], 'w').close()
with open(out_paths['index'], 'w') as f:
    f.write('transcript_id,transcript_position,start,end\n')
with open(out_paths['readcount'], 'w') as f:
    f.write('transcript_id,transcript_position,n_reads\n')
open(out_paths['log'], 'w').close()

df_eventalign_index = pd.read_csv(os.path.join(out_dir, 'eventalign.index'))
tx_ids = df_eventalign_index['transcript_id'].values.tolist()
tx_ids = list(dict.fromkeys(tx_ids))  # ['ENST00000361055.8', ... ]
df_eventalign_index = df_eventalign_index.set_index('transcript_id')

with open(eventalign_filepath, 'r') as eventalign_result:
    for tx_id in tx_ids:
        data_dict = dict()
        readcount = 0
        # loop through different reads from same transcript
        # each read may correspond to multiple rows in eventalign.txt, covering multiple position on that transcript
        for _, row in df_eventalign_index.loc[[tx_id]].iterrows():
            read_index, pos_start, pos_end = row['read_index'], row['pos_start'], row['pos_end']
            # online reads the lines related to tx_id - read_index
            eventalign_result.seek(pos_start, 0)
            events_str = eventalign_result.read(pos_end - pos_start)

            # data = combine(events_str) -- start
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

        if readcount >= readcount_min:
            features_arrays = []
            reference_kmer_arrays = []
            transcriptomic_positions_arrays = []

            for _, events_per_read in data_dict.items():
                # events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS) -- start

                # step 1 events = partition_into_continuous_positions(events, window_size=1)

                window_size = 1
                arr = events_per_read
                arr = arr[np.argsort(arr["transcriptomic_position"])]
                float_features = ['dwell_time', 'norm_std', 'norm_mean']
                float_dtypes = [('norm_mean', '<f8'), ('norm_std', '<f8'), ('dwell_time', '<f8')]

                float_arr = arr[float_features].astype(float_dtypes).view('<f8').reshape(-1, 3)
                kmer_arr = arr["reference_kmer"].reshape(-1, 1)
                tx_pos_arr = arr["transcriptomic_position"]
                tx_id_arr = arr["transcript_id"]

                """
                enumerate(tx_pos_arr) will index transcript position like [(0, 553), (1, 554), (2, 555)]
                groupby using funciton lambda x: x[0] - x[1] will do subtraction between the two values of each element
                e.g., 0-533=-553, 1-554=-553, 2-555=-553, resulting in the same key value -553 for all three elements
                Thus these three continuous positions will be grouped into one partition [(0, 553), (1, 554), (2, 555)]
                Finally, list(map(itemgetter(0), g)) will take the first value in each element, that is [0, 1, 2]
                
                tx_pos_arr = [553 554 555] ==> partitions = [[0, 1, 2]] i.e., index of the continuous positions
                """

                partitions = [list(map(itemgetter(0), g)) for k, g in groupby(enumerate(tx_pos_arr),
                                                                              lambda x: x[0] - x[1])]

                # only partitions with at least 3 positions reserved
                step1_output = [(float_arr[partition], kmer_arr[partition], tx_id_arr[partition], tx_pos_arr[partition])
                                for partition in partitions if len(partition) >= 2 * window_size + 1]

                # step 1 -- end

                # step 2 events = filter_partitions(events, window_size, kmers)

                def _roll(to_roll, window_size=1):
                    nex = np.concatenate([np.roll(to_roll, i, axis=0) for i in range(-1, - window_size - 1, -1)],
                                         axis=1)
                    prev = np.concatenate([np.roll(to_roll, i, axis=0) for i in range(window_size, 0, -1)], axis=1)
                    return np.concatenate((prev, to_roll, nex), axis=1)[window_size: -window_size, :]

                def _create_features(partition, window_size=1):
                    float_arr, kmer_arr, tx_id_arr, tx_pos_arr = partition
                    return _roll(float_arr, window_size), _roll(kmer_arr, window_size), \
                        tx_id_arr[window_size: -window_size], tx_pos_arr[window_size: -window_size]

                """
                the code here might look a bit complicated, 
                but all the author wants is to cut the partition into preset neighbors
                i.e., window_size = 1 => range(-1, - window_size - 1, -1) = -1 & range(window_size, 0, -1) = 1
                if you have a partition [1, 2, 3], nex = [2, 3, 1] and prev [3, 1, 2]
                np.concatenate((prev, to_roll, nex) = [[3, 1, 2], [1, 2, 3], [2, 3, 1]]
                and the index here, [window_size: -window_size, :] => [1: -1, :],
                which means from the row 1 to the last row (-1) but not include the last row, and all columns reserved
                that is, np.array([[3, 1, 2], [1, 2, 3], [2, 3, 1]])[1, :] = [1, 2, 3]
                
                What if more complicated? Let's try partition = [1, 2, 3, 4, 5] and again window_size = 1,
                nex = [2, 3, 4, 5, 1] and prev = [5, 1, 2, 3, 4]
                After concatenate, we have [[5, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 1]].
                You see, what we want is the middle three, which actually is the index [1:-1, :].
                Start from row 1 means not include the row 0, and end with row -1 means not include the last row.
                Therefore, two boundaries [5, 1, 2] and [4, 5, 1], which we also do not want, were removed.
                We have [[1, 2, 3], [2, 3, 4], [3, 4, 5]] => cutting the partition into groups of three (continously)
                """

                window_size = n_neighbors
                windowed_partition = [_create_features(partition, window_size) for partition in step1_output]

                def _filter_by_kmer(partition, kmers, window_size):
                    feature_arr, kmer_arr, tx_id_arr, tx_pos_arr = partition
                    kmers_5 = kmer_arr[:, (2 * window_size + 1) // 2]  # take the middle one of the three kmers
                    mask = np.isin(kmers_5, kmers)  # check if the middle one in the targeted kmers
                    filtered_feature_arr = feature_arr[mask, :]
                    filtered_kmer_arr = kmer_arr[mask, :]
                    filtered_tx_pos_arr = tx_pos_arr[mask]
                    filtered_tx_id_arr = tx_id_arr[mask]

                    if len(filtered_kmer_arr) == 0:
                        return []
                    else:
                        return filtered_feature_arr, filtered_kmer_arr, filtered_tx_id_arr, filtered_tx_pos_arr

                filtered_by_kmers = [_filter_by_kmer(partition, M6A_KMERS, window_size)
                                     for partition in windowed_partition]
                final_partitions = [x for x in filtered_by_kmers if len(x) > 0]

                """
                example element in final_partitions:
                (array([[2.990e-03, 2.364e+00, 1.252e+02, 7.640e-03, 4.523e+00, 1.225e+02, 
                         4.980e-03, 2.740e+00, 9.260e+01]]), 
                 array([['TGGAC', 'GGACT', 'GACTC']], dtype='<U5'), 
                 array(['ENST00000361055.8'], dtype='<U17'), 
                 array([554]))
                """

                events_per_read = step2_output = final_partitions

                # step 2 -- end

                # events_per_read = filter_events(events_per_read, n_neighbors, M6A_KMERS) -- end

                def _combine_sequence(kmers):
                    kmer = kmers[0]
                    for _kmer in kmers[1:]:
                        kmer += _kmer[-1]
                    return kmer

                """
                kmers = ['TGGAC', 'GGACT', 'GACTC'] => _combine_sequence(kmers) = ['TGGACTC']
                """

                for event_per_read in events_per_read:
                    features_arrays.append(event_per_read[0])
                    reference_kmer_arrays.append([_combine_sequence(kmer) for kmer in event_per_read[1]])
                    transcriptomic_positions_arrays.append(event_per_read[3])

            if len(features_arrays) == 0:
                continue  # end the processing of this transcript, start the processing for the next one
            else:
                features_arrays = np.concatenate(features_arrays)
                reference_kmer_arrays = np.concatenate(reference_kmer_arrays)
                transcriptomic_positions_arrays = np.concatenate(transcriptomic_positions_arrays)
                assert (len(features_arrays) == len(reference_kmer_arrays) == len(transcriptomic_positions_arrays))

            # So far, we looped through the data_dict whose keys are read index,
            # which means we process the events read by read.
            # Since different reads can cover the same position, we then sort and split by transcript position
            idx_sorted = np.argsort(transcriptomic_positions_arrays)
            positions, index = np.unique(transcriptomic_positions_arrays[idx_sorted], return_index=True, axis=0)
            # index here is the start position of next transcript position
            features_arrays = np.split(features_arrays[idx_sorted], index[1:])
            reference_kmer_arrays = np.split(reference_kmer_arrays[idx_sorted], index[1:])

            # Prepare
            # print('Reformating the data for each genomic position ...')
            data = defaultdict(dict)

            for position, features_array, reference_kmer_array in zip(positions, features_arrays,
                                                                      reference_kmer_arrays):
                kmer = set(reference_kmer_array)
                assert (len(kmer) == 1)
                if (len(set(reference_kmer_array)) == 1) and ('XXXXX' in set(reference_kmer_array)) or (
                        len(features_array) == 0):
                    continue
                if len(features_array) >= min_segment_count:
                    data[int(position)] = {kmer.pop(): features_array.tolist()}

            # write to file.
            log_str = '%s: Data preparation ... Done.' % (tx_id)
            with open(out_paths['json'], 'a') as f, \
                    open(out_paths['index'], 'a') as g, \
                    open(out_paths['readcount'], 'a') as h:

                for pos, dat in data.items():
                    """
                    Example:
                    pos: 549, dat: {'TGGACTG': [[feature_from_read_1], ..., [feature_from_read_n]]
                    pos: 554, dat: {'TGGACTG': [[feature_from_read_1], ..., [feature_from_read_n]]
                    """
                    pos_start = f.tell()  # returns the current position of the file read/write pointer within the file.
                    f.write('{')
                    f.write('"%s":{"%d":' % (tx_id, pos))  # write transcript ID and position
                    ujson.dump(dat, f)  # write kmer and features
                    f.write('}}\n')
                    pos_end = f.tell()  # returns the position after writing all the information above

                    # data.index store the txt position of each tx_id-pos in data.json
                    g.write('%s,%d,%d,%d\n' % (tx_id, pos, pos_start, pos_end))

                    # count the number of reads/features for each tx_id-pos and store it in data.readcount
                    n_reads = 0
                    for kmer, features in dat.items():
                        n_reads += len(features)
                    h.write('%s,%d,%d\n' % (tx_id, pos, n_reads))

            # Finally, write the log info
            with open(out_paths['log'], 'a') as f:
                f.write(log_str + '\n')

# parallel_preprocess_tx(...) -- end

# dataprep -- end