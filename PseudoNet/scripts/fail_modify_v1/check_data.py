import os
import json
import numpy as np
import pandas as pd
# check data ,found index problem(0) in self-preprocessed data.json

fpath = './'
# data_path = './data.json'

data_index = pd.read_csv(os.path.join(fpath, "data.index"))
read_count = pd.read_csv(os.path.join(fpath, "data.readcount.labelled"))
read_count = read_count[read_count["set_type"] == 'Train'].reset_index(drop=True)
# print(data_index)
# print(read_count)
data_info = data_index.merge(read_count, on=["transcript_id", "transcript_position"])
# print(data_info)
data_fpath = os.path.join(fpath, "data.json")
data_info = data_info[data_info["n_reads"] >= 20].reset_index(drop=True)
# print(data_info)

# for idx in np.arange(data_info.shape[0]):
#     print(idx) # find error at idx = 525
#     with open(data_fpath, 'r') as f:
#         tx_id, tx_pos, start_pos, end_pos = data_info.iloc[idx][["transcript_id", "transcript_position",
#                                                                  "start", "end"]]
#         f.seek(start_pos, 0)
#         json_str = f.read(end_pos - start_pos)
#         pos_info = json.loads(json_str)[tx_id][str(tx_pos)]
#
#         assert (len(pos_info.keys()) == 1)
#
#         kmer, features = list(pos_info.items())[0]

idx = 525
with open(data_fpath, 'r') as f:
    tx_id, tx_pos, start_pos, end_pos = data_info.iloc[idx][["transcript_id", "transcript_position",
                                                             "start", "end"]]
    print(tx_id, tx_pos, start_pos, end_pos)
    f.seek(start_pos, 0)
    json_str = f.read(end_pos - start_pos)
    print(json_str)
    pos_info = json.loads(json_str)[tx_id][str(tx_pos)]

    assert (len(pos_info.keys()) == 1)

    kmer, features = list(pos_info.items())[0]
