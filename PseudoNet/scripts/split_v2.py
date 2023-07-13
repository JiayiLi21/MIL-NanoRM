import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

data_path = '/data/jiayi/IVT_results/root/data.readcount'
output = pd.read_csv(data_path)
min_reads = 1
print(output)
idx = output['n_reads'] >= min_reads
print(output)