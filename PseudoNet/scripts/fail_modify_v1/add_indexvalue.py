import numpy as np
import pandas as pd

#data_path = '/data/jiayi/IVT_results/pseudoU/Net_PseudoU_results/data.index'
data_path='/data/jiayi/IVT_results/pseudoU/PseudoU_v2_prep_T/data.index'
df = pd.read_csv(data_path)
df.loc[:, 'start'] += 26912203
df.loc[:,'end']+= 26912203

df_path='/data/jiayi/IVT_results/pseudoU/PseudoU_v2_prep_T/data_addindex'
df.to_csv(df_path,index=False)
