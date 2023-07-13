import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

data_path = '/data/jiayi/IVT_results/Root/data.readcount'
df = pd.read_csv(data_path)
min_reads = 1
print(df)
read_idx = df['n_reads'] >= min_reads

df['set_type'] = 'N'
#pos=pd.DataFrame(output[output['modification_status']== 1])
#neg =pd.DataFrame(output[output['modification_status']==0])
#print(pos)
num_neg = len(df[df['modification_status']==0])
num_pos = len(df[df['modification_status']==1])
print(num_neg)


split_neg = np.random.permutation(num_neg)
df.iloc[split_neg[:int(0.8*num_neg)],df.columns.get_loc('set_type')]='Train'
#pos_Train = pos.iloc[np.array(split_pos[:int(0.8*num_pos)]),]
#pos_Train['set_type'] = pos_Train['set_type'].replace({'N': 'Train'})
#pos[split_pos[:int(0.8*num_pos)]]['set_type'] = 'Train'
#print(pos_Train)

df.iloc[split_neg[int(0.8*num_neg):(int(0.8*num_neg)+int(0.1*num_neg))],df.columns.get_loc('set_type')]='Val'
#pos_Val=pos.iloc[np.array(split_pos[int(0.8*num_pos):(int(0.8*num_pos)+int(0.1*num_pos))]),]
#pos_Val['set_type'] = pos['set_type'].replace({'N': 'Val'})
#print(pos_Val)
#pos[split_pos[int(0.8*num_pos):(int(0.8*num_pos)+int(0.1*num_pos))]]['set_type']= 'Val'



df.iloc[split_neg[int(0.9*num_neg):(int(0.9*num_neg)+int(0.1*num_neg))],df.columns.get_loc('set_type')]='Test'
#pos_Test=pos.iloc[np.array(split_pos[int(0.9*num_pos):(int(0.9*num_pos)+int(0.1*num_pos))]),]
#pos_Test['set_type'] = pos_Test['set_type'].replace({'N': 'Test'})
#pos[split_pos[int(0.9*num_pos):(int(0.9*num_pos)+int(0.1*num_pos))]]['set_type']= 'Test'
#print(pos_Test)

#pos=(pos_Train.append(pos_Val)).append(pos_Test)
#pos_frames=[pos_Train,pos_Test,pos_Val]
#pos=pd.concat(pos_frames,axis=0)


#pos_path='/data/jiayi/IVT_results/pos.csv'
#pos.to_csv(pos_path)

split_pos = np.random.permutation(num_pos)
add = list(map(lambda x : x + 1372, split_pos))
df.iloc[add[:int(0.8*num_pos)],df.columns.get_loc('set_type')]='Train'
df.iloc[add[int(0.8*num_pos):(int(0.8*num_pos)+int(0.1*num_pos))],df.columns.get_loc('set_type')]='Val'
df.iloc[add[int(0.9*num_pos):(int(0.9*num_pos)+int(0.1*num_pos))],df.columns.get_loc('set_type')]='Test'
print(df)



#df.iloc[]
#pos_Train=neg.iloc[np.array(split_neg[:int(0.8*num_neg)]),]
#pos_Train['set_type']=neg_Train['set_type'].replace({'N': 'Train'})

#for i in split_neg[int(0.8*num_neg):(int(0.8*num_neg)+int(0.1*num_neg))]:
#neg_Val=neg.iloc[np.array(split_neg[int(0.8*num_neg):(int(0.8*num_neg)+int(0.1*num_neg))])]
#neg_Val['set_type'] = neg_Val['set_type'].replace({'N': 'val'})

#for i in split_neg[int(0.9*num_neg):(int(0.9*num_neg)+int(0.1*num_neg))] :
#neg_Test=neg.iloc[np.array(split_neg[int(0.9*num_neg):(int(0.9*num_neg)+int(0.1*num_neg))])]
#neg_Test['set_type'] = neg_Test['set_type'].replace({'N': 'Test'})




#neg=(neg_Train.append(neg_Val)).append(neg_Test)
#print(neg)
#neg_path='/data/jiayi/IVT_results/neg.csv'
#neg.to_csv(neg_path)

#outcome=pd.concat([neg,pos])
# outcome=outcome.resetindex()
df_path='/data/jiayi/IVT_results/Root/data.readcount.labelled'
df.to_csv(df_path,index=False)


#print(output['modification_status'].mean())