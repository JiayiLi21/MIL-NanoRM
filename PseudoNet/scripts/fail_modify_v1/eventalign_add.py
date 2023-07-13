import pandas as pd
df_1 = pd.read_csv('/data/jiayi/IVT_results/check/Normal/eventalign.txt')
df_2 = pd.read_csv('/data/jiayi/IVT_results/check/Pseudo/eventalign.txt')
df_1.loc[:,'modification_status'] = 0
df_2.loc[:, 'modification_status'] = 1

df_frames=[df_1,df_2]
df=pd.concat(df_frames,axis=0)

#df_1 = df_1.append(df_2)
df.to_csv('/data/jiayi/IVT_results/check/Check_2*1000/eventalign.csv',index=False)

df=pd.read_csv('/data/jiayi/IVT_results/check/Check_2*1000/eventalign.csv')

with open('/data/jiayi/IVT_results/check/Check_2*1000/eventalign.txt','w') as f:
    for i in df:
        #i = i.replace(',', ' ')
        f.write(str(i))


