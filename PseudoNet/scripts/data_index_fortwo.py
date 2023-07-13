origin = open('/data/jiayi/IVT_results/normalU/normalU_v2_prep_T/data.index','r')
origin2 = open('/data/jiayi/IVT_results/pseudoU/PseudoU_v2_prep_T/data_addindex','r')
new = open('/data/jiayi/IVT_results/Root/data.index','w')


# combine eventalign.index ,data.index from dataprep for normal and Pseudo
for idx, i in enumerate(origin):
    if idx == 0:
        #i = i.replace("\n",',')
        new.write(i)
        print(i)

        continue
    else:
        #i = i.replace("\n", ',')
        new.write(i)

for idx, i in enumerate(origin2):
    if idx == 0:
        continue
    else:
        #i = i.replace("\n", ",")
        new.write(i)
new.close()
origin.close()
origin2.close()