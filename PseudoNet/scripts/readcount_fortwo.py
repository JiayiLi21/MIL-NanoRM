origin = open('/data/jiayi/IVT_results/normalU/normalU_v2_prep_T/data.readcount','r')
origin2 = open('/data/jiayi/IVT_results/pseudoU/PseudoU_v2_prep_T/data.readcount','r')
new = open('/data/jiayi/IVT_results/Root/data.readcount','w')

for idx, i in enumerate(origin):
    if idx == 0:
        i = i.replace("\n", ',') + 'modification_status\n'
        new.write(i)
        print(i)

        continue
    else:
        i = i.replace("\n", ',')
        new.write(i  + "0\n")

for idx, i in enumerate(origin2):
    if idx == 0:
        continue
    else:
        i = i.replace("\n", ",")
        new.write(i + "1\n")
new.close()
origin.close()
origin2.close()