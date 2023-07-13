# add modification_status (0 for normalU results, 1 for PseudoU results)
# Combine two files into one
#Shuffle and distribute rows , add set_type(Train, test, Val)


import re
import pandas as pd
# https://stackoverflow.com/questions/32120949/converting-text-files-to-pandas-dataframe
# open  the file and seperate every line like below:
df = open('data.readcount.txt', "r")
lines = df.readlines()
df.close()
# remove /n at the end of each line
for index, line in enumerate(lines):
      lines[index] = line.strip()
#creating a dataframe(consider u want to convert your data to 2 columns)

df = pd.DataFrame(columns=('first_col', 'second_col'))



df['modification_status'] = df.apply(0, axis = 1) # for normalU
df['modification_status'] = df.apply(1, axis = 1) #for PseudoU
