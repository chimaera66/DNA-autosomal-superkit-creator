##############################################################################################
# 
##############################################################################################

import os
import pandas as pd

file = '23andMe.csv'

#df = pd.read_csv(file, dtype=str, sep='\t', comment="#")
df = pd.read_csv(file, dtype=str, sep='\t', comment="#", index_col=False, engine='python')
df.columns = ['rsid', 'chromosome', 'position', 'result']
df['company'] = '23andME'

if '23andMe' in file:
    print("23andMe")
#    ancestry = clean_ancestry(f)
#    result_files.append(ancestry)

df.info()

print( df )


##########################################
# Read file partially and search for text
#

N = 19
with open(file) as myfile:
    head = [next(myfile) for x in range(N)]
print(head)

mystring = ' '

for x in head:
    mystring += ' ' + x

if '23andMe' in mystring:
	print("shoes are in stock.")

#print(type(head))

#test ' '.join(head)




#x = head.find("23andMe")

#print(x)

