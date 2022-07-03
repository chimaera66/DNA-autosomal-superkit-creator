##############################################################################################
# 
##############################################################################################

import os
import pandas as pd


##########################################
# Customizations
# change to fit your need
#

# Sorting order for company column
companyPriorityList = [ '23andMe', 'FamilyTreeDNA', 'Living DNA', 'MyHeritage', 'ancestry' ]
# Sorting order for chromosome column
chromosomePriorityList = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' ]

# Input file
#inputFile = 'LivingDNA2.csv'
#inputFile = '23andMe.csv'
#inputFile = 'LivingDNA.csv'
#inputFile = 'MyHeritage.csv'
inputFile = 'FamilyTreeDNA.csv'
#inputFile = 'ancestry.csv'

##########################################


##########################################
# Normalization tables
# 

# 23andMe chromosome numbering and order:       Living DNA chromosome numbering and order:
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18                11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22                                19, 20, 21, 22
# X, Y, MT                                      X
# Alleles = paired                              Alleles = paired
# Nocalls = --                                  Nocalls = not included

# MyHeritage chromosome numbering and order:    FamilyTreeDNA chromosome numbering and order:
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18                11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22                                19, 20, 21, 22
# X, Y                                          XY, MT, X
# Alleles = paired                              Alleles = paired
# Nocalls = --                                  Nocalls = --

# ancestry chromosome numbering and order:
#
# Alleles = 
# Nocalls = 

# Normalized chromosome numbering and order:
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22
# X, Y, MT
# Alleles = paired
# Nocalls = --

# Table for normalizing genotype
genotype_table = {
    '00': '--',
    'CA': 'AC',
    'TA': 'AT',
    'TC': 'CT',
    'ID': 'DI',
    'GA': 'AG',
    'GC': 'CG',
    'TG': 'GT',
    'A': 'AA',
    'T': 'TT',
    'D': 'DD',
    'C': 'CC',
}

# Table for normalizing chromosome names
chromosome_table = {
#    '23': 'X',
#    '25': 'X',
    'XY': 'X',
#    '26': 'MT',
#    '24': 'Y'
}

##########################################


##########################################
# Pre-screen file to determin DNA company
#

def prescreenDNAFile( inputDNAFile ):
    ##############################
    #  n = max number of comment lines.
    #       23andMe = 19
    #       Living DNA = 11
    #       MyHeritage = 6
    #       FamilyTree DNA = 0
    #       ancestry = ?

    n = 19 
    mystring = ' '

    with open( inputDNAFile ) as myfile:
        head = [ next( myfile ) for x in range( n ) ]
    #print(head)

    for x in head:
       mystring += ' ' + x

    return mystring

##########################################


##########################################
# Prepare DNA file differently depending
# on which DNA testing company the file
# comes from
#

# 23andMe
def prepare23andMe( inputFile ):
    print( inputFile + ' contains data from 23andMe' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = '23andMe'

    return df

# LivingDNA
def prepareLivingDNA( inputFile ):
    print( inputFile + ' contains data from Living DNA' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'Living DNA'

    return df

# MyHeritage
def prepareMyHeritage( inputFile ):
    print( inputFile + ' contains data from MyHeritage' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, comment='#' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'MyHeritage'

    return df

# FamilyTreeDNA
def prepareFamilyTreeDNA(inputFile):
    print( inputFile + ' contains data from FamilyTreeDNA' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, comment='#' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'FamilyTreeDNA'

    return df

# Ancestry
def prepareAncestry( inputFile ):
    print( inputFile + ' contains data from ancestry' )

    # Load input file into pandas and create the proper columns
    data = pd.read_csv( inputFile, dtype=str, sep='\t', comment='#' )
    data[ 'genotype' ] = data[ 'allele1' ] + data[ 'allele2' ]
    del data[ 'allele1' ]
    del data[ 'allele2' ]
    df[ 'company' ] = 'ancestry'

    return df

##########################################


##########################################
# Clean file and normalize chromosome
# and genotype

def cleanDNAFile( inputFile ):

    # Normalize genotype with custom genotype_table
    inputFile[ 'genotype' ].replace( to_replace=genotype_table, inplace=True )
    # Normalize chromosome order with custom chromosome_table
    inputFile[ 'chromosome' ].replace( to_replace=chromosome_table, inplace=True )

    # Drop nocalls
    indexNames = inputFile[ inputFile[ 'genotype' ] == '--' ].index
    inputFile.drop( indexNames, inplace = True )
    #inputFile[ inputFile.genotype != '--' ] 

    # Drop duplicates
    inputFile.drop_duplicates( subset=[ 'chromosome','position', 'genotype'], keep='first', inplace=True )

    return inputFile

##########################################


##########################################
# Sort file based on custom chromosome
# order, position and custom genotype
# order
def sortDNAFile( inputFile ):
    # Custom sorting order on chromosome column. Modify at top of file.
    #df[ 'chromosome' ] = pd.Categorical(df[ 'chromosome' ], chromosomePriorityList )
    inputFile[ 'chromosome' ] = pd.Categorical( inputFile[ 'chromosome' ], chromosomePriorityList )


    # Custom sorting order on company column. Modify at top of file.
    inputFile[ 'company' ] = pd.Categorical( inputFile[ 'company' ], companyPriorityList )
    inputFile.sort_values( [ 'chromosome', 'position', 'company' ] , ascending=(True, True, True), inplace=True )

    return inputFile


##########################################


####################################################################################
# MAIN LOOP
####################################################################################

fileScreening = prescreenDNAFile( inputFile )

resultFiles = []

if '23andMe' in fileScreening:
    DNAFile23andME = prepare23andMe( inputFile )
    DNAFile23andME = cleanDNAFile( DNAFile23andME )
    DNAFile23andME = sortDNAFile( DNAFile23andME )

    resultFiles.append( DNAFile23andME )

elif 'Living DNA' in fileScreening:
    DNAFileLivingDNA = prepareLivingDNA( inputFile )
    DNAFileLivingDNA = cleanDNAFile( DNAFileLivingDNA )
    DNAFileLivingDNA = sortDNAFile( DNAFileLivingDNA )

    resultFiles.append( DNAFileLivingDNA )

elif 'MyHeritage' in fileScreening:
    DNAFileMyHeritage = prepareMyHeritage( inputFile )
    DNAFileMyHeritage = cleanDNAFile( DNAFileMyHeritage )
    DNAFileMyHeritage = sortDNAFile( DNAFileMyHeritage )

    resultFiles.append( DNAFileMyHeritage )

elif 'RSID,CHROMOSOME,POSITION,RESULT' in fileScreening:
    DNAFileFamilyTreeDNA = prepareFamilyTreeDNA( inputFile )
    DNAFileFamilyTreeDNA = cleanDNAFile( DNAFileFamilyTreeDNA )
    DNAFileFamilyTreeDNA = sortDNAFile( DNAFileFamilyTreeDNA )

    resultFiles.append( DNAFileFamilyTreeDNA )

elif 'ancestry' in fileScreening:
    DNAFileAncestry = prepareAncestry( inputFile )
    DNAFileAncestry = cleanDNAFile( DNAFileAncestry )
    DNAFileAncestry = sortDNAFile( DNAFileAncestry )

    resultFiles.append( DNAFileAncestry )

else:
    print( 'Source file is unknown' )

for f in resultFiles:
    df = f
    #df.info()
    #print( df )

# TESTING
df = resultFiles[0]
print( df )

# clean and sort Raw data file
#df = cleanDNAFile( df )
#df = sortDNAFile( df )

# Write results to file
df.to_csv( 'test.csv', sep='\t', index=None )

####################################################################################
####################################################################################

#print( df )




##########################################
# DEBUGGING

#df.info()
#print( df )

######################################








######################################
# Help
#

# ---------------------------------------------------------------------
# Drop duplicate value and preserver first entry
#import pandas as pd

#df = pd.DataFrame()
#df.insert(loc=0,column='Column1',value=['cat',     'toy',    'cat'])
#df.insert(loc=1,column='Column2',value=['bat',    'flower',  'bat'])
#df.insert(loc=2,column='Column3',value=['xyz',     'abc',    'lmn'])

#df = df.drop_duplicates(subset=['Column1','Column2'],keep='first')
#print(df)
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# You can remove duplicates based on duplicates in the Name and Age 
# df = df.drop_duplicates(subset=['Name', 'Age'])
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Sort pandas dataframe by column
# df.sort_values(by = 'Name')
# df = df.sort_values(['column1', 'column2'], ascending=(False, True))
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Sort pandas dataframe by custom order on string index

# First make the month column a categorical and specify the ordering to use.
#df['m'] = pd.Categorical(df['m'], ['March', 'April', 'Dec'])
# Now, when you sort the month column it will sort with respect to that list:
#df.sort_values('m')

#df['Company'] = pd.Categorical(df['Company'], ['23andMe', 'Living DNA', 'MyHeritage', 'FamilyTreeDNA', 'ancestry'])
#df.sort_values('Company')
# ---------------------------------------------------------------------

#reorderlist = [
#    'Maurice Baker',
#    'Adrian Caldwell',
#    'Ratko Varda',
#    'Ryan Bowen',
#    'Cedric Hunter'
#]
