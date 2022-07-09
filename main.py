##############################################################################################
# 
##############################################################################################

import os
from typing import List
import pandas as pd



##########################################
# Customizations
# change to fit your need
#

# Sorting order for company column
companyPriorityList = [ '23andMe',
                        'FamilyTreeDNA',
                        'Living DNA',
                        'MyHeritage',
                        'ancestry'
                        ]

# Sorting order for chromosome column
#chromosomePriorityList = [ '0', # Only in FamilyTreeDNA, probably contains genotypes that could not be matched to a rsid.
#                           '1',
chromosomePriorityList = [ '1',
                           '2',
                           '3',
                           '4',
                           '5',
                           '6',
                           '7',
                           '8',
                           '9',
                           '10',
                           '11',
                           '12',
                           '13',
                           '14',
                           '15',
                           '16',
                           '17',
                           '18',
                           '19',
                           '20',
                           '21',
                           '22',
                           'X',
                           'Y',
                           'MT'
                        ]

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)

# Input/output file directory
inputFileDir = './input/'
outputFileDir = './output/'

outputFileName = 'DNASuperKit.csv'

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
# 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22
# X, Y, MT
# Alleles = paired
# Nocalls = --

# Most of these conversions are not really necessary
# Still makes it easier to normalize noncalls to drop them
# Table for normalizing genotype
genotype_table = {
    '00': '--',
#    'CA': 'AC',
#    'TA': 'AT',
#    'TC': 'CT',
#    'ID': 'DI',
#    'GA': 'AG',
#    'GC': 'CG',
#    'TG': 'GT',
#    'A': 'AA',
#    'T': 'TT',
#    'D': 'DD',
#    'C': 'CC',
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
# Find all files in directory with
# the desired file ending from 'fileEndings'

def findDNAFiles( inputFiles ) -> List:
    # Get script directory
    scriptDir = os.path.dirname( os.path.realpath( __file__ ) )
    # Add inputFileDir to directory to get subdir
    scriptDir = scriptDir + inputFileDir

    # List all files in subdir and add files with fileEndings
    # and append to result list
    fileList = [ f for f in os.listdir( path=scriptDir ) ]
    result = []
    for f in fileList:
        if f.lower().endswith( inputFiles ):
            result.append( inputFileDir + f )
    return result

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

    n = 0
    mystring = ' '

    # count lines in the file
    with open( inputDNAFile, 'r') as fp:
        for n, line in enumerate(fp):
            pass

    # cutoff at 19 lines
    if n > 19:
        n = 19

    # look at the first n lines
    with open( inputDNAFile ) as myfile:
        head = [ next( myfile ) for x in range( n ) ]

    # put all lines in a string
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
    
    # Drop 0 value, likely nocalls or lacking information
    df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )

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

    # Convert position to numerical
    inputFile[ 'position' ] = pd.to_numeric( inputFile[ 'position' ] )

    return inputFile

##########################################


##########################################
# Drop duplicates on genotype, keeping
# only genotype according to priority list
# in companyPriorityList

def dropDuplicatesDNAFile( inputFile ):

    # I guess this step isn't really necessary since we drop all duplicates except the first
    # according to the order in 'companyPriorityList'. With or without this creates identical data and the same nr of rows
    # First drop genotype duplicates in the same position
    #inputFile.drop_duplicates( subset=[ 'chromosome','position', 'genotype'], keep='first', inplace=True )

    # If genotype is different on the same position, then only keep the genotype from the company according to the order in companyPriorityList
    inputFile.drop_duplicates( subset=[ 'chromosome','position' ], keep='first', inplace=True )

    return inputFile

##########################################


##########################################
# Sort file based on custom chromosome order,
# position and custom genotype order

def sortDNAFile( inputFile ):
    # Custom sorting order on chromosome and company column. Modify at top of file.
    inputFile[ 'chromosome' ] = pd.Categorical( inputFile[ 'chromosome' ], chromosomePriorityList )
    inputFile[ 'company' ] = pd.Categorical( inputFile[ 'company' ], companyPriorityList )

    # Sort frame based on custom sorting orders and position
    inputFile.sort_values( [ 'chromosome', 'position', 'company' ], ascending=( True, True, True ), inplace=True )

    return inputFile


##########################################


####################################################################################
# MAIN LOOP
####################################################################################

#TEST

# Find files in dir with the correct file endings
rawDNAFiles = findDNAFiles( fileEndings )

# Check if there are any files in the directory
if not rawDNAFiles:
    print ("There is no files in the directory")

    exit()

# empty array to put results in
resultFiles = []

for f in rawDNAFiles:
    # Screening file to determin company
    fileScreening = prescreenDNAFile( f )

    # Living DNA
    if '23andMe' in fileScreening:
        DNAFile23andME = prepare23andMe( f )
        DNAFile23andME = cleanDNAFile( DNAFile23andME )
        DNAFile23andME = sortDNAFile( DNAFile23andME )

        resultFiles.append( DNAFile23andME )

    # Living DNA
    elif 'Living DNA' in fileScreening:
        DNAFileLivingDNA = prepareLivingDNA( f )
        DNAFileLivingDNA = cleanDNAFile( DNAFileLivingDNA )
        DNAFileLivingDNA = sortDNAFile( DNAFileLivingDNA )

        resultFiles.append( DNAFileLivingDNA )

    # MyHeritage
    elif 'MyHeritage' in fileScreening:
        DNAFileMyHeritage = prepareMyHeritage( f )
        DNAFileMyHeritage = cleanDNAFile( DNAFileMyHeritage )
        DNAFileMyHeritage = sortDNAFile( DNAFileMyHeritage )

        resultFiles.append( DNAFileMyHeritage )

    # FamilyTreeDNA
    elif 'RSID,CHROMOSOME,POSITION,RESULT' in fileScreening:
        DNAFileFamilyTreeDNA = prepareFamilyTreeDNA( f )
        DNAFileFamilyTreeDNA = cleanDNAFile( DNAFileFamilyTreeDNA )
        DNAFileFamilyTreeDNA = sortDNAFile( DNAFileFamilyTreeDNA )

        resultFiles.append( DNAFileFamilyTreeDNA )

    # Ancestry
    elif 'ancestry' in fileScreening:
        DNAFileAncestry = prepareAncestry( f )
        DNAFileAncestry = cleanDNAFile( DNAFileAncestry )
        DNAFileAncestry = sortDNAFile( DNAFileAncestry )

        resultFiles.append( DNAFileAncestry )

    # Else file is unknown
    else:
        print( 'File ' + f + ' is unknown' )

# Check if there are lists in DNASuperKit
# if not, then quit script
if not resultFiles:
    print ("No compatible files has been found")

    exit()

# Concatenate all DNA files into one list
DNASuperKit = pd.concat(resultFiles, sort=False, ignore_index=True)

# Sort DNA according to order provided in customization
DNASuperKit = sortDNAFile( DNASuperKit )

# Drop duplicates
DNASuperKit = dropDuplicatesDNAFile( DNASuperKit )

# Delete 'company' column
del DNASuperKit[ 'company' ]

# Write result to file
DNASuperKit.to_csv( outputFileDir + outputFileName, sep='\t', index=None )
print ("DNA SuperKit successfully created!")

####################################################################################
####################################################################################

# TODO
# * Add comments on top of superkit file
# * Choosable output format (23andMe, FamilyTreeDNA, LivingDNA, MyHeritage, Ancestry)

##########################################
# DEBUGGING

#df.info()
#print( df )

######################################
