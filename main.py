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

# Output format
outputFormat = '23andMe'
#outputFormat = 'ancestry'
#outputFormat = 'FamilyTreeDNA'
#outputFormat = 'MyHeritage'
#outputFormat = 'Living DNA'


# Sorting order for company column
companyPriorityList = [ '23andMe',
                        'ancestry',
                        'FamilyTreeDNA',
                        'MyHeritage',
                        'Living DNA'
                        'MyHeritage (OLD)',
                        ]

# Input/output file directory
inputFileDir = './input/'
outputFileDir = './output/'

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)

# Output File name
outputFileName = 'DNASuperKit'
outputFileEnding = '.csv'

##########################################


##########################################
# Normalization tables
# 

# 23andMe chromosome numbering and order:       Living DNA chromosome numbering and order:
## rsid	chromosome	position	genotype        # rsid	chromosome	position	genotype
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18                11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22                                19, 20, 21, 22
# X, Y, MT                                      X
# Alleles = paired                              Alleles = paired
# Nocalls = --                                  Nocalls = not included
# Tabulated                                     Tabulated

# 23andMe Alleles                               Living DNA Alleles
# ['DD' 'II' 'DI']
# ['D' 'I']

# ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']               ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'GT' 'CG' 'AT']                         ['AC' 'GT' 'CG' 'AT']
#                                               ['CA' 'TG' 'GC' 'TA' 'TC' 'GA']
#
# ['A' 'C' 'G' 'T']                             ['A' 'C' 'G' 'T']
#
# ['--']


# BEFORE MARCH 1 2019
# MyHeritage chromosome numbering and order:    FamilyTreeDNA chromosome numbering and order:
#RSID,CHROMOSOME,POSITION,RESULT                RSID,CHROMOSOME,POSITION,RESULT
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18                11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22                                19, 20, 21, 22
# X, Y                                          XY, MT, X (XY has a overlap with the top and bottom part of 23andMe X CHROMOSOME)
# Alleles = paired                              Alleles = paired
# Nocalls = --                                  Nocalls = --
# comma                                         comma
# "rsid"

# MyHeritage Alleles                            FamilyTreeDNA Alleles
# ['AA' 'CC' 'GG' 'TT' 'AG']                    ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'CG' 'AT']                              ['AC' 'GT' 'CG' 'AT']
# ['TG' 'GC' 'TA' 'TC']
#
# ['GG' 'CC' 'AA' 'TT']                         ['-G' '-C' '-A' '-T']
#
# ['--']                                        ['--']


# ancestry chromosome numbering and order:
#
# Alleles = 
# Nocalls = 

#    Its normal. 23 is the X, 24 is Y, 25 is the PAR region, and 26 is mtDNA.

#    If the Ancestry.com data file says a SNP is from chromosome 23, it's actually from the X chromosome.
#    If it indicates chromosome 24, it's from the portion of the Y chromosome that is not part of the pseudoautosomal region.
#    If it indicates chromosome 25, the designated SNP is from the pseudoautosomal region of the Y chromosome.
#    If it indicates chromosome 26, it's mitochondrial data (which is present in at least some Ancestry data produced since May 2016).


# Normalized chromosome numbering and order:
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22
# X, Y, MT
# Alleles = paired
# Nocalls = --

# Sorting order for chromosome column
chromosomePriorityList = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' ]

# Table for normalizing chromosome names
chromosomeTable = { '23': 'X', '25': 'X', 'XY': 'X',
                    '26': 'MT',
                    '24': 'Y'
}

# Normalized allele naming convention
# ['DD' 'II' 'DI']
# ['D' 'I']
# ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'GT' 'CG' 'AT']
# ['A' 'C' 'G' 'T']
# ['--']

# Table for normalizing genotype
genotypeTable = {
    '00': '--',
    'CA': 'AC',
    'TG': 'GT',
    'GC': 'CG',
    'TA': 'AT',
    'TC': 'CT',
    'GA': 'AG',

    'ID': 'DI',

#    'A': 'AA',
#    'T': 'TT',
#    'D': 'DD',
#    'C': 'CC',

    '-G': 'G',
    '-C': 'C',
    '-A': 'A',
    '-T': 'T'
}

# Table for normalizing genotype on X Y MT
genotypeTableXYMT = {
    'GG': 'G',
    'CC': 'C',
    'AA': 'A',
    'TT': 'T'
}

# Nocall list
noCalls = [ 'DD', 'II', 'DI', 'D', 'I', '--' ]

##########################################


##########################################
# Company tables
# 

# Sorting order for FamilyTreeDNA chromosome column
chromosomePriorityListFamilyTreeDNA = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'XY', 'MT', 'X' ]

genotypeTableFamilyTreeDNA = {
    'G': '-G',
    'C': '-C',
    'A': '-A',
    'T': '-T'
}

noCallsFamilyTreeDNA = [ 'DD', 'II', 'DI', 'D', 'I' ]

# Sorting order for ancestry chromosome column
chromosomePriorityListAncestry =      [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '26' ]

# Table for converting file to ancestry format
chromosomeTableAncestry = { 'X':  '23',
                            'MT': '26',
                            'Y':  '24'
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
    #  n = number of comment lines.
    #       23andMe = 19
    #       Living DNA = 11
    #       MyHeritageBef1March2019 = 6
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
    #print( inputFile + ' contains data from 23andMe' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = '23andMe'

    return df

# LivingDNA
def prepareLivingDNA( inputFile ):
    #print( inputFile + ' contains data from Living DNA' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'Living DNA'

    return df

# MyHeritage (Before 1 March, 2019)
def prepareMyHeritageBef1March2019( inputFile ):
    #print( inputFile + ' contains data from MyHeritage (Before 1 March, 2019)' )

    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, comment='#' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'MyHeritage (OLD)'

    return df

# FamilyTreeDNA
def prepareFamilyTreeDNA(inputFile):
    # Load input file into pandas and create the proper columns
    df = pd.read_csv( inputFile, dtype=str, comment='#' )
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]
    df[ 'company' ] = 'FamilyTreeDNA'
    
    # Drop 0 value, likely nocalls or incomplete information
    df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )

    return df

# Ancestry
def prepareAncestry( inputFile ):
    #print( inputFile + ' contains data from ancestry' )

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

    # Normalize genotype with custom genotypeTable
    inputFile[ 'genotype' ].replace( to_replace=genotypeTable, inplace=True )

    # Normalize chromosome order with custom chromosomeTable
    inputFile[ 'chromosome' ].replace( to_replace=chromosomeTable, inplace=True )

    # Normalize genotypes on X, Y & MT where double char becomes single, eg. AA = A
    rep = inputFile[ 'genotype' ].replace( genotypeTableXYMT )
    inputFile[ 'genotype' ] = inputFile[ 'genotype' ].mask( inputFile[ 'chromosome' ].isin(['X','Y','MT']), rep)
    print( inputFile )


    # Drop nocalls (-, --, 00, DD, II, I, D, DI)
    #indexNames = inputFile[ inputFile[ 'genotype' ] == '--' ].index
    #inputFile.drop( indexNames, inplace = True )

    # Convert position to numerical
    inputFile[ 'position' ] = pd.to_numeric( inputFile[ 'position' ] )

    return inputFile

##########################################


##########################################
# Drop duplicates on genotype, keeping
# only genotype according to priority list
# in companyPriorityList

def dropDuplicatesDNAFile( inputFile ):

    # First drop genotype duplicates in the same position
    #inputFile.drop_duplicates( subset=[ 'chromosome','position', 'genotype'], keep='first', inplace=True )

    # Drop nocalls only if there are duplicate with calls
    # is the result not a "--"?
    m = inputFile.loc[ :, 'genotype' ].ne( '--' )
    # is there at least a non "--" in the group?
    m2 = (m
        .groupby( [ inputFile[ 'chromosome' ], inputFile[ 'position' ] ] )
        .transform( 'max' )
        )
    # perform dropping
    inputFile.loc[ m|~m2 ]

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


##########################################
# prepare database for company specific
# output format

# 23andMe format
def formatDNAFile23andMe( inputFile ):
    # Change column names according to 23andMe format
    inputFile.rename( columns = { 'rsid':'# rsid' }, inplace = True )

    return inputFile

# LivingDNA format
def formatDNAFileLivingDNA( inputFile ):
    # Drop nocalls
    for f in noCalls:
        indexNames = inputFile[ inputFile[ 'genotype' ] == f ].index
        inputFile.drop( indexNames, inplace = True )

    # Change column names according to LivingDNA format
    inputFile.rename( columns = { 'rsid':'# rsid' }, inplace = True )

    return inputFile

# MyHeritage format
def formatDNAFileMyHeritage( inputFile ):
    # Change column names according to  (Before 1 March, 2019) format
    inputFile.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )

    return inputFile

# FamilyTreeDNA format
def formatDNAFileFamilyTreeDNA( inputFile ):
    # Normalize genotype with custom genotypeTable
    inputFile[ 'genotype' ].replace( to_replace=genotypeTableFamilyTreeDNA, inplace=True )

    # Drop nocalls
    for f in noCallsFamilyTreeDNA:
        indexNames = inputFile[ inputFile[ 'genotype' ] == f ].index
        inputFile.drop( indexNames, inplace = True )

    # Custom sorting order on chromosome and company column. Modify at top of file.
    inputFile[ 'chromosome' ] = pd.Categorical( inputFile[ 'chromosome' ], chromosomePriorityListFamilyTreeDNA )

    # Sort frame based on custom sorting orders and position
    inputFile.sort_values( [ 'chromosome', 'position', ], ascending=( True, True ), inplace=True )

    # Change column names according to FamilyTreeDNA format
    inputFile.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )

    return inputFile

# Ancestry format
def formatDNAFileAncestry( inputFile ):
    # Normalize chromosome order with custom chromosomeTable
    #inputFile[ 'chromosome' ].replace( to_replace=chromosomeTableAncestry, inplace=True )

    # Custom sorting order on chromosome and company column. Modify at top of file.
    #inputFile[ 'chromosome' ] = pd.Categorical( inputFile[ 'chromosome' ], chromosomePriorityListAncestry )

    # Sort frame based on custom sorting orders and position
    #inputFile.sort_values( [ 'chromosome', 'position', ], ascending=( True, True ), inplace=True )

    # Change column names according to Ancestry format
    #inputFile.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )



    return inputFile

##########################################


####################################################################################
# MAIN LOOP
####################################################################################

print()

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
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from 23andMe' )

        DNAFile23andME = prepare23andMe( f )
        DNAFile23andME = cleanDNAFile( DNAFile23andME )

        # IF position contains genotype larger than two aleles, drop row (clean dirty information)
        DNAFile23andME = DNAFile23andME.drop( DNAFile23andME[ DNAFile23andME['genotype'].str.len() > 2 ].index )

        #print ( DNAFile23andME.genotype.unique() )

        resultFiles.append( DNAFile23andME )
        print( 'Done!' )
        print()

    # Living DNA
    elif 'Living DNA' in fileScreening:
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from Living DNA' )

        DNAFileLivingDNA = prepareLivingDNA( f )
        DNAFileLivingDNA = cleanDNAFile( DNAFileLivingDNA )

        # IF position contains genotype larger than two aleles, drop row (clean dirty information)
        DNAFileLivingDNA = DNAFileLivingDNA.drop( DNAFileLivingDNA[ DNAFileLivingDNA['genotype'].str.len() > 2 ].index )

        #print ( DNAFileLivingDNA.genotype.unique() )

        resultFiles.append( DNAFileLivingDNA )
        print( 'Done!' )
        print()

    # MyHeritage (Before 1 March, 2019)
    elif 'MyHeritage' in fileScreening:
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from MyHeritage (Before 1 March, 2019)' )

        DNAFileMyHeritageBef1March2019 = prepareMyHeritageBef1March2019( f )
        DNAFileMyHeritageBef1March2019 = cleanDNAFile( DNAFileMyHeritageBef1March2019 )

        # IF position contains genotype larger than two aleles, drop row (clean dirty information)
        DNAFileMyHeritageBef1March2019 = DNAFileMyHeritageBef1March2019.drop( DNAFileMyHeritageBef1March2019[ DNAFileMyHeritageBef1March2019['genotype'].str.len() > 2 ].index )

        #print ( DNAFileMyHeritageBef1March2019.genotype.unique() )

        resultFiles.append( DNAFileMyHeritageBef1March2019 )
        print( 'Done!' )
        print()

    # FamilyTreeDNA
    elif 'RSID,CHROMOSOME,POSITION,RESULT' in fileScreening:
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from FamilyTreeDNA' )

        DNAFileFamilyTreeDNA = prepareFamilyTreeDNA( f )
        DNAFileFamilyTreeDNA = cleanDNAFile( DNAFileFamilyTreeDNA )

        # IF position contains genotype larger than two aleles, drop row (clean dirty information)
        DNAFileFamilyTreeDNA = DNAFileFamilyTreeDNA.drop( DNAFileFamilyTreeDNA[ DNAFileFamilyTreeDNA['genotype'].str.len() > 2 ].index )

        #print ( DNAFileFamilyTreeDNA.genotype.unique() )

        resultFiles.append( DNAFileFamilyTreeDNA )
        print( 'Done!' )
        print()

    # Ancestry
    elif 'ancestry' in fileScreening:
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from Ancestry' )

        DNAFileAncestry = prepareAncestry( f )
        DNAFileAncestry = cleanDNAFile( DNAFileAncestry )
        #DNAFileAncestry = sortDNAFile( DNAFileAncestry )

        # IF position contains genotype larger than two aleles, drop row (clean dirty information)
        DNAFileAncestry = DNAFileAncestry.drop( DNAFileAncestry[ DNAFileAncestry['genotype'].str.len() > 2 ].index )

        print ( DNAFileAncestry.genotype.unique() )

        resultFiles.append( DNAFileAncestry )
        print( 'Done!' )
        print()

    # Else file is unknown
    else:
        print( 'File ' + f.replace( inputFileDir, '' ) + ' is unknown' )
        print()


# Check if there are objects in DNASuperKit
# if not, then quit script
if not resultFiles:
    print ("No compatible files has been found")

    exit()

print()
print( 'Processing results' )
print()

# Concatenate all DNA files into one list
DNASuperKit = pd.concat(resultFiles, sort=False, ignore_index=True)

# Sort DNA according to order provided in customization
DNASuperKit = sortDNAFile( DNASuperKit )

#DEBUG
#print ( DNASuperKit.genotype.unique() )

# Drop duplicates
DNASuperKit = dropDuplicatesDNAFile( DNASuperKit )

# Delete 'company' column
del DNASuperKit[ 'company' ]

#Format .csv to a specific company format
print( 'Correcting data to correspond with ' + outputFormat + ' format and saving to ' + outputFileEnding )
# 23andMe
if outputFormat == '23andMe':
    DNASuperKit = formatDNAFile23andMe( DNASuperKit )
    DNASuperKit.to_csv( outputFileDir + outputFileName + '-23andMe' + outputFileEnding, sep='\t', index=None )

# Living DNA
elif outputFormat == 'Living DNA':
    DNASuperKit = formatDNAFileLivingDNA( DNASuperKit )
    DNASuperKit.to_csv( outputFileDir + outputFileName + '-LivingDNA' + outputFileEnding, sep='\t', index=None )

# FamilyTreeDNA
elif outputFormat =='FamilyTreeDNA':
    DNASuperKit = formatDNAFileFamilyTreeDNA( DNASuperKit )
    DNASuperKit.to_csv( outputFileDir + outputFileName + '-FamilyTreeDNA' + outputFileEnding, sep='\t', index=None )

# MyHeritage
elif outputFormat == 'MyHeritage':
    DNASuperKit = formatDNAFileMyHeritage( DNASuperKit )
    DNASuperKit.to_csv( outputFileDir + outputFileName + '-MyHeritage' + outputFileEnding, sep='\t', index=None )

# Ancestry
elif outputFormat == 'Ancestry':
    print ( outputFormat )

# Success!   
print()
print( 'DNA SuperKit successfully created!' )
print()

####################################################################################
####################################################################################

# TODO
# * Add comments on top of superkit file
# * Choosable output format (23andMe, FamilyTreeDNA, LivingDNA, MyHeritage, Ancestry). Partially done, not entirely possible because of the need to normalize
# * Improve code

##########################################
# DEBUGGING

#df.info()
#print( df )

######################################
