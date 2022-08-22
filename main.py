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
outputFormat = 'SuperKit'
#outputFormat = '23andMe V5'
#outputFormat = 'Ancestry'
#outputFormat = 'FamilyTreeDNA'
#outputFormat = 'MyHeritage (Old)'
#outputFormat = 'LivingDNA v1.0.2'


# Sorting order for company column
companyPriorityList = [ '23andMe V5',
                        '23andMe',
                        'Ancestry',
                        'FamilyTreeDNA',
                        'MyHeritage',
                        'LivingDNA v1.0.2',
                        'LivingDNA',
                        'MyHeritage (Old)',
                        ]

# Input/output file directory
inputFileDir = './input/'
outputFileDir = './output/'
# Output File name and file ending
outputFileName = 'DNASuperKit'
outputFileEnding = '.csv'

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)

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
# X, Y                                          XY, MT, X (XY has a overlap with the top and bottom part of 23andMe X CHROMOSOME) PAR
#                                               XY = position 153977 - 2697868 & 8503715 - 155234707
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


# Ancestry chromosome numbering and order:
# rsid	chromosome	position	allele1	allele2
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18,
# 19, 20, 21, 22
# 23, 24, 25, 26
#
# Alleles = not paired
# Nocalls = 00
# Tabulated

# Ancestry Alleles
# ['DD' 'II' 'DI' 'ID']
#
# ['AA' 'CC' 'GG' 'TT' 'AG']
# ['AC' 'CG' 'AT']
# ['CA' 'TG' 'GC' 'TA' 'TC' 'GA']
#
# ['GG' 'CC' 'AA' 'TT']
#
# ['00']

#    Its normal. 23 is the X, 24 is Y, 25 is the PAR region, and 26 is mtDNA.

#    If the Ancestry.com data file says a SNP is from chromosome 23, it's actually from the X chromosome.
#    If it indicates chromosome 24, it's from the portion of the Y chromosome that is not part of the pseudoautosomal region.
#    If it indicates chromosome 25, the designated SNP is from the pseudoautosomal region of the Y chromosome. (PAR/XY)
#    If it indicates chromosome 26, it's mitochondrial data (which is present in at least some Ancestry data produced since May 2016).


# Normalized chromosome numbering and order:
# 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22
# X, Y, XY, MT
# Alleles = paired
# Nocalls = --

# Sorting order for chromosome column
chromosomePriorityList = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'XY', 'MT' ]

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
noCallsHyphen = [ 'DD', 'II', 'DI', 'D', 'I' ]

##########################################


##########################################
# Company tables
# 

# ====================
chromosomePriorityList23andMe = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' ]
# ====================


# ====================
# Sorting order for FamilyTreeDNA chromosome column
chromosomePriorityListFamilyTreeDNA = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'XY', 'MT', 'X' ]

# FamilyTreeDNA genotypes for Y, X, MT
genotypeTableFamilyTreeDNA = {
    'G': '-G',
    'C': '-C',
    'A': '-A',
    'T': '-T'
}
# ====================


# ====================
# Sorting order for MyHeritage (Before 1 March, 2019) chromosome column
chromosomePriorityListMyHeritage = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y' ]

# MyHeritage (Before 1 March, 2019) genotypes for Y, X, MT
genotypeTableMyHeritage = {
    'G': 'GG',
    'C': 'CC',
    'A': 'AA',
    'T': 'TT'
}
# ====================


# ====================

# ====================
# Sorting order for LivingDNA chromosome column
chromosomePriorityListLivingDNA = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X' ]
# ====================

# ====================
# Sorting order for ancestry chromosome column
chromosomePriorityListAncestry =  [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '26' ]

# Table for normalizing chromosome names
chromosomeTableAncestryIn = { '23': 'X',
                              '24': 'Y',
                              '25': 'XY',
                              '26': 'MT',
}

# Table for converting file to ancestry format
chromosomeTableAncestryOut = { 'X':  '23',
                               'Y':  '24',
                               'XY': '25',
                               'MT': '26',
}

# Ancestry genotypes for 00
genotypeTableAncestry = {
    '--': '00'
}

# ====================

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
    #       ancestry = 18

    #n = 0
    n = 1
    mystring = ' '

    # count lines in the file
    with open( inputDNAFile, 'r') as fp:
        for n, line in enumerate(fp):
            pass

    # cutoff at 19 lines
#    if n > 19:
#        n = 19

    # look at the first n lines
    with open( inputDNAFile ) as myfile:
        head = [ next( myfile ) for x in range( n ) ]

    # put all lines in a string
    for x in head:
       mystring += ' ' + x

    return mystring

##########################################


##########################################
# Try to determin what DNA testing
# company the file originates from
#

def determinDNACompany( t, f ):

    #print( f )
    #print( t )
    if '_v5_Full_' in f:
        company = '23andMe V5'

    elif '23andMe' in t:
        company = '23andMe'

    elif '# Living DNA customer genotype data download file version: 1.0.2' in t:
        company = 'LivingDNA v1.0.2'

    elif 'Living DNA' in t:
        company = 'LivingDNA'

#    elif 'MyHeritage' in t:
    elif '# MyHeritage DNA raw data.' in t:
        company = 'MyHeritage (Old)'

    elif 'RSID,CHROMOSOME,POSITION,RESULT' in t:
        company = 'FamilyTreeDNA'

    elif 'AncestryDNA' in t:
        company = 'Ancestry'

    else:
        company = 'unknown'

    return company

##########################################


##########################################
# Prepare DNA file differently depending
# on which DNA testing company the file
# comes from
#

def prepareDNAFile( f, c ):

#    print ( c )
    # 23andMe and Living DNA
    if c == '23andMe V5' or c == '23andMe' or c == 'LivingDNA' or c == 'LivingDNA v1.0.2':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )

    # MyHeritage (Before 1 March, 2019) and FamilyTreeDNA
    elif c == 'MyHeritage (Old)' or c == 'FamilyTreeDNA':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, comment='#' )

    # Ancestry
    elif c == 'Ancestry':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, sep='\t', comment='#' )

        # Normalize chromosome order with custom chromosomeTable
        df[ 'chromosome' ].replace( to_replace=chromosomeTableAncestryIn, inplace=True )

        # Merge allele1 and allele2 to genotype column
        df[ 'genotype' ] = df[ 'allele1' ] + df[ 'allele2' ]
        del df[ 'allele1' ]
        del df[ 'allele2' ]

    # Normalize columnnames
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    # and add company column
    df[ 'company' ] = c

# DEBUG
#    print ( c )
#    print ( df )
    print()
    print( 'Unique chromosomes:' )
    print( df.chromosome.unique() )
    print()
    print( 'Unique genotypes:')
    print( df.genotype.unique() )
    print()

    return df


##########################################


##########################################
# Clean file and normalize chromosome
# and genotype

def cleanDNAFile( df, c):

    # Normalize genotype with custom genotypeTable
    df[ 'genotype' ].replace( to_replace=genotypeTable, inplace=True )

    # Normalize chromosome order with custom chromosomeTable
    #df[ 'chromosome' ].replace( to_replace=chromosomeTable, inplace=True )

    # Normalize genotypes on X, Y & MT where double char becomes single, eg. AA = A
    rep = df[ 'genotype' ].replace( genotypeTableXYMT )
    df[ 'genotype' ] = df[ 'genotype' ].mask( df[ 'chromosome' ].isin( [ 'X','Y','XY','MT' ] ), rep )

    # Drop nocalls (-, --, 00, DD, II, I, D, DI)
    #indexNames = df[ df[ 'genotype' ] == '--' ].index
    #df.drop( indexNames, inplace = True )

    # Convert position to numerical
    df[ 'position' ] = pd.to_numeric( df[ 'position' ] )

    # FamilyTreeDNA, additional preparations
    if c == 'FamilyTreeDNA':
        # Drop chromosome 0 rows, likely nocalls or incomplete information
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )

    # IF position contains genotype larger than two aleles, drop row (clean dirty information from LivingDNA and more?)
    df = df.drop( df[ df['genotype'].str.len() > 2 ].index )

    return df

##########################################


##########################################
# Drop duplicates on genotype, keeping
# only genotype according to priority list
# in companyPriorityList

def dropDuplicatesDNAFile( df ):

    # First drop genotype duplicates in the same position
    #df.drop_duplicates( subset=[ 'chromosome','position', 'genotype'], keep='first', inplace=True )

    # Drop nocalls only if there are duplicate with calls
    # is the result not a "--"?
    m = df.loc[ :, 'genotype' ].ne( '--' )
    # is there at least a non "--" in the group?
    m2 = (m
        .groupby( [ df[ 'chromosome' ], df[ 'position' ] ] )
        .transform( 'max' )
        )
    # perform dropping
    df.loc[ m|~m2 ]

    # If genotype is different on the same position, then only keep the genotype from the company according to the order in companyPriorityList
    df.drop_duplicates( subset=[ 'chromosome','position' ], keep='first', inplace=True )

    return df

##########################################


##########################################
# Sort file based on custom chromosome order,
# position and custom genotype order

def sortDNAFile( df ):
    # Custom sorting order on chromosome and company column. Modify at top of file.
    df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList )
    df[ 'company' ] = pd.Categorical( df[ 'company' ], companyPriorityList )

    # Sort frame based on custom sorting orders and position
    df.sort_values( [ 'chromosome', 'position', 'company' ], ascending=( True, True, True ), inplace=True )

    return df

##########################################


##########################################
# prepare database for company specific
# output format

def formatDNAFile( df, c ):

    # Extra operations for 23andMe format
    if c == '23andMe V5':
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )

        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList23andMe )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )


    # Extra operations for FamilyTreeDNA format
    if c == 'FamilyTreeDNA':
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == 'Y' ].index )

        #Replace genotypes for X, Y, MT( A = -A )
        df[ 'genotype' ].replace( to_replace=genotypeTableFamilyTreeDNA, inplace=True )

        # Drop nocalls according to noCallsHyphen
        for f in noCallsHyphen:
            indexNames = df[ df[ 'genotype' ] == f ].index
            df.drop( indexNames, inplace = True )

        # Concat dataframe with previously dropped chromosome 0
        df = pd.concat( [df, chromosomeZero] , sort=False, ignore_index=True)

        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListFamilyTreeDNA )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )


    # Extra operations for MyHeritage (Before 1 March, 2019) format
    if c == 'MyHeritage (Old)':
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'MT' ].index )

        #Replace genotypes for X, Y, MT( A = AA )
        df[ 'genotype' ].replace( to_replace=genotypeTableMyHeritage, inplace=True )

        # Drop nocalls according to noCallsHyphen
        for f in noCallsHyphen:
            indexNames = df[ df[ 'genotype' ] == f ].index
            df.drop( indexNames, inplace = True )

        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListMyHeritage )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )


    # Extra operations for Living DNA format
    if c == 'LivingDNA v1.0.2':
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'MT' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'Y' ].index )
        
        # Drop nocalls according to noCalls
        for f in noCalls:
            indexNames = df[ df[ 'genotype' ] == f ].index
            df.drop( indexNames, inplace = True )

        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListLivingDNA )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

    # Extra operations for SuperKit format
    if c == 'SuperKit':
        # Concat dataframe with previously dropped chromosome 0
        df = pd.concat( [df, chromosomeZero] , sort=False, ignore_index=True)

        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

    # 23andMe and Living DNA
    if c == '23andMe V5' or c == 'LivingDNA v1.0.2':
        df.rename( columns = { 'rsid':'# rsid' }, inplace = True )

    # FamilyTreeDNA and MyHeritage (Before 1 March, 2019)
    elif c == 'FamilyTreeDNA' or c == 'MyHeritage (Old)':
        # Change column names according to  (Before 1 March, 2019) format
        df.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )

    elif c == 'Ancestry':
        # Normalize chromosome order with custom chromosomeTable
        df[ 'chromosome' ].replace( to_replace=chromosomeTableAncestryOut, inplace=True )

        # Custom sorting order on chromosome and company column. Modify at top of file.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListAncestry )

        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position', ], ascending=( True, True ), inplace=True )

        # Changed nocalls from -- to 00
        df[ 'genotype' ].replace( to_replace=genotypeTableAncestry, inplace=True )

        # Split genotype into allele1 and allele2
        df[ 'allele1' ] = df[ 'genotype' ].str[ :1 ]
        df[ 'allele2' ] = df[ 'genotype' ].str[ -1: ]
        del df[ 'genotype' ]

        # Change column names according to Ancestry format
        #df.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )


    return df

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


########################
# Preparing and cleaning DNA files

# empty array to put results in
resultFiles = []
chromosomeZero = pd.DataFrame()

for f in rawDNAFiles:
    # Screening file to determin company
    fileScreening = prescreenDNAFile( f )

    # Get DNA company from file comment
    company = determinDNACompany( fileScreening , f)
#    print( company )

    if company != 'unknown':
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from ' + company )
        
        # Normalize and clean DNA file in preparation for concatenation
        df = prepareDNAFile( f, company )

        if company == 'FamilyTreeDNA':
            chromosomeZero = df.loc[ df[ 'chromosome' ] == '0' ]
            del chromosomeZero[ 'company' ]
#            print( chromosomeZero )

        df = cleanDNAFile( df, company )
        
        # Append dataframe
        resultFiles.append( df )

        # Let user know processing is completed successfully
        print( 'Done!' )
        print()

    # If file is unknown
    else:
        print( 'File ' + f.replace( inputFileDir, '' ) + ' is unknown' )
        print()

# Check if there are objects in DNASuperKit
# if not, then quit script
if not resultFiles:
    print ("No compatible files has been found")

    exit()

########################


########################
# Concat and remove duplicates

print()
print( 'Processing results' )
print()

# Concatenate all DNA files into one list
DNASuperKit = pd.concat(resultFiles, sort=False, ignore_index=True)

print()
print( 'Sorting results' )
print()

# Sort DNA according to order provided in customization
DNASuperKit = sortDNAFile( DNASuperKit )

print()
print( 'Dropping duplicates' )
print()

# Drop duplicates
DNASuperKit = dropDuplicatesDNAFile( DNASuperKit )

# DEBUG
print()
print( 'Unique chromosomes:' )
print( DNASuperKit.chromosome.unique() )
print()
print( 'Unique genotypes:')
print( DNASuperKit.genotype.unique() )
print()
print( 'Total count for each company:' )
for f in companyPriorityList:
    print( f + ': ' + DNASuperKit['company'].value_counts()[f].astype(str))
print()

# Delete 'company' column
del DNASuperKit[ 'company' ]

########################


########################
# 

#Format .csv to a specific company format


DNASuperKit = formatDNAFile( DNASuperKit, outputFormat )

if outputFormat == 'SuperKit' or outputFormat == '23andMe V5' or outputFormat == 'LivingDNA v1.0.2' or outputFormat == 'Ancestry':
    outputFileEnding = '.txt'
    print( 'Correcting data to correspond with ' + outputFormat + ' format and saving to ' + outputFileEnding )
    DNASuperKit.to_csv( outputFileDir + outputFileName + '-' + outputFormat + outputFileEnding, sep='\t', index=None )
elif outputFormat == 'FamilyTreeDNA' or outputFormat == 'MyHeritage (Old)':
    outputFileEnding = '.csv'
    print( 'Correcting data to correspond with ' + outputFormat + ' format and saving to ' + outputFileEnding )
    if outputFormat == 'MyHeritage (Old)':
        DNASuperKit.to_csv( outputFileDir + outputFileName + '-' + outputFormat + outputFileEnding, index=None, quoting=2 )
    else:
        DNASuperKit.to_csv( outputFileDir + outputFileName + '-' + outputFormat + outputFileEnding, index=None )

# Success!   
print()
print( 'DNA SuperKit successfully created!' )
print()

####################################################################################
####################################################################################

# TODO
# * Add comments on top of superkit file
# * Improve company detection "algorithm"
# * Improve genotype count output

##########################################
# DEBUGGING

#df.info()
#print( df )

######################################
