##############################################################################################
# Create DNA superkit
#

import os                   # For findDNAFiles
from typing import List
import pandas as pd
import re                   # For determineDNACompany

import argparse             # Command line argument parser
import sys                  # sys.exit(1)
import time
import datetime             # Get time

import random


####################################################################################
# VARIABLES
####################################################################################

# Input/output file directory
inputFileDir = './input/'
outputFileDir = './output/'

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)

# Output File name and file ending
outputFileName = 'DNASuperKit'
outputFileEnding = '.csv'

# Sorting order for company column
companyPriorityList = [ '23andMe v5',
                        'AncestryDNA v2',
                        'FamilyTreeDNA v3',
                        'MyHeritage v2',
                        'LivingDNA v1.0.2',
                        'MyHeritage v1',
                        ]

# Trim SNPs to be within output format range
trimSNP = False
# Set majority vote off as default
majorityVote = False

# Kit statistics
superkitUniqueChromosomes = []
superkitUniqueGenotypes = []


####################################################################################
####################################################################################



####################################################################################
# COMMAND LINE ARGUMENT PARSER
####################################################################################

# Parser arguments
parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter )
parser.add_argument('-o', '--outputFormat', type=str, required=False,
                    help='''
                    Sets the template for the resulting DNA file. Valid formats:
                    SuperKit (Default)
                    23andMe v5
                    AncestryDNA v2
                    FamilyTreeDNA v3
                    LivingDNA v1.0.2
                    MyHeritage v1
                    MyHeritage v2
                    ''')
parser.add_argument('-mv', '--majorityVote', action='store_true', help='Drops duplicate genotype based on a majority vote. Considerably slower than regular keep first row drop. Only resonable if you want to merge 3 kits or more.', required=False)
parser.add_argument('-t', '--trimSNP', action='store_true', help='Trims the SNPs to fit within the tested ranges of the different companys that the outputFormat tries to emulate.', required=False)

# Get arguments from command line
args = parser.parse_args()

# Save the arguments to variables
outputFormat = args.outputFormat
trimSNP = args.trimSNP
majorityVote = args.majorityVote

# If outputFormat are not set, default to SuperKit
if not outputFormat:
    outputFormat = 'SuperKit'

# Allowed outputFormats
allowed_outputFormats = ['SuperKit', '23andMe v5', 'AncestryDNA v2', 'FamilyTreeDNA v3', 'LivingDNA v1.0.2', 'MyHeritage v1', 'MyHeritage v2']

# Check if outputFormat are valid, if not then exit
if outputFormat and outputFormat not in allowed_outputFormats:
    print(f'Invalid output format: {outputFormat}. Allowed formats are: {", ".join(allowed_outputFormats)}.')
    sys.exit(1)

####################################################################################
####################################################################################



####################################################################################
# DNA file information
####################################################################################

# 23andMe v5 chromosome numbering and order:    Living DNA v1.0.2 chromosome numbering and order:
## rsid	chromosome	position	genotype        # rsid	chromosome	position	genotype
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18                11, 12, 13, 14, 15, 16, 17, 18
# 19, 20, 21, 22                                19, 20, 21, 22
# X, Y, MT                                      X
# Alleles = paired                              Alleles = paired
# Nocalls = --                                  Nocalls = not included
# Tabulated                                     Tabulated

# 23andMe v5 Alleles                            Living DNA v1.0.2 Alleles
# ['DD' 'II' 'DI']
# ['D' 'I']

# ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']               ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'GT' 'CG' 'AT']                         ['AC' 'GT' 'CG' 'AT']
#                                               ['CA' 'TG' 'GC' 'TA' 'TC' 'GA']
#
# ['A' 'C' 'G' 'T']                             ['A' 'C' 'G' 'T']
#
# ['--']


# MyHeritage v1 chrom numbering and order:      FamilyTreeDNA v3 chromosome numbering and order:
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

# MyHeritage v1 Alleles                         FamilyTreeDNA v3 Alleles
# ['AA' 'CC' 'GG' 'TT' 'AG']                    ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'CG' 'AT']                              ['AC' 'GT' 'CG' 'AT']
# ['TG' 'GC' 'TA' 'TC']
#
# ['GG' 'CC' 'AA' 'TT']                         ['-G' '-C' '-A' '-T']
#
# ['--']                                        ['--']


# AncestryDNA v2 chromosome numbering and order:
# rsid	chromosome	position	allele1	allele2
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
# 11, 12, 13, 14, 15, 16, 17, 18,
# 19, 20, 21, 22
# 23, 24, 25, 26
#
# Alleles = not paired
# Nocalls = 00
# Tabulated

# AncestryDNA v2 Alleles
# ['DD' 'II' 'DI' 'ID']
#
# ['AA' 'CC' 'GG' 'TT' 'AG']
# ['AC' 'CG' 'AT']
# ['CA' 'TG' 'GC' 'TA' 'TC' 'GA']
#
# ['GG' 'CC' 'AA' 'TT']
#
# ['00']

#    23 is X, 24 is Y, 25 is the PAR/XY region, and 26 is mtDNA.

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

# Normalized allele naming convention
# ['DD' 'II' 'DI']
# ['D' 'I']
# ['AA' 'CC' 'GG' 'TT' 'CT' 'AG']
# ['AC' 'GT' 'CG' 'AT']
# ['A' 'C' 'G' 'T']
# ['--']


####################################################################################
####################################################################################


####################################################################################
# Normalization tables
####################################################################################

# Sorting order for chromosome column
chromosomePriorityList = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'XY', 'MT' ]


# Table for normalizing genotypes
genotypeTable = {
    '00': '--',
    'CA': 'AC',
    'TG': 'GT',
    'GC': 'CG',
    'TA': 'AT',
    'TC': 'CT',
    'GA': 'AG',

    'ID': 'DI',

    '-G': 'G',
    '-C': 'C',
    '-A': 'A',
    '-T': 'T'
}


# Table for normalizing genotype on X and Y chromosome (Males)
genotypeTableXYMales = {
    'GG': 'G',
    'CC': 'C',
    'AA': 'A',
    'TT': 'T',

    'DD': 'D',
    'II': 'I',

    'AC': '--',
    'GT': '--',
    'CG': '--',
    'AT': '--',
    'CT': '--',
    'AG': '--'
}


# Nocall list
noCalls = [ 'DD', 'II', 'DI', 'D', 'I', '--' ]
noCallsHyphen = [ 'DD', 'II', 'DI', 'D', 'I' ]


####################################################################################
####################################################################################


####################################################################################
# Output format tables
####################################################################################


##########################################
# 23AndMe v5

chromosomePriorityList23andMe = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' ]


##########################################


##########################################
# AncestryDNA v2

# Sorting order for ancestry chromosome column
chromosomePriorityListAncestry =  [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26' ]

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


##########################################


##########################################
# FamilyTreeDNA v3

# Sorting order for FamilyTreeDNA chromosome column
chromosomePriorityListFamilyTreeDNA = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'XY', 'MT', 'X' ]

# FamilyTreeDNA genotypes for Y, X, MT
genotypeTableFamilyTreeDNA = {
    'G': 'GG',
    'C': 'CC',
    'A': 'AA',
    'T': 'TT'
}


##########################################


##########################################
# LivingDNA v1.0.2

# Sorting order for LivingDNA chromosome column
chromosomePriorityListLivingDNA = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X' ]


##########################################


##########################################
# MyHeritage v1

# Sorting order for MyHeritage v1 chromosome column
chromosomePriorityListMyHeritage = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y' ]

# MyHeritage v1 genotypes for Y, X, MT
genotypeTableMyHeritage = {
    'G': 'GG',
    'C': 'CC',
    'A': 'AA',
    'T': 'TT',
    'D': 'DD',
    'I': 'II'
}


##########################################


####################################################################################
# FUNCTIONS
####################################################################################

##########################################
# Find all files in directory with
# the desired file ending from 'fileEndings'

def findDNAFiles( fileEndings: List ) -> List:
    # Get script directory
    scriptDir = os.path.dirname( os.path.realpath( __file__ ) )
    # Add inputFileDir to directory to get subdir
    scriptDir = scriptDir + inputFileDir

    # List all files in subdir and add files with fileEndings
    # and append to result list
    fileList = [ f for f in os.listdir( path=scriptDir ) ]
    result = []
    for f in fileList:
        if f.lower().endswith( fileEndings ):
            result.append( inputFileDir + f )


    return result

##########################################


##########################################
# Pre-screen file to determine DNA company
#

def prescreenDNAFile( inputDNAFile: str ) -> str:

    ##############################
    #  n = number of comment lines.
    #       23andMe v5        = 19
    #       AncestryDNA v2    = 18
    #       FamilyTreeDNA v3  = 0
    #       Living DNA v1.0.2 = 11
    #       MyHeritage v1     = 6
    #       MyHeritage v2     = 12


    #n = 0
    n = 1
    mystring = ' '

    # count lines in the file
    with open( inputDNAFile, 'r') as fp:
        for n, line in enumerate(fp):
            pass

    # look at the first n lines
    with open( inputDNAFile ) as myfile:
        head = [ next( myfile ) for x in range( n ) ]

    # put all lines in a string
    for x in head:
       mystring += ' ' + x


    return mystring

##########################################


##########################################
# Try to determine what DNA testing
# company the file originates from
#

#### NEEDS IMPROVMENT ####
def determineDNACompany( text: str, filename: str ) -> str:

    # List of company and patterns
    company_patterns = {
        '23andMe v5': r'_v5_full_',
        'AncestryDNA v2': r'ancestrydna',
        'LivingDNA v1.0.2': r'# living dna customer genotype data download file version: 1\.0\.2',
        'MyHeritage v2': r'##format=mhv1\.0',
        'MyHeritage v1': r'# myheritage dna raw data\.',
        'FamilyTreeDNA v3': r'rsid,chromosome,position,result'
    }

    # Convert to lowercase to make it easier
    filename = filename.lower()
    text = text.lower()

    # Search for pattern in file
    for company, pattern in company_patterns.items():
        if re.search(pattern, filename) or re.search(pattern, text):
            # Return company according to list depending on textpattern
            return company

    # search for pattern in filename
    # 23andMe v5
#    if '_v5_full_' in filename.lower():
#        return '23andMe v5'


    # If it cannot find a pattern, return unknown
    return 'unknown'


##########################################


##########################################
# Load DNA file into pandas dataframe
#

def loadDNAFile( file: str, company: str ) -> pd.DataFrame:

    # Create a dictionary with the file reading options for each company
    company_options = {
        '23andMe v5': {'dtype': str, 'sep': '\t', 'comment': '#', 'index_col': False, 'header': None, 'engine': 'python'},
        'AncestryDNA v2': {'dtype': str, 'sep': '\t', 'comment': '#'},
        'FamilyTreeDNA v3': {'dtype': str, 'comment': '#'},
        'LivingDNA v1.0.2': {'dtype': str, 'sep': '\t', 'comment': '#', 'index_col': False, 'header': None, 'engine': 'python'},
        'MyHeritage v1': {'dtype': str, 'comment': '#'},
        'MyHeritage v2': {'dtype': str, 'comment': '#'}
    }
    
    # Check if the company name is valid
    if company not in company_options:
        raise ValueError(f"Invalid company name: {company}")
    
    # Load input file into pandas using the company-specific options
    df = pd.read_csv(file, **company_options[company])


    return df


##########################################


##########################################
# Normalize the  DNA file
#

def normalizeDNAFile( df: pd.DataFrame, company: str ) -> pd.DataFrame:

    # AncestryDNA v2
    if company == 'AncestryDNA v2':
        # Merge allele1 and allele2 to genotype column
        df[ 'genotype' ] = df[ 'allele1' ] + df[ 'allele2' ]
        df = df.drop( [ 'allele1', 'allele2' ], axis=1 )

    # Normalize column names
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    # Add company column to kit
    df[ 'company' ] = company

    # Normalize chromosome order with custom chromosomeTable
    df[ 'chromosome' ] = df[ 'chromosome' ].replace( to_replace=chromosomeTableAncestryIn )

    # Normalize genotype with custom genotypeTable
    df[ 'genotype' ] = df[ 'genotype' ].replace( to_replace=genotypeTable )

    df = df.astype( {'rsid': str, 'chromosome': str, 'position': int, 'genotype': str, 'company': str} )


    return df


##########################################


##########################################
# Guesses the gender based on the genotype data in chromosome X/23.
#

def guessGenderFromDataframe( df: pd.DataFrame, company: str ) -> str:

    # Filter the DataFrame to include only chromosome X/23
    df_chr23 = df[df['chromosome'] == 'X']
    
    # Count the number of heterozygous SNPs on chromosome 23
    hetero_count = df_chr23['genotype'].str.contains('/').sum()
    
    # Guess the gender based on the proportion of homozygous SNPs on chromosome 23
    if hetero_count / len(df_chr23) < 0.05:
        gender = 'Male'
    else:
        gender = 'Female'


    return gender


##########################################


##########################################
# Clean file and normalize chromosome
# and genotype

def cleanDNAFile( df: pd.DataFrame, company: str, gender: str ) -> pd.DataFrame:

    # IF position contains genotype larger than two alleles, replace with nocall '--' (clean dirty information from LivingDNA and more?)
    df.loc[ df[ 'genotype' ].str.len() > 2, 'genotype' ] = '--'

    # If position are 0 on other chromosome than 0, assume wrong read from chip and move that row to "junk" chromosome 0
    df.loc[(df['chromosome'] != '0') & (df['position'] == 0), 'chromosome'] = '0'

    # FamilyTreeDNA v3, additional step to keep chromosome 0
    if company == 'FamilyTreeDNA v3':
        # Drop chromosome 0 rows, likely nocalls or incomplete information
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )


    return df


##########################################


##########################################
# Sort file based on custom chromosome order,
# position and custom genotype order

def sortDNAFile( df: pd.DataFrame ) -> pd.DataFrame:
    # Custom sorting order on chromosome and company column. Modify at top of file.
    df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList )
    df[ 'company' ] = pd.Categorical( df[ 'company' ], companyPriorityList )

    # Sort dataframe based on custom sorting orders and position
    df.sort_values( [ 'chromosome', 'position', 'company' ], ascending=( True, True, True ), inplace=True )


    return df


##########################################


##########################################
# Drop duplicates on genotype, keeping
# only genotype according to priority list
# in companyPriorityList

def dropDuplicatesDNAFile( df: pd.DataFrame ) -> pd.DataFrame:


##### STEP 1 - Drop NoCalls only if there are duplicate rows with atleast one genotype that is not a nocall #####
    print()
    print( 'Drop nocall if there is a non nocall genotype on duplicate position' )
    print()
    # Create a boolean mask for rows where position is not a duplicate within each chromosome
    mask = df.groupby(['chromosome', 'position']).genotype.transform('nunique') == 1

    # Select rows where the mask is True or genotype != '--'
    df = df[mask | df.genotype.ne('--')]
    print( 'DONE!' )
    print()



##### STEP 2 - Do a majority vote on the duplicates and choose the genotype that has the most of the same #####

    if majorityVote == True:
        print()
        print( 'Drop based on majority vote' )
        print()
        # Group the data frame by chromosome and position columns and apply a lambda function to each group
        df = df.groupby([df['chromosome'], 'position']).apply(
            lambda x: 
                # If the majority genotype count is greater than or equal to 50%, filter the group to keep only rows with the majority genotype and return the first row
                x[x['genotype'] == x['genotype'].mode().iloc[0]].iloc[:1] if (x['genotype'] == x['genotype'].mode().iloc[0]).sum()/x.shape[0] >= 0.5 
                # Otherwise, return the original group
                else x
            # Reset the index of the resulting data frame after grouping and filtering
        ).reset_index(drop=True)
        print( 'DONE!' )
        print()



##### STEP 3 - Drop duplicates and save only the first row. Which genotype that is first are determined by companyPriorityList #####

    # If genotype is different on the same position, then only keep the genotype from the company according to the order in companyPriorityList (which got sorted earlier)
    print()
    print( 'Keep first duplicate, drop the rest' )
    print()
    df = df.drop_duplicates(subset=['chromosome', 'position'], keep='first')
    print( 'DONE!' )
    print()


    return df


##########################################


##########################################
# Restore original RSID names to the
# DNA Superkit file, based on outputFormat

def restoreRSIDtoDNAFile( df: pd.DataFrame, outputFormat: str ) -> pd.DataFrame:

    # Load the data file
    df_replace = pd.read_csv( './data/' + outputFormat + '.df', dtype=str, sep='\t' )

    # Set the index of both dataframes to (chromosome, position)
    df.set_index( ['chromosome', 'position'], inplace=True )
    df_replace.set_index( ['chromosome', 'position'], inplace=True )

    # Update values in the original dataframe with values from the replacement dataframe
    df.update( df_replace )

    # Reset the index of the dataframe to default integer index
    df.reset_index( inplace=True )

    # Reorder columns to original format
    df = df.reindex(columns=['rsid', 'chromosome', 'position', 'genotype'])

    return df

##########################################


##########################################
# prepare database for company specific
# output format

def formatDNAFile( df: pd.DataFrame, company: str ) -> pd.DataFrame:

    # 23andMe v5
    if company == '23andMe v5':

######### DROP #########
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )

######### SORTING #########
        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList23andMe )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

######### RENAME COLUMNS #########
        # Rename columns
        df.rename( columns = { 'rsid':'# rsid' }, inplace = True )



    # AncestryDNA v2
    elif company == 'AncestryDNA v2':

######### DROP #########
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )

######### NORMALIZE #########
        # Normalize chromosome order with custom chromosomeTable
        df[ 'chromosome' ].replace( to_replace=chromosomeTableAncestryOut, inplace=True )

######### SORTING #########
        # Custom sorting order on chromosome and company column. Modify at top of file.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListAncestry )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position', ], ascending=( True, True ), inplace=True )

######### NOCALLS #########
        # Changed nocalls from -- to 00
        df[ 'genotype' ].replace( to_replace=genotypeTableAncestry, inplace=True )

######### ALLELES #########
        # Split genotype into allele1 and allele2
        df[ 'allele1' ] = df[ 'genotype' ].str[ :1 ]
        df[ 'allele2' ] = df[ 'genotype' ].str[ -1: ]
        del df[ 'genotype' ]



    # FamilyTreeDNA v3
    elif company == 'FamilyTreeDNA v3':

######### ADD CHROMOSOME 0 #########
        # Concat dataframe with previously dropped chromosome 0
        df = pd.concat( [df, chromosomeZero] , sort=False, ignore_index=True)

######### DROP #########
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == 'Y' ].index )

######### NORMALIZE #########
        #Replace genotypes for X, Y, MT( A = -A )
        df[ 'genotype' ].replace( to_replace=genotypeTableFamilyTreeDNA, inplace=True )

######### NOCALLS #########
        # Drop nocalls according to noCallsHyphen
        for f in noCallsHyphen:
            indexNames = df[ df[ 'genotype' ] == f ].index
            df.drop( indexNames, inplace = True )
#        df = df[ ~df[ 'genotype' ].isin( noCallsHyphen ) ]

######### SORTING #########
        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListFamilyTreeDNA )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

######### RENAME COLUMNS #########
        # Rename columns
        df.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )



    # Living DNA v1.0.2
    elif company == 'LivingDNA v1.0.2':

######### DROP #########
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'MT' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'Y' ].index )

######### NOCALLS #########
        # Drop nocalls according to noCalls
        for f in noCalls:
            indexNames = df[ df[ 'genotype' ] == f ].index
            df.drop( indexNames, inplace = True )
#        df = df[ ~df[ 'genotype' ].isin( noCalls ) ]

######### SORTING #########
        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListLivingDNA )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

######### RENAME COLUMNS #########
        # Rename columns
        df.rename( columns = { 'rsid':'# rsid' }, inplace = True )



    # MyHeritage v1 & MyHeritage v2
    elif company == 'MyHeritage v1' or company == 'MyHeritage v2':

######### DROP #########
        # Drop chromosomes that arent used
        df = df.drop( df[ df[ 'chromosome' ] == '0' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'XY' ].index )
        df = df.drop( df[ df[ 'chromosome' ] == 'MT' ].index )

######### NORMALIZE #########
        #Replace genotypes for X, Y, MT( A = AA )
        df[ 'genotype' ].replace( to_replace=genotypeTableMyHeritage, inplace=True )

######### NOCALLS #########
        if company == 'MyHeritage v1':
        # Drop nocalls according to noCallsHyphen
            for f in noCallsHyphen:
                indexNames = df[ df[ 'genotype' ] == f ].index
                df.drop( indexNames, inplace = True )
#            df = df[ ~df[ 'genotype' ].isin( noCallsHyphen ) ]

######### SORTING #########
        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityListMyHeritage )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )

        # Make all columns strings
        df = df.astype( {'rsid': str, 'chromosome': str, 'position': str, 'genotype': str} )

######### RENAME COLUMNS #########
        # Rename columns
        df.rename( columns = { 'rsid':'RSID', 'chromosome':'CHROMOSOME', 'position':'POSITION', 'genotype':'RESULT' }, inplace = True )



    # SuperKit format
    elif company == 'SuperKit':

######### ADD CHROMOSOME 0 #########
        # Concat dataframe with previously dropped chromosome 0
        df = pd.concat( [df, chromosomeZero] , sort=False, ignore_index=True)

######### SORTING #########
        # Custom sorting order on chromosome and company column.
        df[ 'chromosome' ] = pd.Categorical( df[ 'chromosome' ], chromosomePriorityList )
        # Sort frame based on custom sorting orders and position
        df.sort_values( [ 'chromosome', 'position' ], ascending=( True, True ), inplace=True )


    return df

##########################################


##########################################
# Prepare database for company specific
# output format

def trimSNPs( df: pd.DataFrame, company: str ) -> pd.DataFrame:

    # Define trim mask
    mask = None

    # 23andMe v5
    if company == '23andMe v5':
        # Tested ranges
        chromosome_ranges = {
            '1': (69869, 249222527),
            '2': (11944, 243041411),
            '3': (66206, 197884212),
            '4': (73071, 190937862),
            '5': (10735, 180715140),
            '6': (100114, 170904808),
            '7': (43748, 159124173),
            '8': (162363, 146300855),
            '9': (116581, 141101939),
            '10': (95844, 135473014),
            '11': (128623, 134945120),
            '12': (68276, 133838353),
            '13': (19025153, 115106996),
            '14': (19000060, 107283150),
            '15': (20004966, 102462479),
            '16': (68031, 90161959),
            '17': (1389, 81162706),
            '18': (11358, 78010620),
            '19': (110925, 59097933),
            '20': (61098, 62954925),
            '21': (9418848, 48099610),
            '22': (16055122, 51214796),

            'X': (167699, 155234707),
            'Y': (1726773, 58997068),
            'MT': (10, 16545)
        }
        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['chromosome'] == chromosome) & (df['position'] >= start) & (df['position'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask


    # AncestryDNA v2
    elif company == 'AncestryDNA v2':
        # Tested ranges
        chromosome_ranges = {
            '1': (752721, 249208772),
            '2': (15672, 243068403),
            '3': (61044, 197837967),
            '4': (71566, 190985642),
            '5': (38139, 180690937),
            '6': (149636, 170919470),
            '7': (42920, 159123700),
            '8': (164984, 146292734),
            '9': (46587, 141026318),
            '10': (93244, 135475096),
            '11': (193808, 134939292),
            '12': (188285, 133810854),
            '13': (19121950, 115091346),
            '14': (19120274, 107282973),
            '15': (20028058, 102398060),
            '16': (101263, 90170495),
            '17': (8547, 81047721),
            '18': (13034, 78007784),
            '19': (249357, 59097752),
            '20': (63231, 62960292),
            '21': (10862919, 48090629),
            '22': (16869887, 51213613),

            '23': (2703633, 154916845),
            '24': (2655180, 58883690),
            '25': (170770, 59330613),
            '26': (523, 16391)
        }
        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['chromosome'] == chromosome) & (df['position'] >= start) & (df['position'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask


    # FamilyTreeDNA v3
    elif company == 'FamilyTreeDNA v3':
        # Tested ranges
        chromosome_ranges = {
            '1': (72526, 249236761),
            '2': (11944, 243178150),
            '3': (66206, 197901799),
            '4': (37263, 191014647),
            '5': (14782, 180715140),
            '6': (63979, 170919470),
            '7': (43748, 159124173),
            '8': (33142, 146300855),
            '9': (46587, 141103966),
            '10': (95844, 135490953),
            '11': (124655, 134945120),
            '12': (68276, 133820694),
            '13': (19025153, 115103150),
            '14': (19000060, 107283150),
            '15': (20004966, 102516235),
            '16': (89659, 90274695),
            '17': (2220, 81125758),
            '18': (11358, 78010620),
            '19': (247265, 59097933),
            '20': (63231, 62960292),
            '21': (9432346, 48099610),
            '22': (16057310, 51214796),

            'XY': (153977, 155234707),
            'MT': (56, 16482),
            'X': (2699898, 154854338)
        }
        mask = (df['CHROMOSOME'] == '0')

        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['CHROMOSOME'] == chromosome) & (df['POSITION'] >= start) & (df['POSITION'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask


    # Living DNA v1.0.2
    elif company == 'LivingDNA v1.0.2':
        # Tested ranges
        chromosome_ranges = {
            '1': (752721, 249222527),
            '2': (11944, 243044147),
            '3': (63411, 197863379),
            '4': (71566, 190904950),
            '5': (14782, 180715140),
            '6': (165391, 170919470),
            '7': (43748, 159124173),
            '8': (164984, 146292681),
            '9': (126903, 141098720),
            '10': (95844, 135434303),
            '11': (200948, 134945120),
            '12': (190980, 133777645),
            '13': (19020095, 115103150),
            '14': (19264875, 107282024),
            '15': (20044342, 102397317),
            '16': (89659, 90153346),
            '17': (13905, 81096636),
            '18': (43621, 78001282),
            '19': (267039, 59097933),
            '20': (61098, 62957652),
            '21': (10971951, 48099610),
            '22': (16061342, 51211383),

            'X': (170193, 155230724)
        }
        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['chromosome'] == chromosome) & (df['position'] >= start) & (df['position'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask


    # MyHeritage v1
    elif company == 'MyHeritage v1':
        # Tested ranges
        chromosome_ranges = {
            '1': (82154, 249218992),
            '2': (18674, 243048760),
            '3': (61495, 197838262),
            '4': (71566, 190915650),
            '5': (25328, 180693127),
            '6': (203878, 170919470),
            '7': (44935, 159119486),
            '8': (164984, 146293414),
            '9': (46587, 141066491),
            '10': (98087, 135477883),
            '11': (198510, 134934063),
            '12': (191619, 133777645),
            '13': (19058717, 115103529),
            '14': (19255726, 107287663),
            '15': (20071673, 102461162),
            '16': (88165, 90163275),
            '17': (8547, 81046413),
            '18': (13034, 78015180),
            '19': (260912, 59097160),
            '20': (63244, 62912463),
            '21': (10827533, 48100155),
            '22': (16114244, 51211392),

            'X': (1410495, 154916845),
            'Y': (2655180, 59032197)
        }
        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['CHROMOSOME'] == chromosome) & (df['POSITION'] >= start) & (df['POSITION'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask


    # MyHeritage v2
    elif company == 'MyHeritage v2':
        # Tested ranges
        chromosome_ranges = {
            '1': (72526, 249222527),
            '2': (11944, 243048760),
            '3': (66206, 197884574),
            '4': (73071, 190963766),
            '5': (14782, 180715140),
            '6': (174798, 170919470),
            '7': (43748, 159124173),
            '8': (170692, 146292681),
            '9': (46587, 141068637),
            '10': (95844, 135434303),
            '11': (198986, 134945120),
            '12': (190980, 133820694),
            '13': (19025153, 115103150),
            '14': (19108923, 107282024),
            '15': (20004966, 102434534),
            '16': (89659, 90234618),
            '17': (2220, 81104884),
            '18': (11358, 78010620),
            '19': (247265, 59097933),
            '20': (63231, 62960292),
            '21': (9922018, 48099610),
            '22': (16139442, 51214796),

            'X': (2699898, 154913173),
            'Y': (2654333, 58997092)
        }
        # Trim chromosomes according to list above
        for chromosome, (start, end) in chromosome_ranges.items():
            current_mask = (df['CHROMOSOME'] == chromosome) & (df['POSITION'] >= start) & (df['POSITION'] <= end)
            if mask is None:
                mask = current_mask
            else:
                mask = mask | current_mask



    # Do actual trimming
    if company != 'SuperKit':
        df = df[ mask ]


    return df

##########################################


##########################################
# Adds comments on top of the DNA file
# that mimics the original comments of the output format

def addCommentsToFile( tmpFileName: str, outputFormat: str):

    ### ADD COMMENTS
    print()
    print("Add comments to file")
    # Add comments to file
    if outputFormat == '23andMe v5':
        # Get current time
        now = datetime.datetime.utcnow().strftime('%a %b %d %H:%M:%S %Y')

        # Generate a random integer
        rand_int = random.randint(0, 2**64-1)
        hex_str = format(rand_int, '016x')

        # Top comment
        data = (
            f"# This data file generated by 23andMe at: {now}\n"
            "#\n"
            "# This file contains raw genotype data, including data that is not used in 23andMe reports.\n"
            "# This data has undergone a general quality review however only a subset of markers have been \n"
            "# individually validated for accuracy. As such, this data is suitable only for research, \n"
            "# educational, and informational use and not for medical or other use.\n"
            "# \n"
            "# Below is a text version of your data.  Fields are TAB-separated\n"
            "# Each line corresponds to a single SNP.  For each SNP, we provide its identifier \n"
            "# (an rsid or an internal id), its location on the reference human genome, and the \n"
            "# genotype call oriented with respect to the plus strand on the human reference sequence.\n"
            "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
            "# Note that it is possible that data downloaded at different times may be different due to ongoing \n"
            "# improvements in our ability to call genotypes. More information about these changes can be found at:\n"
            f"# https://you.23andme.com/p/{hex_str}/tools/data/download/\n"
            "# \n"
            "# More information on reference human assembly builds:\n"
            "# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/\n"
            "#\n"
        )

        # Read the existing contents of the file
        with open(tmpFileName, "r") as f:
            existing_data = f.read()

        # Write the new data followed by the original contents
        with open(tmpFileName, "w", newline="\r\n") as f:
            f.write(data + existing_data)


    elif outputFormat == 'AncestryDNA v2':
        # Get current time
        now = datetime.datetime.utcnow().strftime('%m/%d/%Y %H:%M:%S UTC')

        # Top comment
        data = (
            "#AncestryDNA raw data download\n"
            f"#This file was generated by AncestryDNA at: {now}\n"
            "#Data was collected using AncestryDNA array version: V2.0\n"
            "#Data is formatted using AncestryDNA converter version: V1.0\n"
            "#Below is a text version of your DNA file from Ancestry.com DNA, LLC.  THIS\n"
            "#INFORMATION IS FOR YOUR PERSONAL USE AND IS INTENDED FOR GENEALOGICAL RESEARCH\n"
            "#ONLY.  IT IS NOT INTENDED FOR MEDICAL, DIAGNOSTIC, OR HEALTH PURPOSES.  THE EXPORTED DATA IS\n"
            "#SUBJECT TO THE ANCESTRY TERMS AND CONDITIONS, BUT PLEASE BE AWARE THAT THE\n"
            "#DOWNLOADED DATA WILL NO LONGER BE PROTECTED BY OUR SECURITY MEASURES.\n"
            "#WHEN YOU DOWNLOAD YOUR RAW DNA DATA, YOU ASSUME ALL RISK OF STORING,\n"
            "#SECURING AND PROTECTING YOUR DATA.  FOR MORE INFORMATION, SEE ANCESTRYDNA FAQS.\n"
            "#\n"
            "#Genetic data is provided below as five TAB delimited columns.  Each line\n"
            "#corresponds to a SNP.  Column one provides the SNP identifier (rsID where\n"
            "#possible).  Columns two and three contain the chromosome and basepair position\n"
            "#of the SNP using human reference build 37.1 coordinates.  Columns four and five\n"
            "#contain the two alleles observed at this SNP (genotype).  The genotype is reported\n"
            "#on the forward (+) strand with respect to the human reference.\n"
        )

        # Read the existing contents of the file
        with open(tmpFileName, "r") as f:
            existing_data = f.read()

        # Write the new data followed by the original contents
        with open(tmpFileName, "w", newline="\r\n") as f:
            f.write(data + existing_data)


    elif outputFormat == 'LivingDNA v1.0.2':
        # Get current time
        now = datetime.datetime.utcnow()

        # Format date as month-day-year
        now = datetime.datetime.utcnow()
        if os.name == "nt":  # Windows
            month = now.strftime("%#m")
            day = now.strftime("%#d")
        else: # Unix/Linux/Mac
            month = now.strftime("%-m")
            day = now.strftime("%-d")
        year = now.strftime("%Y")
        formatted_date = f"{month}-{day}-{year}"


        # Top comment
        data = (
            "# Living DNA customer genotype data download file version: 1.0.2\n"
            f"# File creation date: {formatted_date}\n"
            "# Genotype chip: Sirius\n"
            "# The content of this file is subject to updates and changes depending on the time of download.\n"
            "# This genotype data should be treated as personal information.\n"
            "# This genotype data is not suitable for clinical/medical research or diagnosis.\n"
            "# The user assumes all responsibility for the security of this file.\n"
            "# Please refer to the Living DNA Terms and Conditions on our website (www.livingdna.com) for more information.\n"
            "# Human Genome Reference Build 37 (GRCh37.p13).\n"
            "# Genotypes are presented on the forward strand.\n"
            "#\n"
        )

        # Read the existing contents of the file
        with open(tmpFileName, "r") as f:
            existing_data = f.read()

        # Write the new data followed by the original contents
        with open(tmpFileName, "w", newline="\n") as f:
            f.write(data + existing_data)


    elif outputFormat == 'MyHeritage v1':
        # Top comment
        data = (
            "RSID,CHROMOSOME,POSITION,RESULT\n"
        )

        # Read the existing contents of the file
        with open(tmpFileName, "r") as f:
            existing_data = f.read()

        # Write the new data followed by the original contents
        with open(tmpFileName, "w", newline="\r\n") as f:
            f.write(data + existing_data)


    elif outputFormat == 'MyHeritage v2':
        # Get current time
#        now = datetime.datetime.utcnow().strftime('%m/%d/%Y %H:%M:%S UTC')
        now = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')

        # Top comment
        data = (
            "##fileformat=MyHeritage\n"
            "##format=MHv1.0\n"
            "##chip=GSA\n"
            f"##timestamp={now}\n"
            "##reference=build37\n"
            "#\n"
            "# MyHeritage DNA raw data.\n"
            "# For each SNP, we provide the identifier, chromosome number, base pair position and genotype.\n"
            "# The genotype is reported on the forward (+) strand with respect to the human reference build 37.\n"
            "# THIS INFORMATION IS FOR YOUR PERSONAL USE AND IS INTENDED FOR GENEALOGICAL RESEARCH\n"
            "# ONLY. IT IS NOT INTENDED FOR MEDICAL OR HEALTH PURPOSES. PLEASE BE AWARE THAT THE\n"
            "# DOWNLOADED DATA WILL NO LONGER BE PROTECTED BY OUR SECURITY MEASURES.\n"
            "RSID,CHROMOSOME,POSITION,RESULT\n"
        )

        # Read the existing contents of the file
        with open(tmpFileName, "r") as f:
            existing_data = f.read()

        # Write the new data followed by the original contents
        with open(tmpFileName, "w", newline="\r\n") as f:
            f.write(data + existing_data)


    return


####################################################################################
####################################################################################



####################################################################################
# MAIN LOOP
####################################################################################

# Get the start time
start_time = time.time()

# Find files in dir with the correct file endings
rawDNAFiles = findDNAFiles( fileEndings )

# Check if there are any files in the directory
if not rawDNAFiles:
    print()
    print ("There is no files in the directory")
    exit()


########################
# Preparing and cleaning DNA files

# empty array to put results in
resultFiles = []
chromosomeZero = pd.DataFrame()
DNACount = 0


##########################################
# Look for files and process them


for file in rawDNAFiles:

#    print( type(file) )

    # Screening file to determine company
    fileScreening = prescreenDNAFile( file )

    # Get DNA company from file comment
    company = determineDNACompany( fileScreening , file)

    if company != 'unknown':

        DNACount = DNACount + 1

        # Load the DNA file into pandas and get columns
        df = loadDNAFile( file, company )
        # Normalize the DNA file
        df = normalizeDNAFile( df, company )
        # Guess gender in kit
        guessGender = guessGenderFromDataframe( df, company )



        ################################################################
        ##### SAVE RSID, CHROMOSOME AND POSITION to ./data/ folder #####
#        df_rsid = df.copy()
#        df_rsid = df_rsid.drop_duplicates(subset=['chromosome', 'position'], keep='first')
#        df_rsid.drop(columns=['company', 'genotype'], inplace=True)
#        df_rsid.to_csv('./data/' + company + '.df', index=None, sep='\t', encoding='ascii', lineterminator='\r\n')
#        df_rsid.to_csv('./data/' + company + '-full' + '.df', index=None, sep='\t', encoding='ascii', lineterminator='\r\n')

# Just spare code
#       df_rsid['genotype'] = '--'

        ################################################################
        ################################################################



        # Normalize genotypes on X and Y (MT?) chromosomes where heterozygous calls are defined as nocalls '--'
        # (as males only have one X and one Y), and the rest are changed to a single letter
        #
        # Add commandline to bypass this check to handle mutations?
        # HANDLE DI calls?
        if guessGender == 'Male':
            rep = df[ 'genotype' ].replace( genotypeTableXYMales )
            df[ 'genotype' ] = df[ 'genotype' ].mask( df[ 'chromosome' ].isin( [ 'X','Y', 'MT' ] ), rep )

        # Workaround to keep Chromosome 0 (nocalls? bad data?)
        if company == 'FamilyTreeDNA v3':
            chromosomeZero = df.loc[ df[ 'chromosome' ] == '0' ]
            del chromosomeZero[ 'company' ]

        # Clean dataframe
        df = cleanDNAFile( df, company, guessGender )

        # Concatenate DNA data
        resultFiles.append( df )

        # Presenting results
        print()
        print( '######################################################################')
        print( "#")
        print( f"# Testcompany:            {company}" )
        print( "#" )
        print( f"# File:                   {file.replace( inputFileDir, '' )}")
        print( f"# SNPs tested in kit:     {len(df)}")
        print( f"# Assumed gender in kit:  {guessGender}" )
        print( "#")
        print( '######################################################################')
        print()
        print( f"Chromosomes: {df.chromosome.unique().tolist()}" )
        print()

    # If file is unknown
    else:
        print()
        print( '######################################################################')
        print( "#")
        print( f"# Testcompany:            {company}" )
        print( "#" )
        print( f"# File:                   {file.replace( inputFileDir, '' )}")
        print( "#")
        print( '######################################################################')
        print()


# Check if there are objects in DNASuperKit
# if not, then quit script
if not resultFiles:
    print ("No compatible files has been found")

    exit()


##########################################
##########################################


########################
# Concatenate and remove duplicates

print()
print( '######################################################################' )
print( "#" )
print( "# Concatenating files, sorting list and dropping duplicates" )
print( "#" )
print( '######################################################################' )


print()
print( f"# Concatenating {DNACount} DNA files" )
print()

# Concatenate all DNA files into one list
DNASuperKit = pd.concat(resultFiles, sort=False, ignore_index=True)
print( "DONE!" )
print()


print()
print( "Sorting resulting dataframe" )
print()

# Sort DNA according to order provided in customization
DNASuperKit = sortDNAFile( DNASuperKit )
print( "DONE!" )
print()


print()
print( 'Dropping duplicates' )
print()

# Drop duplicates
DNASuperKit = dropDuplicatesDNAFile( DNASuperKit )
print( "DONE!" )
print()


# Count SNPs per included per company
companySNPCounts = []

for f in companyPriorityList:
    companySNPCount = DNASuperKit['company'].value_counts()[f].astype(str)
    companySNPCounts.append(f + ': ' + companySNPCount)


# Delete 'company' column
del DNASuperKit[ 'company' ]

########################


########################
#Format dataframe to a specific company format

print()
print( '######################################################################' )
print( "#" )
print( "# Formatting, trimming and saving data" )
print( "#" )
print( '######################################################################' )


if outputFormat != 'SuperKit':
    print()
    print( f'Restoring {outputFormat} RSID to DNA file to' )
    print()

    # Restore original RSID according to outputFormat
    DNASuperKit = restoreRSIDtoDNAFile( DNASuperKit, outputFormat)

    print( "DONE!" )
    print()

# Nr of SNPs in kit
superkitLength = len(DNASuperKit)

print()
print( f'Formatting DNA file to {outputFormat} format' )
print()

# Format DNA file to match desired output structure
DNASuperKit = formatDNAFile( DNASuperKit, outputFormat )
print( "DONE!" )
print()


if trimSNP == True:
    print()
    print( f'Trimming SNPs to match {outputFormat} SNP ranges' )
    print()

    # Trim SNPs
    DNASuperKit = trimSNPs( DNASuperKit, outputFormat)
    print( "DONE!" )
    print()


# Line terminators: 
# \n = LF (Linux), \r\n = CRLF (Windows)
formats = {
    '23andMe v5': { 'sep': '\t', 'encoding': 'ascii', 'lineterminator': '\r\n' },
    'AncestryDNA v2': { 'sep': '\t', 'encoding': 'ascii', 'lineterminator': '\r\n' },
    'FamilyTreeDNA v3': { 'sep': ',', 'encoding': 'ascii', 'lineterminator': '\n' },
    'LivingDNA v1.0.2': { 'sep': '\t', 'encoding': 'ascii', 'lineterminator': '\n' },
    'MyHeritage v1': { 'header': False, 'sep': ',', 'encoding': 'ascii', 'lineterminator': '\n', 'quoting': 2 },
    'MyHeritage v2': { 'header': False, 'sep': ',', 'encoding': 'ascii', 'lineterminator': '\n', 'quoting': 2 },
    'SuperKit': { 'sep': '\t', 'encoding': 'ascii', 'lineterminator': '\r\n' }
}

# Handle unsupported format (shouldn't be possible though)
if outputFormat not in formats:
    raise ValueError(f"Unsupported format: {outputFormat}")

# Set correct file ending
if outputFormat in ['AncestryDNA v2', 'LivingDNA v1.0.2', '23andMe v5', 'SuperKit']:
    ext = 'txt'
else:
    ext = 'csv'

# File directory + filename to one string variable
tmpFileName = f"{outputFileDir}{outputFileName}-{outputFormat}.{ext}"


# Save to file
print( f'Saving DNA Superkit to {outputFormat} format.' )
DNASuperKit.to_csv(tmpFileName, index=None, **formats[outputFormat])
print( "DONE!" )
print()


# Add comments to file
#if outputFormat == '23andMe v5' or outputFormat == 'AncestryDNA v2' or outputFormat == 'LivingDNA v1.0.2' or outputFormat == 'MyHeritage v1' or outputFormat == 'MyHeritage v2':
companiesWithCommentsToAdd = ['23andMe v5',
                              'AncestryDNA v2',
                              'LivingDNA v1.0.2',
                              'MyHeritage v1',
                              'MyHeritage v2']

if outputFormat in companiesWithCommentsToAdd:

    print()
    print( f'Adding comments to top of file according to {outputFormat} format' )
    print()

    # Adding comment
    addCommentsToFile( tmpFileName, outputFormat )
    print( "DONE!" )
    print()


# Get the end time
end_time = time.time()
# Calculate the elapsed time
elapsed_time = end_time - start_time


print()
print( 'DNA SuperKit successfully created!' )
print( f'time elapsed since start of script: {elapsed_time} seconds')
print()


print()
# Presenting results
print()
print( '######################################################################')
print( '#')
print( '# DNA SuperKit Statistics')
print( '#')
print( f'# Outputformat:           {outputFormat}' )
print( '#' )
print( f'# File:                   {tmpFileName}')
print( f'# Total nr of SNPs:       {superkitLength}')
print( '#')
print( '######################################################################')
print()
print( f'Total SNP per company: {companySNPCounts}')


if outputFormat == 'AncestryDNA v2':
    superkitUniqueChromosomes = DNASuperKit.chromosome.unique().tolist()
    # Merge allele1 and allele2 to genotype column
    DNASuperKit[ 'genotype' ] = DNASuperKit[ 'allele1' ] + DNASuperKit[ 'allele2' ]
    DNASuperKit = DNASuperKit.drop( [ 'allele1', 'allele2' ], axis=1 )
    superkitUniqueGenotypes = DNASuperKit.genotype.unique().tolist()

elif outputFormat == 'FamilyTreeDNA v3' or outputFormat == 'MyHeritage v1' or outputFormat == 'MyHeritage v2':
    superkitUniqueChromosomes = DNASuperKit.CHROMOSOME.unique().tolist()
    superkitUniqueGenotypes = DNASuperKit.RESULT.unique().tolist()

else:
    superkitUniqueChromosomes = DNASuperKit.chromosome.unique().tolist()
    superkitUniqueGenotypes = DNASuperKit.genotype.unique().tolist()

# Information about the results
print()
print( f'Chromosome List: {superkitUniqueChromosomes}')
print()
print( f'Genotype List: {superkitUniqueGenotypes}')
print()



####################################################################################
# EOF #
####################################################################################