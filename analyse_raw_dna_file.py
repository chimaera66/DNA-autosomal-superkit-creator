##############################################################################################
# Analyze DNA files
#

import os                   # For findDNAFiles
from typing import List
import pandas as pd
import re                   # For determineDNACompany


####################################################################################
# VARIABLES
####################################################################################


inputFileDir = './input/'
outputFileDir = './output/'

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)


####################################################################################
####################################################################################



####################################################################################
# FUNCTIONS
####################################################################################

##########################################
# Find all files in directory with
# the desired file ending from 'fileEndings'

def findDNAFiles( fileEndings ) -> List:
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

def prescreenDNAFile( inputDNAFile ):

    ##############################
    #  n = number of comment lines.
    #       23andMe v5 = 19
    #       Living DNA v1.0.2 = 11
    #       MyHeritage v1 = 6
    #       FamilyTreeDNA v3 = 0
    #       AncestryDNA v2 = 18

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

def determineDNACompany(text, filename):

    company_patterns = {
        '23andMe v5': r'_v5_Full_',
        'LivingDNA v1.0.2': r'# living dna customer genotype data download file version: 1\.0\.2',
        'MyHeritage v1': r'# myheritage dna raw data\.',
        'FamilyTreeDNA v3': r'rsid,chromosome,position,result',
        'AncestryDNA v2': r'ancestrydna'
    }

    filename = filename.lower()
    text = text.lower()

    for company, pattern in company_patterns.items():
        if re.search(pattern, filename) or re.search(pattern, text):
            return company

    return 'unknown'


##########################################


##########################################
# Load DNA file into pandas dataframe
#

def loadDNAFile( file, company: str ) -> pd.DataFrame:

    # Create a dictionary with the file reading options for each company
    company_options = {
        '23andMe v5': {'dtype': str, 'sep': '\t', 'comment': '#', 'index_col': False, 'header': None, 'engine': 'python'},
        'LivingDNA v1.0.2': {'dtype': str, 'sep': '\t', 'comment': '#', 'index_col': False, 'header': None, 'engine': 'python'},
        'MyHeritage v1': {'dtype': str, 'comment': '#'},
        'FamilyTreeDNA v3': {'dtype': str, 'comment': '#'},
        'AncestryDNA v2': {'dtype': str, 'sep': '\t', 'comment': '#'}
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

def normalizeDNAFile(df: pd.DataFrame, company: str) -> pd.DataFrame:

    # AncestryDNA v2
    if company == 'AncestryDNA v2':
        # Merge allele1 and allele2 to genotype column
        df[ 'genotype' ] = df[ 'allele1' ] + df[ 'allele2' ]
        del df[ 'allele1' ]
        del df[ 'allele2' ]

    # Normalize column names
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    return df


##########################################


####################################################################################
####################################################################################



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
# Prepare DNA files

# empty array to put results in
resultFiles = []
chromosomeZero = pd.DataFrame()

for file in rawDNAFiles:
    # Screening file to determine company
    fileScreening = prescreenDNAFile( file )

    # Get DNA company from file comment
    company = determineDNACompany( fileScreening , file)

    if company != 'unknown':
        print( 'Analyzing file:')
        print( file.replace( inputFileDir, '' ) )
        print( 'Which contains data from: ' + company )
        
        # Load the DNA file into pandas and get columns
        df = loadDNAFile( file, company )

        print()
        print( 'Columns:')
        print( df.columns.tolist() )

        # Normalize the DNA file
        df = normalizeDNAFile( df, company )

        print()
        print( 'Chromosomes:' )
        print( df.chromosome.unique() )
        print()
        print( 'Unique genotypes (on all chromosomes):')
        print( df.genotype.unique() )
        print()
        print( 'Unique genotypes on chromosome 0' )
        filtered_df = df[df['chromosome'] == '0']
        unique_genotypes = filtered_df['genotype'].unique()
        print( unique_genotypes )
        print()
        print( 'Unique genotypes on chromosomes 1 - 22' )
        excluded_chromosomes = ['0', 'X', 'Y', 'MT', 'XY']
        filtered_df = df[~df['chromosome'].isin(excluded_chromosomes)]
        unique_genotypes = filtered_df['genotype'].unique()
        print ( unique_genotypes )
        print()
        print( 'Unique genotypes on chromosome X' )
        filtered_df = df[df['chromosome'] == 'X']
        unique_genotypes = filtered_df['genotype'].unique()
        print( unique_genotypes )
        print()
        print( 'Unique genotypes on chromosome Y' )
        filtered_df = df[df['chromosome'] == 'Y']
        unique_genotypes = filtered_df['genotype'].unique()
        print( unique_genotypes )
        print()
        print( 'Unique genotypes on chromosome MT' )
        filtered_df = df[df['chromosome'] == 'MT']
        unique_genotypes = filtered_df['genotype'].unique()
        print( unique_genotypes )
        print()
        print( 'Unique genotypes on chromosome XY' )
        filtered_df = df[df['chromosome'] == 'XY']
        unique_genotypes = filtered_df['genotype'].unique()
        print( unique_genotypes )
        print()
        # Let user know processing is completed successfully
        print( 'Done analyzing the file: ' + file.replace( inputFileDir, '' ) )
        print()

    # If file is unknown
    else:
        print()
        print( 'File ' + file.replace( inputFileDir, '' ) + ' is unknown' )
        print()


####################################################################################
# EOF #
####################################################################################