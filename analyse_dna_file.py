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

def determineDNACompany(text: str, filename: str) -> str:

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

    if '_v5_full_' in filename.lower():
        return '23andMe v5'

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
        df.drop(['allele1', 'allele2'], axis=1, inplace=True)
#        del df[ 'allele1' ]
#        del df[ 'allele2' ]

    # Normalize column names
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    df['chromosome'] = df['chromosome'].astype(str)
    df['position'] = df['position'].astype(int)

    return df


##########################################


##########################################
# Guesses the gender based on the genotype data in chromosome X/23.
#

def guessGenderFromDataframe(df: pd.DataFrame, company: str) -> str:

    # Filter the DataFrame to include only chromosome X/23
    if company == 'AncestryDNA v2':
        df_chr23 = df[df['chromosome'] == '23']
    else:
        df_chr23 = df[df['chromosome'] == 'X']
    
    # Count the number of heterozygous SNPs on chromosome 23
    hetero_count = df_chr23['genotype'].str.contains('/').sum()
    
    # Guess the gender based on the proportion of homozygous SNPs on chromosome 23
    if hetero_count / len(df_chr23) < 0.05:
        gender = 'Male'
    else:
        gender = 'Female'
        
    return gender


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
        # Load the DNA file into pandas and get columns
        df = loadDNAFile( file, company )
        # Normalize the DNA file
        df = normalizeDNAFile( df, company )
        # Guess gender in kit
        guessGender = guessGenderFromDataframe( df, company )

        # Presenting results
        print()
        print( '######################################################################')
        print( f"#")
        print( f"# Testcompany:            {company}" )
        print( f"#" )
        print( f"# File:                   {file.replace( inputFileDir, '' )}")
        print( f"# SNPs tested in kit:     {len(df)}")
        print( f"# Assumed gender in kit:  {guessGender}" )
#        print( f"# Line terminator {print(df._engine.data.dialect.lineterminator)}")
        print( f"#")
        print( '######################################################################')
        print()
        print( f"Chromosomes: {df.chromosome.unique().tolist()}" )

        # Chromosome 0 data
        filtered_df = df[df['chromosome'] == '0']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome 0")
            print( f"Unique genotypes: {unique_genotypes}")
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range is between {snp_min} and {snp_max}" )

        # Chromosome 1-22 data
        if company == 'AncestryDNA v2':
            excluded_chromosomes = ['0', '23', '24', '25', '26']
        else:
            excluded_chromosomes = ['0', 'X', 'Y', 'XY', 'MT']
        filtered_df = df[~df['chromosome'].isin(excluded_chromosomes)]
        unique_genotypes = filtered_df['genotype'].unique()
        print()
        print( "Chromosomes 1 - 22" )
        print( f"Unique genotypes: {unique_genotypes}" )
        for chromosome in range(1, 23):
            filtered_df = df[df['chromosome'] == str(chromosome)]
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range for chromosome {chromosome} is between {snp_min} and {snp_max}" )
        
        # Chromosome X (23) data
        if company == 'AncestryDNA v2':
            filtered_df = df[df['chromosome'] == '23']
        else:
            filtered_df = df[df['chromosome'] == 'X']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome X (23)")
            print( f"Unique genotypes: {unique_genotypes}" )
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range is between {snp_min} and {snp_max}" )
        
        # Chromosome Y (24) data
        if company == 'AncestryDNA v2':
            filtered_df = df[df['chromosome'] == '24']
        else:
            filtered_df = df[df['chromosome'] == 'Y']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome Y (24)")
            print( f"Unique genotypes: {unique_genotypes}" )
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range is between {snp_min} and {snp_max}" )
        
        # Chromosome XY (25) data
        if company == 'AncestryDNA v2':
            filtered_df = df[df['chromosome'] == '25']
        else:
            filtered_df = df[df['chromosome'] == 'XY']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome XY (25)")
            print( f"Unique genotypes: {unique_genotypes}" )
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range is between {snp_min} and {snp_max}" )
        
        # Chromosome MT (26) data
        if company == 'AncestryDNA v2':
            filtered_df = df[df['chromosome'] == '26']
        else:
            filtered_df = df[df['chromosome'] == 'MT']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome MT (25)")
            print( f"Unique genotypes: {unique_genotypes}" )
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            print( f"SNP range is between {snp_min} and {snp_max}" )
        
        print()
        # Let user know processing is completed successfully
        print( 'Done analyzing the file: ' + file.replace( inputFileDir, '' ) )
        print()
        print()
        print()

    # If file is unknown
    else:
        print()
        print( 'File ' + file.replace( inputFileDir, '' ) + ' is unknown' )
        print()


####################################################################################
# EOF #
####################################################################################
