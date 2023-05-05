##############################################################################################
# Analyze DNA files
#

import os                   # For findDNAFiles
from typing import List
import pandas as pd
import re                   # For determineDNACompany

import argparse             # Command line argument parser

####################################################################################
# COMMAND LINE ARGUMENT PARSER
####################################################################################

saveStructure = False
saveDuplicates = False

# Parser arguments
parser = argparse.ArgumentParser( formatter_class=argparse.RawTextHelpFormatter )
parser.add_argument('-ss', '--saveStructure', action='store_true', help='Save DNA file structure (without genotype) to a .df file in the ./data/ directory.', required=False)
parser.add_argument('-sd', '--saveDuplicates', action='store_true', help='Save DNA file duplicate rows to a .df file in the ./data/ directory.', required=False)

# Get arguments from command line
args = parser.parse_args()

# Save the arguments to variables
saveStructure = args.saveStructure
saveDuplicates = args.saveDuplicates


####################################################################################
####################################################################################



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

# Genotype list
genotypeList = [
                'AA', 'CC', 'GG', 'TT',

                'AC', 'AG', 'AT',
                'CA', 'GA', 'TA',

                'CG', 'CT',
                'GT', 'TC',

                'GC',
                'TG',

                'A', 'C', 'G', 'T',
                '-G', '-C', '-A', '-T',

                'DD', 'II',
                'DI', 'ID',
                'D', 'I',

                '--',
                '00'
                ]

####################################################################################
####################################################################################



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
    #       23andMe v5 = 19
    #       AncestryDNA v2 = 18
    #       FamilyTreeDNA v3 = 0
    #       Living DNA v1.0.2 = 11
    #       MyHeritage v1 = 6
    #       MyHeritage v2 = 12



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

def determineDNACompany( text: str, filename: str ) -> str:

    # List of company and patterns
    company_patterns = {
        '23andMe v5': r'_v5_full_',
#        'AncestryDNA v2': r'ancestrydna',
        'AncestryDNA v2': r'ancestrydna array version: v2\.0',
        'LivingDNA v1.0.2': r'# living dna customer genotype data download file version: 1\.0\.2',
        'MyHeritage v2': r'##format=mhv1\.0',
        'MyHeritage v1': r'# myheritage dna raw data\.',
        'FamilyTreeDNA v3': r'rsid,chromosome,position,result',
        'tellmeGen v4': r'# rsid	chromosome	position	genotype'
    }

    # Convert to lowercase to make it easier
    filename = filename.lower()
    text = text.lower()

    # Search for pattern in file
    for company, pattern in company_patterns.items():
        if re.search(pattern, filename) or re.search(pattern, text):
            return company


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
        'MyHeritage v2': {'dtype': str, 'comment': '#'},
        'tellmeGen v4': {'dtype': str, 'sep': '\t', 'comment': '#', 'index_col': False, 'header': None, 'engine': 'python'},
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
        df.drop(['allele1', 'allele2'], axis=1, inplace=True)

    # Normalize column names
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    df['chromosome'] = df['chromosome'].astype(str)
    df['position'] = df['position'].astype(int)


    return df


##########################################


##########################################
# Guesses the gender based on the genotype data in chromosome X/23.
#

def guessGenderFromDataframe( df: pd.DataFrame, company: str ) -> str:

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


##########################################
# Function to save the structure of the DNA file to a .df template file
#
def saveDNAStructureToFile( df: pd.DataFrame, company: str ):

    # Make a copy of dataframe
    df_rsid = df.copy()

    df_rsid.drop(columns=['genotype'], inplace=True)
    df_rsid.to_csv('./data/' + company + '.df', index=None, sep='\t', encoding='ascii', lineterminator='\r\n')


    return

####################################################################################
####################################################################################


##########################################
# Function to save the duplicate rows in the DNA file to a .df template file
#

def saveDNAFileDuplicates( df: pd.DataFrame, company: str ):

    # Make a copy of dataframe
    df_duplicates = df.copy()

    # group the rows based on 'chromosome' and 'position'
    groups = df_duplicates.groupby(['chromosome', 'position'])

    # filter the groups that have more than one row (i.e., duplicates)
    duplicates = groups.filter(lambda x: len(x) > 1)

    # concatenate the filtered groups into a new dataframe
    duplicates_df = pd.concat([duplicates])

    # Save to file
    duplicates_df.to_csv('./data/' + company + '-duplicates' + '.df', index=None, sep='\t', encoding='ascii', lineterminator='\r\n')


    return

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


        # Save DNA file structure as a .df template to ./data/ folder
        if saveStructure == True:
            saveDNAStructureToFile( df, company )

        # Save DNA file duplicate rows to a .df file in ./data/ folder
        if saveDuplicates == True:
            saveDNAFileDuplicates( df, company )


        # Get a count of nocalls
        # These companies have 00 as nocalls
        if company == 'AncestryDNA v2':
            Nocall_count = df['genotype'].value_counts()['00']
        # There is no noccalls in these companys kits
        elif company == 'LivingDNA v1.0.2':
            Nocall_count = 0
        # The rest of the companies use -- as nocalls
        else:
            Nocall_count = df['genotype'].value_counts()['--']
        # Calculate percentage of nocalls
        Nocalls_Percentage = round( Nocall_count / len(df) * 100, 2 )

        # Get nr of duplicate positions
        # group the rows based on 'chromosome' and 'position'
        duplicate_count_group = df[df['chromosome'] != '0'].groupby(['chromosome', 'position'])
        # filter the groups that have more than one row (i.e., duplicates)
        duplicates_count_filter = duplicate_count_group.filter(lambda x: len(x) > 1)
        # concatenate the filtered groups into a new dataframe, count rows and divide by two (to get the nr of real duplicates)
        duplicates_count = int( len( pd.concat([duplicates_count_filter]) ) / 2 )
        # Calculate percentage of duplicates
        duplicates_percentage = round( duplicates_count / len( df[df['chromosome'] != '0'] ) * 100, 2 )
  

        # This will return a series with the count of each unique value in the genotype column
        genotype_counts_all = df["genotype"].value_counts()
        # Now we can filter the counts for the genotypes in the genotypeList
        filtered_counts_all = genotype_counts_all[genotype_counts_all.index.isin(genotypeList)]


        # Presenting results
        print()
        print( '#' * 70)
        print( f'#')
        print( f'# Testcompany:                {company}' )
        print( f'#' )
        print( f'# File:                       {file.replace( inputFileDir, "" )}' )
        print( f'# Assumed gender in kit:      {guessGender}' )
        print( '#')
        print( f'# SNPs tested in kit:         {len(df)}' )
        print( f'# Number of nocalls in kit:   {Nocall_count} / {Nocalls_Percentage}%' )
        print( f'# Duplicate positions in kit: {duplicates_count} / {duplicates_percentage}%' )
        print( f'#')
        print( '#' * 70)
        print()
        print( f'Chromosomes: {df.chromosome.unique().tolist()}' )
        print( f'Genotypes: {df.genotype.unique().tolist()}' )
        print( 'Occurances of genotypes in the file:')
        for genotype, count in filtered_counts_all.items():
            print(f"{genotype}: {count}")
        print()


        # Chromosome 0 data
        filtered_df = df[df['chromosome'] == '0']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome 0")
            print( f"Unique genotypes: {unique_genotypes}")
            print()
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")


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
        print()
        chromosome_count = 1
        for chromosome in range(1, 23):
            filtered_df = df[df['chromosome'] == str(chromosome)]
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"Chromosome {chromosome_count}")
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")
            chromosome_count = chromosome_count + 1
        

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
            print()
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")
        

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
            print()
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")
        

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
            print()
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")
        

        # Chromosome MT (26) data
        if company == 'AncestryDNA v2':
            filtered_df = df[df['chromosome'] == '26']
        else:
            filtered_df = df[df['chromosome'] == 'MT']
        unique_genotypes = filtered_df['genotype'].unique()
        if len(unique_genotypes) > 0:
            print()
            print( "Chromosome MT (26)")
            print( f"Unique genotypes: {unique_genotypes}" )
            print()
            snp_min = filtered_df['position'].min()
            snp_max = filtered_df['position'].max()
            snp_count = filtered_df['position'].count()
            print( f"SNP range is between {snp_min} and {snp_max}" )
            print( f"Total tested SNPs: {snp_count}")


        # Analysing nr of occurances of each genotype on chromosomes x, Y and MT (to measure errors)
        if company =='AncestryDNA v2':
            # filter rows with chromosomes X, Y, and MT
            df_unique_genotype_filter = df[df['chromosome'].isin(['23', '24', '26'])]
        else:
            # filter rows with chromosomes X, Y, and MT
            df_unique_genotype_filter = df[df['chromosome'].isin(['X', 'Y', 'MT'])]

        # get unique genotypes and count their occurrences
        genotype_counts_XYMT = df_unique_genotype_filter.groupby('genotype').size()


        print()
        print( f'Unique genotypes and their occurance on X, Y and MT')
        print(genotype_counts_XYMT)
        print()


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
