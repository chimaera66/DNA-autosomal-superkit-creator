##############################################################################################
# Analyze DNA files
#

import pandas as pd
import chardet               # For detecting file encoding

import os                   # For findDNAFiles
from typing import List
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

#    print( f"{result}" )

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

    # Add company column to kit
    df[ 'company' ] = company

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


##########################################
# Function to check which file termination the file has
#

def getLineTerminator( filePath: str ) -> str:

    # Open file
    with open(filePath, 'rb') as file:
        data = file.read()
    
    # Windows line terminator
    if b'\r\n' in data:
        return 'CRLF \\r\\n (Windows)'
    
    # Unix/Linux line terminator
    elif b'\n' in data:
        return 'LF \\n (Unix/Linux) '
    
    # Mac line terminator
    elif b'\r' in data:
        return 'CR \\r (Mac)'
    
    # Unknown line terminator
    else:
        return None


####################################################################################
####################################################################################


##########################################
# Function to check what file encoding the file has
#

def getFileEncoding( filePath: str ) -> str:
    with open( filePath, 'rb' ) as file:
        rawData = file.read()

        # use chardet to detect encoding
        result = chardet.detect(rawData)

        # saving encoding
        encoding = result['encoding']

        # saving confidence (not really needed)
        confidence = result['confidence']
        
        return encoding, confidence

####################################################################################
####################################################################################


##########################################
# Function to compare two dataframes and count nr of common rows
#

def compareDataframes(df1, df2):

    # Replace ancestrys chromosome numbering with standard to enable dropping
    # 23 is X, 24 is Y, 25 is the PAR/XY region, and 26 is mtDNA.
    chromosome_mapping = { '23': 'X', '24': 'Y', '25': 'XY', '26': 'MT' }
    df1['chromosome'] = df1['chromosome'].replace(chromosome_mapping)
    df2['chromosome'] = df2['chromosome'].replace(chromosome_mapping)

    # Merge dataframes on chromosome and position
    mergedDf = df1.merge(df2, on=['chromosome', 'position'], how='inner', suffixes=['_df1', '_df2'])
    mergedDf = mergedDf.drop_duplicates(subset=['chromosome', 'position'])

    # Drop chromosomes not used in autosomal dna testing
    chromosomes_to_drop = ['XY', 'MT', 'Y', 'X']
    mergedDfAutosomal = mergedDf[~mergedDf['chromosome'].isin(chromosomes_to_drop)]

    # Count the number of rows where chromosome and position match
    matchCountAutosomal = len( mergedDfAutosomal )


    # Drop all chromosomes except X
    chromosomes_to_drop = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'XY', 'MT', 'Y']
    mergedDfX = mergedDf[~mergedDf['chromosome'].isin(chromosomes_to_drop)]
    
    # Count the number of rows where chromosome and position match
    matchCountX = len( mergedDfX )

    
    return matchCountAutosomal, matchCountX

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
dataframeList = []
companyList = []
chromosomeZero = pd.DataFrame()

for file in rawDNAFiles:

    # Display current file
    print( f'Analysing file: {file.replace( inputFileDir, "" )}' )

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

        # Add dataframes to list
        dataframeList.append( df )
        # Add companies to list
        companyList.append( company )


    ########################
    # Save file structure and duplicates to file
    # for debugging reasons. Use command line arguments

        # Save DNA file structure as a .df template to ./data/ folder
        if saveStructure == True:
            saveDNAStructureToFile( df, company )

        # Save DNA file duplicate rows to a .df file in ./data/ folder
        if saveDuplicates == True:
            saveDNAFileDuplicates( df, company )


    ########################
    # Get data and statistics from DNA files

        # Check File Encoding
        fileEncoding, fileEncodingConfidence = getFileEncoding( file )

        # Check Line Terminator
        lineTerminator = getLineTerminator( file )

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


    ########################
    # Print DNA file statistics

        # number of #
        fenceNr = 70

        # Presenting results
        print()
        print( '#' * fenceNr)
        print( f'#')
        print( f'# Testcompany:                {company}' )
        print( f'#' )
        print( f'# Filename:                   {file.replace( inputFileDir, "" )}' )
        print( f'# File encoding:              {fileEncoding}, ({fileEncodingConfidence})' )
        print( f'# Line terminator:            {lineTerminator}' )
        print( f'#' )
        print( f'# Assumed gender in kit:      {guessGender}' )
        print( '#')
        print( f'# SNPs tested in kit:         {len(df)}' )
        print( f'# Number of nocalls in kit:   {Nocall_count} / {Nocalls_Percentage}%' )
        print( f'# Duplicate positions in kit: {duplicates_count} / {duplicates_percentage}%' )
        print( f'#')
        print( '#' * fenceNr)
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
            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )
            print()


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

            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )

            chromosome_count = chromosome_count + 1

            print()


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

            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )

            print()


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

            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )

            print()


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

            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )

            print()


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

            # Unique genotypes
            genotypeCount = filtered_df.groupby('genotype').size()
            print( f'Unique genotypes and their occurance' )
            print( genotypeCount )

            print()


        print()
        # Let user know processing is completed successfully
        print( 'Done analyzing file: ' + file.replace( inputFileDir, '' ) )
        print()
        print()
        print()

    # If file is unknown
    else:
        print()
        print( 'File ' + file.replace( inputFileDir, '' ) + ' is unknown' )
        print()


########################
# Compare overlapping SNPs

countCompanies = len(companyList)
companyIndex = 0

if countCompanies > 1:
    print( '#' * fenceNr )
    print( '#')
    print( '# Compare nr of SNPs that overlaps between companies' )
    print()

    # Get nr of companies analysed
    print( f'Nr of DNA files: {countCompanies}' )
    print()


    # Compare SNPs between companies
    for c in companyList:

        # Print current company and set index for dataframes to 0
        print( f'Nr {companyIndex}: {companyList[companyIndex]} overlap:' )
        dataFrameIndex = 0

        # Compare each dataframe in the dataframeList, not comparing companies to themselves
        for d in dataframeList:
            if companyIndex != dataFrameIndex:
                matchingCount, matchingCountX = compareDataframes( dataframeList[companyIndex], d )
                matchingCountAll = matchingCount + matchingCountX
                print( f'{companyList[companyIndex]} and {companyList[dataFrameIndex]}: 1-22 {matchingCount} SNPs, X {matchingCountX} SNPs. Sum: {matchingCountAll}' )

            dataFrameIndex = dataFrameIndex + 1
        
        companyIndex = companyIndex + 1
        print()

    print( '#' * fenceNr )
    print()



####################################################################################
# EOF #
####################################################################################
