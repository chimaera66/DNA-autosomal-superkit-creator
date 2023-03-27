


import os                   # For findDNAFiles
from typing import List
import pandas as pd
import re                   # For determineDNACompany



# Input/output file directory
inputFileDir = './input/'
outputFileDir = './output/'
# Output File name and file ending
#outputFileName = 'DNASuperKit'
#outputFileEnding = '.csv'

# Input filetypes
fileEndings = (
    'txt',
    'csv'
)


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
# Pre-screen file to determine DNA company
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
# Try to determine what DNA testing
# company the file originates from
#

def determineDNACompany(text, filename):

    company_patterns = {
        '23andMe V5': r'_v5_Full_',
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

def prepareDNAFile( f, c ):

#    print ( c )
    # 23andMe and Living DNA
    if c == '23andMe v5' or c == 'LivingDNA v1.0.2':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, sep='\t', comment='#', index_col=False, header=None, engine='python' )

    # MyHeritage (Before 1 March, 2019) and FamilyTreeDNA v3
    elif c == 'MyHeritage v1' or c == 'FamilyTreeDNA v3':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, comment='#' )

    # AncestryDNA v2
    elif c == 'AncestryDNA v2':
        # Load input file into pandas and create the proper columns
        df = pd.read_csv( f, dtype=str, sep='\t', comment='#' )

        # Normalize chromosome order with custom chromosomeTable
        #df[ 'chromosome' ].replace( to_replace=chromosomeTableAncestryIn, inplace=True )

        # Merge allele1 and allele2 to genotype column
        df[ 'genotype' ] = df[ 'allele1' ] + df[ 'allele2' ]
        del df[ 'allele1' ]
        del df[ 'allele2' ]

    # Normalize columnnames
    df.columns = [ 'rsid', 'chromosome', 'position', 'genotype' ]

    # and add company column
    df[ 'company' ] = c

# DEBUG ?
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

for f in rawDNAFiles:
    # Screening file to determine company
    fileScreening = prescreenDNAFile( f )

    # Get DNA company from file comment
    company = determineDNACompany( fileScreening , f)
#    print( company )


    if company != 'unknown':
        print( 'Processing file:')
        print( f.replace( inputFileDir, '' ) )
        print( 'Which contains data from ' + company )
        
        # Normalize and clean DNA file in preparation for concatenation
        df = prepareDNAFile( f, company )

        # Handle FamilyTreeDNA v3 a bit differently
        # delete all data in chromosome 0 (nocalls? bad data?)
#        if company == 'FamilyTreeDNA v3':
#            chromosomeZero = df.loc[ df[ 'chromosome' ] == '0' ]
#            del chromosomeZero[ 'company' ]

        # Clean dataframe
#        df = cleanDNAFile( df, company )
        
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


####################################################################################
# EOF #
####################################################################################