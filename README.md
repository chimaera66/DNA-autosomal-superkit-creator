# DNA autosomal "Super Kit" creator

## Description:
This is a autosomal DNA super kit creator, based on this video by [Family History Fanatics & Genealogists](https://www.youtube.com/watch?v=IJmAHNSODuw).
More information can be found at [Do GEDmatch Superkits Give You More DNA Matches Than Non-Superkits?](https://www.familyhistoryfanatics.com/gedmatch-superkits)

This program will merge your DNA raw datafiles from different companies to one "super kit".
A super kit will supposedly increase overlap in matches and reduce false matches in DNA matching sites like [GEDMatch](https://www.gedmatch.com).

This script does this by searching `./input/` for `.csv` or `.txt` files and then screen the files for information about what test company that produced the file.
If the testing company is known, then the script cleans and standardizes and finally concatenates the data.
It will then sort the data by chromosome, position and then testcompany.
Lastly it will drop duplicates according to the priority list and format the data to the desired outputformat.
The data are then saved to the `./output/` folder for you to use.



## Currently supported companies (for both the creator and the analyser script):
* 23andMe v5 (since 08/2017)
* Ancestry v2 (since 05/2016)
* FamilyTreeDNA v3 (since 04/2019)
* MyHeritage v1 (before 03/2019)
* LivingDNA v1.0.2 (since 10/2018)



## Requirements:
* Python
* Pandas



## How to use create_superkit.py:

1. Put your raw autosomal DNA files into the ./input/ folder and make sure they are unpacked (in .txt or .csv format)

2. Open create_superkit.py with a text editor and change your prefered options under "Customizations"
Currently, the only things you can change are the following:

    * Output format
        - Changes the format to fit the format used by one of these companies:
            - Superkit (a kind of standardized format, close to 23andMe format)
            - 23andMe v5 (since 08/2017)
            - Ancestry v2 (since 05/2016)
            - FamilyTreeDNA v3 (since 04/2019)
            - MyHeritage v1 (before 03/2019)
            - LivingDNA v1.0.2 (since 10/2018)

    * Company priority
        - Decides which company genome that are preferred if several companies has called the same position in the same chromosome (calls get priority over nocalls by default

3. run python create_superkit.py and the program will parse the DNA files in the default directory ./input/ and merge them together to a "SuperKit"


## How to use analyse_dna_file.py:

1. Put your raw autosomal DNA files into the ./input/ folder and make sure they are unpacked (in .txt or .csv format)

2. run python analyse_dna_file.py and the program will parse the DNA files in the default directory ./input/ and merge them together to a "SuperKit"


## TODO list
- [ ] Add comments on top of superkit file
- [ ] Add support for the latest MyHeritage file format
- [ ] Add support for command line arguments
- [ ] Improve company detection "algorithm"
- [ ] Improve the genotype count per company produced in the end
- [ ] Add algorithm to decide which duplicate is more correct by "majority count"
- [ ] Add information about what SNP span on each chromosome each company tests
- [ ] Improve README!
