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



## How it works
1. The script will determine what company that are used based the filename/comments. This is far from perfect.
2. It will then "normalize" the testkit to a standard format.
3. The gender of the kit will also be guessed, since it changes how the scrip handles X/Y/MT chromosomes (males have only one X and Y chromosome and cannot have heterozygous calls on these chromosomes)
4. If the kit is determined to be of male origin, then it will change heterozygous calls to nocalls.
5. The file will be somewhat cleaned by removing genotypes larger than two alleles and move calls on position 0 to "junk" chromosome 0.
6. Then it will concatenate the dna files and sort according a predetermined order.
7. Lastly it will format the superkit to the desired format and save it to ./output/




## How to use create_superkit.py:
This script will combine kits from compatible versions described above to a superkit

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
This script will analyse kits from compatible versions described above.
It will present the following data:
* Total SNPs
* Assumed gender of kit
* Unique chromosomes
* Unique genotypes
* SNP ranges

1. Put your raw autosomal DNA files into the ./input/ folder and make sure they are unpacked (in .txt or .csv format)

2. run python analyse_dna_file.py and the program will parse the DNA files in the default directory ./input/ and merge them together to a "SuperKit"


## TODO list
### Superkit Creator
- [ ] Add comments on top of superkit file
- [ ] Add support for the latest MyHeritage file format
- [ ] Add support for the latest TellMeGen file format
- [ ] Add support for command line arguments?
- [ ] Improve company detection "algorithm"
- [ ] Improve the genotype count per company produced in the end
- [ ] Improve each ouput format to more closely be like the originals
- [ ] Add algorithm to decide which duplicate is more correct by "majority count"

### DNA File Analyzer
- [X] Add information about what SNP span on each chromosome each company tests
- [ ] Add support for the latest MyHeritage file format
- [ ] Add support for the latest TellMeGen file format

### Other
- [X] Improve README!
