# DNA autosomal "Super Kit" creator

## Description:
This is a autosomal DNA super kit creator, based on this video by [Family History Fanatics & Genealogists](https://www.youtube.com/watch?v=IJmAHNSODuw).
More information can be found at [Do GEDmatch Superkits Give You More DNA Matches Than Non-Superkits?](https://www.familyhistoryfanatics.com/gedmatch-superkits)

This script will merge your DNA raw datafiles from different companies to one "super kit".
A super kit will supposedly increase overlap in matches and reduce false matches in DNA matching sites like [GEDMatch](https://www.gedmatch.com).

The script can also act like a converter, to convert from one format to another. Use command line parameters -t -o "DESIREDFORMAT".
It is important to know that different companies tests different SNP range, and the script cannot create SNPs that doesn't alreade exist in your provided files, so it is possible that some part of the output formats SNP ranges are missing if you convert from another company that doesn't test these ranges.

This script works by searching `./input/` for `.csv` or `.txt` files and then screen the files for information about what test company that produced the file.
If the testing company is known, then the script cleans and standardizes and finally concatenates the data.
It will then sort the data by chromosome, position and then testcompany.
Lastly it will drop duplicates in two or three steps (depending on if you choose majority vote or not), trim SNP ranges (if you choose that option in the command line argument) and format the data to the desired outputformat.
The data are then saved to the `./output/` folder for you to use.



## Currently supported companies:
* 23andMe v5 (since 08/2017)
* Ancestry v2 (since 05/2016)
* FamilyTreeDNA v3 (since 04/2019)
* MyHeritage v1 (before 03/2019)
* MyHeritage v2 (after 03/2019)
* LivingDNA v1.0.2 (since 10/2018)



## Requirements:
* Python
* Pandas



## How it works
1. The script will determine what company that are used based the filename/comments. This is far from perfect.
2. It will then "normalize" the testkit to a standard format.
3. The gender of the kit will also be guessed, since it changes how the script handles X/Y/MT chromosomes (males only have one X and Y chromosome and cannot have heterozygous calls on these chromosomes)
4. If the kit is determined to be of male origin, then it will change heterozygous calls to nocalls.
5. The file will be somewhat cleaned by removing genotypes larger than two alleles and move calls on position 0 to "junk" chromosome 0.
6. Then it will concatenate the dna files and sort according to a predetermined order.
7. Lastly it will format the superkit to the desired format with correct top commments and save it to `./output/`

* At the moment, MyHeritage does not accept the resulting superkit, not in any output format




## How to use create_superkit.py:
This script will combine kits from compatible versions described above to a superkit

1. Put your raw autosomal DNA files into the `./input/` folder and make sure they are unpacked (in `.txt` or `.csv` format)

2. Open `create_superkit.py` with a text editor and change your prefered options under `Customizations`
Currently, the only things you can change are the following:

    * Company priority
        - Decides which company genome that are preferred if several companies has called the same position in the same chromosome (calls get priority over nocalls)

3. run python `create_superkit.py` and the program will parse the DNA files in the default directory `./input/` and merge them together to a `SuperKit`

4. Currently supported command line arguments are
    * -o, --outputFormat: Sets the template for the formatting of the output file. Valid formats are: SuperKit, "23andMe v5", "AncestryDNA v2", "FamilyTreeDNA v3", "LivingDNA v1.0.2", "MyHeritage v1" and "MyHeritage v2". Defaults to SuperKit.
    * -t, --trimSNP: Trims the SNPs to fit within the tested ranges of the different companys that the outputFormat tries to emulate. Defaults to false.
    * -mv, --majorityVote: Drops genotype based on a majority vote. If there are two AA and one CC on the same position, then one AA is kept and the other rows drops. This is considerably slower than the normal keep first row, but it should be more accurate. Mostly meaningful when merging three kits or more. Defaults to false.




## How to use analyse_dna_file.py:
This script will analyse kits from compatible versions described above.
It will present the following data:
* Total SNPs
* Assumed gender of kit
* Unique chromosomes
* Unique genotypes
* SNP ranges

1. Put your raw autosomal DNA files into the `./input/` folder and make sure they are unpacked (in `.txt` or `.csv` format)

2. run python `analyse_dna_file.py` and the program will parse the DNA files in the default directory `./input/`.


## TODO list
### Superkit Creator
- [X] Add support for the latest MyHeritage file format
- [ ] Add support for the latest TellMeGen file format
- [ ] Improve company detection "algorithm"
- [ ] Add correct filenames for the superkits?


### DNA File Analyzer
- [X] Add support for the latest MyHeritage file format
- [ ] Add support for the latest TellMeGen file format
