# DNA-autosomal-superkit-creator

This is a autosomal DNA super kit creator, based on this video:
https://www.youtube.com/watch?v=IJmAHNSODuw

This program will merge your DNA raw datafiles from different companies to one "super kit".
A super kit will supposedly increase overlap in matches in DNA matching sites like https://www.gedmatch.com

This script does this by searching ./input/ for .csv or .txt files and then screen the files for information about what test company that produced the file.
If the testing company is known, then the script cleans and standardizes and finally concatenates the data.
It will then sort the data by chromosome, position and then testcompany.
Lastly it will drop duplicates according to the priority list and format the data to the desired outputformat.
The data are then saved to the ./output/ folder for you to use.


Currently supported companies:
23andMe (V5)
Ancestry
FamilyTreeDNA
MyHeritage (Old format)
LivingDNA (v1.0.2)

Requirements: Python and Pandas

How to use:

1. Put your raw autosomal DNA files into the ./input/ folder and make sure they are unpacked (in .txt or .csv format)

2. Open main.py with a text editor and change your prefered options under "Customizations"
Currently, the only things you can change are the following:

* Output format
Changes the format to fit the format used by one of these companies:
Superkit (a kind of standardized format, close to 23andMe format)
23andMe (V5)
Ancestry
FamilyTreeDNA
MyHeritage (Old format)
LivingDNA v1.0.2

* Company priority
Decides which company genome that are preferred if several companies has called the same position in the same chromosome (calls get priority over nocalls by default

3. run python main.py and the program will parse the DNA files in the default directory ./input/ and merge them together to a "SuperKit"


TODO list
* Add comments on top of superkit file
* Improve company detection "algorithm"
* Improve the genotype count produced in the end