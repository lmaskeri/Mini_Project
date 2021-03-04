# Computational Biology Mini Project

## Installation
In order to run this code from your working directory, use this git command to clone this repository to your workspace.
```
git clone https://github.com/lmaskeri/Mini_Project.git
```
## Directions
You can run this code with either the whole dataset, or the test dataset. 
Running the test dataset will use the smaller fastq files cloned from the repository instead of downloading and creating them within the wrapper. To run the test dataset, use the following command:
```
python3 mini_wrapper.py test_dataset
```

To run the whole dataset with all paired reads, use the following command:
```
python2 mini_wrapper.py whole_dataset
```

## Software Required
### Packages Imported for Python
**-os** 
https://docs.python.org/3/library/os.html
**-Entrez**
https://biopython.org/docs/1.75/api/Bio.Entrez.html
**-SeqIO**
https://biopython.org/docs/1.75/api/Bio.SeqIO.html
**-argparse**
https://docs.python.org/3/library/argparse.html












