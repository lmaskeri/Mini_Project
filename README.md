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
* **Kallisto**
  - Manual: https://pachterlab.github.io/kallisto/manual
* **Bowtie2**
  - Manual: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
* **SPAdes**
  - Manual: http://cab.spbu.ru/files/release3.12.0/manual.html

### Packages Imported for Python
* **os** 
  - https://docs.python.org/3/library/os.html
* **argparse**
  - https://docs.python.org/3/library/argparse.html
* **Biopython**
  - https://biopython.org/wiki/Documentation
  - From Biopython:
  * **Entrez**
    - https://biopython.org/docs/1.75/api/Bio.Entrez.html
  * **SeqIO**
    - https://biopython.org/docs/1.75/api/Bio.SeqIO.html

### Packages Imported for R
* **sleuth** 
  - https://www.rdocumentation.org/packages/sleuth/versions/0.27.3
* **dplyr**
  - https://www.rdocumentation.org/packages/dplyr/versions/0.7.8













