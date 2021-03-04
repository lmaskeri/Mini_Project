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
## Scripts Included
* **mini_wrapper.py** 
* **sleuth_script.R** 

## MiniProject.log Output Details
1. Number of CDS in the HCMV genome (EF999921).
2. Signficicant (FDR < 0.05) differentially expressed genes between the two timepoints (2pi and 6pi).
3. Number of reads in each transcriptome before and after Bowtie2 mapping for each Donor and timepoint.
4. The number of contigs with a length > 1000 in the SPAdes assembly.
5. The length of the SPAdes assembly (the total number of bp in all of the contigs > 1000 bp in length).
6. The top 10 blastn hits from using the longest contig from the SPAdes assmebly as query against a local database of sequences from the Betaherpesvirinae subfamily.

## Software Required
* **Kallisto**
  - https://pachterlab.github.io/kallisto/manual
* **Bowtie2**
  - http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
* **SPAdes**
  - http://cab.spbu.ru/files/release3.12.0/manual.html

### Packages Imported for Python
* **os** 
  - https://docs.python.org/3/library/os.html
* **argparse**
  - https://docs.python.org/3/library/argparse.html
* **Biopython**
  - https://biopython.org/wiki/Documentation
  - From Biopython:
    - **Entrez**
      - https://biopython.org/docs/1.75/api/Bio.Entrez.html
    - **SeqIO**
      - https://biopython.org/docs/1.75/api/Bio.SeqIO.html

### Packages Imported for R
* **sleuth** 
  - https://www.rdocumentation.org/packages/sleuth/versions/0.27.3
* **dplyr**
  - https://www.rdocumentation.org/packages/dplyr/versions/0.7.8















