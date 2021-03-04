# Computational Biology Mini Project

## Installation
In order to run this code from your working directory, use this git command to clone this repository to your workspace:
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


## Directions for Project:

1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following
transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by
constructing the path based on the SRR numbers for each of these samples).
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI
accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with
kallisto. You will need to extract the CDS features from the GenBank format. Write the following to your log file (replace # with the number of coding sequences in the HCMV genome): The HCMV genome (EF999921) has # CDS.
3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially
expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth. Write the following details for
each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item:
target_id test_stat pval qval
4. It has been proposed that HCMV disease and pathogenesis may be related to the genetic diversity of the virus
(Renzette et al. https://www.ncbi.nlm.nih.gov/pubmed/25154343/). Which publicly available strains are most similar to
these patient samples? To compare to other strains, we will assemble these transcriptome reads. We don’t expect
assembly to produce the entire genome, but enough to be useful in BLAST. Virus sequencing experiments often include
host DNAs. It is difficult to isolate the RNA of just the virus (as it only transcribes during infection of the host cell). Before
assembly, let’s make sure our reads map to HCMV. Using Bowtie2, create an index for HCMV (NCBI accession EF999921).
Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each
transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would
write to the log (numbers here are arbitrary):
Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
5. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write
the SPAdes command used to the log file.
6. Write Python code to calculate the number of contigs with a length > 1000 and write the # to the log file as follows
(replace # with the correct integer):
There are # contigs > 1000 bp in the assembly.
7. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in
length) and write this # to the log file as follows (replace # with the correct integer):
There are # bp in the assembly.
8. Write Python code to retrieve the longest contig from your SPAdes assembly. Use the longest contig as blast+ input to
query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily. You will need to make a local
database of just sequences from the Betaherpesvirinae subfamily. Identify the top 10 hits. For the top 10 hits, write the
following to your log file: Subject accession, Percent identity, Alignment length, Start of alignment in query, End of
alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
Include the following header row in the log file, followed by the top 10 hits, and tab-delimit each item:
sacc pident length qstart qend sstart send bitscore evalue stitle 











