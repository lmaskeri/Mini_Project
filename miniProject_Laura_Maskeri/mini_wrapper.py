import os #importing os in order to access the command line through python
from Bio import Entrez #importing Entrez in order to retrieve fasta and genbank files from ncbi 
from Bio import SeqIO #importing SeqIO for parsing through fasta and fastq files as Seq Objects
import csv #importing csv to read in the blastn csv file 
import argparse



'''
1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following
transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by
constructing the path based on the SRR numbers for each of these samples).

https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1

'''


def get_paired_end_reads(SRR_list): #method to download reads from SRA and split into paired end read files
    for srr in SRR_list:
        path = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/" + srr + "/" + srr + ".1"
        os.system("wget " + path)   #wget from ncbi to get sra files
        os.system("fastq-dump -I --split-files " + srr + ".1") #fastq dump system call

   
##NOTE THE METHOD BELOW WAS USED TO CREATE THE SMALLER TEST FASTQ FILES UPLOADED TO GITHUB 
    ##I WILL NOT BE RUNNING THIS METHOD IN MY PIPELINE FOR THIS REASON
def split_fastq(srr):
    pair_1 = open(srr + ".1_1.fastq")
    pair_1_rec = SeqIO.parse(pair_1, 'fastq')
    count = 0
    with open(srr + "_small.1_1.fastq", "w") as small1:
        for rec in pair_1_rec:
            if count < 20000:
                small1.write(rec.format('fastq'))
                # small1.write(rec.format('fastq'))
                count += 1
            else:
                break
    pair_2 = open(srr + ".1_2.fastq")
    pair_2_rec = SeqIO.parse(pair_2, 'fastq')
    count = 0
    with open(srr + "_small.1_2.fastq", "w") as small2:
        for rec in pair_2_rec:
            if count < 20000:
                small2.write(rec.format('fastq'))
                # small1.write(rec.format('fastq'))
                count += 1
            else:
                break
            
'''
2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI
accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with
kallisto. You will need to extract the CDS features from the GenBank format. Write the following to your log file (replace
# with the number of coding sequences in the HCMV genome):
The HCMV genome (EF999921) has # CDS.

'''

                    
def generate_transcriptome():
    Entrez.email = "lmaskeri@luc.edu"
    # get the number of CDS sections from genbank file
    #I will use efetch with return type as a genbank record
    #So I can parse through the record and find all of the CDS sections present
    gen_handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text") #creating another entrez handle in order to get the genbank record (with all CDS sections) for this acc number
    #searching through the record for this acc number to get the # of CDS sections for the HCMV genome
    #as well as to put each CDS name and sequence into an outfile in fasta format -- this is the transcriptome file
    cds_count = 0 #counter variable to hold the number of cds 
    with open("hcmv_transcriptome.fasta", "w") as out:
        for rec in SeqIO.parse(gen_handle, "genbank"): #makes Seq generator object that holds all of the record components for this acc number
            for rec_feat in rec.features: #every CDS feature after the first 
                feat_type = rec_feat.type
                if feat_type == 'CDS':
                    #I want to find the exact sequence of NUCLEOTIDES (not the ["translation"])
                    #of each CDS. So first I need to pull out the location within the sequence that the CDS feature is found (from what nucelotide: to what nucelotide)
                    #and then extract the record AT that position in the sequence and grab the sequence itself using .seq (since it is still a SeqIO object)
                    cds_seq = rec_feat.location.extract(rec).seq
                    #pull out the name of the specific CDS using .qualifiers
                    cds_name = str(rec_feat.qualifiers["protein_id"]) #pull out protein id name
                    cds_count += 1 #increase count for the number of cds
                    out.write(">" + cds_name[2:-2] + "\n") #write to the transcriptome file excluding the [' '] by using splicing 
                    out.write(str(cds_seq) + "\n") #write the sequence out as well, excluding the [' ']

    with open("miniProject.log", "a") as out: #append to the log
        out.write("The HCMV genome (EF999921) has " + str(cds_count) + "CDS." + "\n")
     
    
'''
3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially
expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth. Write the following details for
each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item:
target_id test_stat pval qval
'''

def create_kallisto_index():
    #creating the kallisto index
    os.system(" time kallisto index -i hcmv_kallisto_index.idx hcmv_transcriptome.fasta")
  
    
def run_kallisto(srr):
    os.system(" time kallisto quant -i hcmv_kallisto_index.idx -o tpm_" + str(srr) + " -b 30 -t 2 " + str(srr) + ".1_1.fastq " + str(srr) + ".1_2.fastq")
    #this produces an abundance.h5 compressed file that contains:
        #quantifications, all bootsftraps, and other information
    #an abundance.tsv tab separated file that holds:
        #abundances and counts for each transcript
    #and a run_info.json formatted file that contains:
        #information about the run 
    #for each srr!
 

def make_sleuth_table(SRR_list):
    #Sleuth requires a table of your samples, conditions, and paths to kallisto results
    #You must generate a tab- or space-delimited .txt file using Unix commands, Python, or R.
    #using the kallisto reads, I need to make a txt table in order to use sleuth to find differentially expressed genes between 2pi and 6pi 
    #2dpi = SRR5660030, SRR5660044
    #6dpi = SRR5660033, SRR5660045
    #formatED each tab delimited line as name -tab- 2dpi_or_6pi -tab- location
    with open("sleuth_table.txt", "w") as out:
        out.write("sample" + "\t" + "condition" + "\t" + "path" + "\n") #writing the header 
        for srr in SRR_list: #for each srr
            file_dir = "./" + "tpm_" + srr #location where folder containing the files needed for sleuth will be (using the current directory)
            if srr == "SRR5660030" or srr == "SRR5660044":
                out.write(srr + ".1" + "\t" + "2dpi" + "\t" + file_dir + "\n")
            if srr == "SRR5660033" or srr == "SRR5660045":
                out.write(srr + ".1" + "\t" + "6dpi" + "\t" + file_dir + "\n")

def run_sleuth():
    os.system("Rscript sleuth_script.R")
    with open("miniProject.log", "a") as out: #append to the log
        sleuth = open("sleuth_results.txt").read().split("\n")
        for s in sleuth:
            out.write(s + "\n")
            
    


'''
4. Using Bowtie2, create an index for HCMV (NCBI accession EF999921).
Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each
transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would
write to the log (numbers here are arbitrary):
Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read
pairs after.
'''

def get_bowtie2_index():
    #access entrez in order to write a fasta file out for HCMV acc EF999921
    #this fasta file is for use in bowtie 2 index
    Entrez.email = "lmaskeri@luc.edu"
    handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype='fasta')
    rec = handle.read()
    with open("hcmv_fasta.fasta", "a") as out:
        out.write(rec)  
    #obtaining the index using this fasta file
    os.system("bowtie2-build hcmv_fasta.fasta hcmv_bowtie")

def run_bowtie2(srr):
    pair_1 = srr + ".1_1.fastq"
    pair_2 = srr + ".1_2.fastq"
    #using --al-conc in order to get all mapped reads 
    #from bowtie documentation: "paired-end reads that align concordantly at least once"
    os.system("bowtie2 --quiet --al-conc mapped_bowtie_" + srr + ".fastq -x hcmv_bowtie -1 " + pair_1 + " -2 " + pair_2 + " -S " + srr + "map.sam")
    
def count_bowtie2_reads(srr):
    #count number of reads BEFORE bowtie within the fastq files
    pair_1 = open(srr + ".1_1.fastq")
    pair_1_rec = SeqIO.parse(pair_1, 'fastq')
    count1 = 0
    for rec in pair_1_rec:
        count1 += 1
    pair_1.close()
    pair_2 = open(srr + ".1_2.fastq")
    pair_2_rec = SeqIO.parse(pair_2, 'fastq')
    count2 = 0
    for rec in pair_2_rec:
        count2 += 1
    pair_2.close()
    total_reads_before = count1 + count2

    #now find the total reads after filtering using bowtie
    bow_pair_1 = open("mapped_bowtie_" + srr + ".1.fastq")
    bow_pair_1_rec = SeqIO.parse(bow_pair_1, 'fastq')
    bow_count1 = 0
    for rec in bow_pair_1_rec:
        bow_count1 += 1
    bow_pair_1.close()
    bow_pair_2 = open("mapped_bowtie_" + srr + ".2.fastq")
    bow_pair_2_rec = SeqIO.parse(bow_pair_2, 'fastq')
    bow_count2 = 0
    for rec in bow_pair_2_rec:
        bow_count2 += 1
    bow_pair_2.close()
    total_reads_after = bow_count1 + bow_count2

    #determine which donor it was and dpi condition to write out to 
    donor = ''
    if srr == "SRR5660030":
        donor = "Donor 1 (2dpi)"
    if srr == "SRR5660033":
        donor = "Donor 1 (6dpi)"
    if srr == "SRR5660044":
        donor = "Donor 3 (2dpi)"
    else:
        donor = "Donor 3 (6dpi)"
        
    with open("miniProject.log", "a") as out: #append to the log
        out.write(donor + " had " + str(int(total_reads_before)) + " read pairs before Bowtie2 filtering and " + str(total_reads_after) + " read pairs after." + "\n")
            

'''
5. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write
the SPAdes command used to the log file.
'''


def SPAdes(SRR_list):

    srr1, srr2, srr3, srr4 = SRR_list[0],SRR_list[1],SRR_list[2],SRR_list[3] #pulling out all srr's from the list to use in the spades command
    #Using this format for multiple read pairs in the spades command
    #--pe1-1 myreads_1.fastq –-pe1-2 myreads_2.fastq --pe2-1 morereads_1.fastq –-pe2-2 morereads_2.fastq
    #Don't have to use all 4 kmer sizes for final run -- just keep consistant
    #For this script, I am just using k sizes of 55 and 127!
    spades_command = " spades -k 55, 127 -t 2 --only-assembler --pe1-1 mapped_bowtie_" + srr1 + ".1.fastq --pe1-2 mapped_bowtie_" + srr1 + ".2.fastq --pe2-1 mapped_bowtie_" + srr2 + ".1.fastq --pe2-2 mapped_bowtie_" + srr2 + ".2.fastq --pe3-1 mapped_bowtie_" + srr3 + ".1.fastq --pe3-2 mapped_bowtie_" + srr3 + ".2.fastq --pe4-1 mapped_bowtie_" + srr4 + ".1.fastq --pe4-2 mapped_bowtie_" + srr4 + ".2.fastq -o spades_assembly/"
    os.system(spades_command)
    with open("miniProject.log", "a") as out: #append to the log
        out.write(spades_command + "\n")
            


'''     
6. Write Python code to calculate the number of contigs with a length > 1000 and write the # to the log file as follows
(replace # with the correct integer):
There are # contigs > 1000 bp in the assembly.    
'''

def count_contigs():
    contigs_file = open("./spades_assembly/contigs.fasta") #opening the fasta file using the wd and the folder name created by SPAdes
    contigs = SeqIO.parse(contigs_file, "fasta") #parsing the fasta file
    count = 0
    with open("longest_contigs_assembly.fasta" , "w") as out: #opening a new fasta file to write the longest contigs over 1000 in length to for the next step
        for contig in contigs:
            if len(contig.seq) > 1000: #if the length of the contig sequence is over 1000
                count += 1 #add one to the count of contigs over 1000 bp long
                out.write(contig.format('fasta')) #write that contig's information out to a new fasta file in fasta format
    with open("miniProject.log", "a") as out: #append to the log
        out.write("There are " + str(count) + " contigs > 1000 bp in the assembly" + "\n")

    
'''
7. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in
length) and write this # to the log file as follows (replace # with the correct integer):
There are # bp in the assembly.
 '''

def length_of_assembly():
    #open back up the longest_contig_assembly.fasta file to calculate the length of the contigs that were > 1000 bp
    contigs_file = open("./longest_contigs_assembly.fasta") #from the current directory, open the fasta file
    contigs = SeqIO.parse(contigs_file, "fasta")
    assembly_length = 0 #variable to hold the total number of bp from ALL of the contigs that were > 1000 bp in length
    for contig in contigs:
        assembly_length = assembly_length + len(contig.seq) #add the # of bp from each contig using the len() function
    with open("miniProject.log","a") as out: #append to the log
        out.write("There are " + str(assembly_length) + " bp in the assembly." + "\n")


'''
8. Write Python code to retrieve the longest contig from your SPAdes assembly. Use the longest contig as blast+ input to
query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily. You will need to make a local
database of just sequences from the Betaherpesvirinae subfamily. Identify the top 10 hits. For the top 10 hits, write the
following to your log file: Subject accession, Percent identity, Alignment length, Start of alignment in query, End of
alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
Include the following header row in the log file, followed by the top 10 hits, and tab-delimit each item:
sacc pident length qstart qend sstart send bitscore evalue stitle 
'''

def longest_contig():
    #retrieving the longest contig from the SPAdes assembly
    contigs_file = open("./longest_contigs_assembly.fasta") #from the current directory, open the fasta file
    contigs = SeqIO.parse(contigs_file, "fasta")
    longest = 0
    longest_contig = ''
    for contig in contigs:
        if len(contig.seq) > longest:
            longest = len(contig.seq)
            longest_contig = contig
    with open("longest_contig.fasta" , "w") as out: #writing the longest contig out in fasta format in order to blast it!
        out.write(longest_contig.format("fasta"))


def make_local_blast_database():
    ##NOTE: I CREATED THIS DATABASE WITH THE DOWNLOADED FASTA FILE WITH ALL SEQS FROM THE BETAHERPESVIRINAE SUBFAMILY (USING NCBI)
        ##HOWEVER, I DID NOT UPLOAD THIS FASTA FILE TO MY REPO BECAUSE IT IS VERY LARGE
        ##HERE IS THE LINK THAT I USED TO IDENTIFY THE FASTA SEQUENCES WITH "txid10357[Organism:exp]" : https://www.ncbi.nlm.nih.gov/nuccore
    #make local database of just sequences from the Betaherpesvirinae subfamily
    os.system("makeblastdb -in Betaherpesvirinae_refseqs.fasta -out blast_db -title blast_db -dbtype nucl")
    
def blast():
    #BLAST!!
    os.system('blastn -query longest_contig.fasta -db blast_db -out blastn_results.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"')
    
def parse_blast_results():
    #now parse through the csv file to pull out the top 10 hits 
    #and print them to the log file
    headers = ["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"]
    blast_results = open("blastn_results.csv", "r")
    rows = csv.DictReader(blast_results, headers, delimiter = ",") #using DictReader to read in each row of csv as dictionary (where header columns are the key and the values are the results for each)
    count = 0
    with open("miniProject.log", "a") as out: #append top 10 hits to the log
        for row in rows:
            if count >= 9: #count used to get only top 10 hits
                break
            else:
                #since each row is a dictionary, I will pull out all values for each key 
                #then write all out tab delimited
                out1 = str(row["sacc"])
                out2 = str(row["pident"])
                out3 = str(row["length"]) 
                out4 = str(row["qstart"]) 
                out5 = str(row["qend"])
                out6 = str(row["sstart"])
                out7 = str(row["send"])
                out8 = str(row["bitscore"]) 
                out9 = str(row["evalue"]) 
                out10 = str(row["stitle"])
                out.write(out1 + "\t" + out2 + "\t" +  out3 + "\t" +  out4 + "\t" +  out5 + "\t" +  out6 + "\t" +  out7 + "\t" +  out8 + "\t" +  out9 + "\t" +  out10 + "\n")
            count += 1


##LIST OF SRRs TO BE USED IN THE PYTHON WRAPPER
SRR = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]


###  RUN SCRIPT ###

#First, setting up argparser in order to take in the argument from the command line 
#to either run the entire dataset or to run the test dataset
parser = argparse.ArgumentParser(description = "Determine which dataset to run through script.")
#adding an argument for the specific flag needed to determine which dataset the user wants to run
parser.add_argument('dataset_to_run', type = str, help = "Specify which dataset you'd like to run: test_dataset or whole_dataset")
#parsing the arguments
argument = parser.parse_args()
#args.dataset_to_run now holds the string from the commandline that will determine if the whole dataset is run or the test dataset is run

#checking the flag with a conditional

if argument.dataset_to_run == "test_dataset":
    #skip step 1 because paired read fastq files are already in github repo (that should be cloned to the user's directory)
    #first, generate the HCMV transcriptome
    generate_transcriptome()
    #First build the kallisto index...
    create_kallisto_index()
    #Then run kallisto using the index and quantify the TPM of each CDS in each SRR transcriptome
    for srr in SRR:
        run_kallisto(srr)
    #using the generated tpm_srr folders, 
    #find differentially expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth
    #First create the sleuth table for input...
    make_sleuth_table(SRR)
    #Then run sleuth using R
    run_sleuth()
    #Using Bowtie2, create an index for HCMV (NCBI accession EF999921).
    get_bowtie2_index()
    #Then run bowtie2 to filter read pairs for each SRR
    for srr in SRR:
        run_bowtie2(srr)
    #count the bowtie2 reads to determine how many read pairs before and after filtering
    for srr in SRR:
        count_bowtie2_reads(srr)
    #assemble all four transcriptomes together to produce 1 assembly using SPAdes
    SPAdes(SRR)
    #calculate the number of contigs with a length > 1000 for the assembly
    count_contigs()
    #calculate the length of the assembly
    length_of_assembly()
    #retrieve the longest contig from your SPAdes assembly
    longest_contig()
    #Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily
    #make_local_blast_database() STEP OMITTED DUE TO THE BLAST DATABASE ALREADY BEING ON MY REPO
    blast() 
    #pull out top ten hits from blast
    parse_blast_results()
    
if argument.dataset_to_run == "whole_dataset":
    #do all steps above, however add in the first step
    #this is where you retrieve the transcriptomes from two patient donors from SRA and convert to paired-end fastq files
    get_paired_end_reads(SRR)
    generate_transcriptome()
    create_kallisto_index()
    for srr in SRR:
        run_kallisto(srr)
    make_sleuth_table(SRR)
    run_sleuth()
    get_bowtie2_index()
    for srr in SRR:
        run_bowtie2(srr)
    for srr in SRR:
        count_bowtie2_reads(srr)
    SPAdes(SRR)
    count_contigs()
    length_of_assembly()
    longest_contig()
    blast()
    parse_blast_results()




                 
                    






































