# Instructions to use HTS pipeline
This pipeline works for raw files coming from the Bauer center sequencing facility.
It works in the HMS Orchestra cluster as of 10/20/2015

Make sure to have in the directory where you plan to analyze the results, aside from the sequencing file,
the following files:
* analysis.sh
* analysis.py
* log_parser.py
* BDGtoWIG.py
* fastq_trimmer.py
* a setup.cfg file (check example file)
** In the setup file don't change the first 3 lines unless you know what you're doing, the next lines should contain the name you want to give your sample
followed by the Input barcode name (eg:BAR40) followed by the IP barcode name, each one separated by a tab. If you don't have input samples or you have very few reads on the input just put as input bar code the same as your sample barcode (the peak calling won't give you any meaningfull results but without inputs you shouldn't be able to peak call, but you would still have the track pileup)

The easiest way to have all the files is going to your orchestra directory and running:
git clone https://github.com/LuisSoares/HTS.git
And then copy the Bauer Center Raw file to the HTS directory.

Update July 2016, the fastest way is now to use the sabre branch of the git repository, using git clone -b sabre https://github.com/LuisSoares/HTS.git

Sugestion: Personal folders in orchestra are restricted to 100GB, at one point the analysis will need more than 30GB so it is better that you run the analysis in a folder in the temporary /n/scratch2 filesystem where you have several terabytes available (beware that scratch is only temporary!).

To Run the script:
* Start your ssh session
* Goto the Directory with your required files
* run: "source analysis.sh "BAUER_CENTER_SEQUENCING_FILE" "
* wait around 6/7 hours (You will receive an email once is over) (less than 1 hour if sabre branch is used)

The script will output the following:

* mochiview directory with wig files ready to load in mochiview
* fastqc directory with one html and one zip file for quality controls for which barcode
* temp directory containing the bam files for each barcode as well as the peak calling files for each barcode
* demultiplex and mapping pdf files with statistics of these two steps
* log.out file with the entire log of the script run

