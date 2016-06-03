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
And then copy the Bauer Center Raw file to the HTS directory

To Run the script:
* Start your ssh session
* Goto the Directory with your required files
* run: "source analysis.sh "BAUER_CENTER_SEQUENCING_FILE" "
* wait around 6/7 hours (You will receive an email once is over)

The output of the script will be in mochiview compatible *.wig files in the mochiview directory
There are also two pdf containing some statistics of the script. If something went wrong send me the log.out file so that I can take a look
