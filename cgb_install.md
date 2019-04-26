# Installing CGB

CGB is distributed as a conda package. The following instructions refer to the installation of CGB using conda.

## Installing conda
To install `miniconda` go to its [website](https://conda.io/docs/user-guide/install/index.html) and [download](https://conda.io/miniconda.html) the appropriate` Python 2.7` miniconda installer version [Linux will be assumed in the following].

For Linux, open up a terminal and execute the bash installer:
`bash Miniconda2-latest-Linux-x86_64`
or
`bash Miniconda2-latest-Linux-x86.sh`

(allow it to prepend to PATH the Miniconda2 install location; close terminal and reopen for changes to go into effect)

## Creating a conda environment
Create the environment
`conda create -n cgb python=2.7`

Activate the environment
`source activate cgb`

Any changes made from now on (while environment is active) will be made only to this environment.

## Setting up the conda environment
Install packages:
- `conda install biopython`
- `conda install cached-property==1.3.0`
- `conda install cycler==0.10.0`
- `conda install scipy==0.17.1`
- `conda install matplotlib==1.5.1`
- `conda install pyparsing==2.1.4`
- `conda install python-dateutil==2.5.3`
- `conda install six==1.10.0`
- `conda install pytz==2016.4`
- `conda install -c conda-forge tqdm`
- `conda install --channel https://conda.anaconda.org/etetoolkit ete3==3.0.0b35`
- `conda install -c bioconda weblogo`
- `conda install networkx`

## Installing git
If you have not installed `git` in your system, install git
`sudo apt-get install git`

Link to your [github](http://www.github.com) account:
- `git config --global user.name "yourusername"`
- `git config --global user.email "youremailaddress"`


### Creating keys
- `mkdir .ssh`
- `ssh-keygen -t rsa -C "youremailaddress"`
- `eval "$(ssh-agent -s)"`
- `ssh-add ~/.ssh/id_rsa`

Open `~/.ssh/id_rsa.pub` file with a text editor and copy the contents of the file manually into [GitHub](https://github.com/settings/keys) --> Account settings --> SSH Keys --> New SSH Key

Test the key:
`ssh -T git@github.com`

Make folder for repo and move there. Then, initialize git and link to and clone GitHub repo:
- `git init`
- `git remote add origin https://github.com/ErillLab/cgb.git`
- `git clone https://github.com/ErillLab/cgb.git`

## Dependencies
CGB requires three external programs to be installed: BLAST, CLUSTALO and BayesTraits. Optionally, HMMER can also be installed to enable querying of the eggNOG database. The following instructions are for installing these programs on Linux.

### BLAST
On a terminal:
`sudo apt-get install ncbi-blast+`

Check that BLAST is properly installed by asking:
`which blastp`
which should return:
`/usr/bin/blastp`

Create a folder where to store BLAST databases (e.g. under your home)
`mkdir BLASTdbs`
and export the environment variable
`export BLASTDB=/home/yourusername/BLASTdbs`

### CLUSTALO
Install simply through:
`sudo apt-get install clustalo`

Check that CLUSTALO is properly installed by asking:
`which clustalo`
which should return:
`/usr/bin/clustalo`

### BayesTrait
BayesTrait should be already present in the `cgb/bin` folder after cloning. If using a 32-bit machine, [download](http://www.evolution.rdg.ac.uk/BayesTraitsV2Beta.html) the BayesTrait 32-bit executable and replace the file `BayesTraitsV2_linux` in `cgb/bin` with the donwloaded executable.

### HMMER
Install simply through:
`sudo apt-get install hmmer`

Check that HMMER is properly installed by asking:
`which hmmsearch`
which should return:
`/usr/bin/hmmsearch`

#### eggNOG
In order to query eggNOG with HMMER to assign NOGs to orthologous groups, you'll need to first download the current version of the [eggNOG database](http://eggnogdb.embl.de/#/app/downloads). You can download any NOG compilation that you deem applicable to your analyses, but for generality the [BactNog](http://eggnogdb.embl.de/download/eggnog_4.5/data/bactNOG/bactNOG.hmm.tar.gz) compilation is recommended.

Donwload and uncompress the ~3 Gb BactNog tarball containing the bacterial NOG Hiden Markov Models (.hmm).
`tar xvzf bactNOG.hmm.tar.gz`
(this can take a while)
(you can remove the bactNOG.hmm.tar.gz file after this)

Then concatenate the extracted files into a single file with:
`for f in *.hmm; do cat "$f" >> bact.hmmer; done`
(this can take _quite_ a while)
(you can delete all  the .hmm files after this)

And use this file to create a HMMER database with:
`hmmpress bact.hmmer`
(this can take a while)
(you can delete the bact.hmmer file after this)

#### PFAM
In order to query PFAM with HMMER to assign PFAM IDs to orthologous groups, you'll need to first download the current version of the [PFAM database](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/). You can download any PFAM compilation that you deem applicable to your analyses, but the latest [PFAM full database](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.full.gz) is recommended.

Donwload and uncompress the ~10 Gb PFAM zip file containing the bacterial NOG Hiden Markov Models (.hmm).
`gzip -d Pfam-A.full.gz`

And use this file to create a HMMER database with:
`hmmpress Pfam-A.full`
