hBayeSSC is a python script that wraps around serial simcoal in order to simulate a multi-taxa community undergoing a coordinated demographic expansion.

## Requirements
The files needed to produce a set of simulations with multitaxa summary statistics for the hABC analysis described in Chan et al. 2014 are: 
 *  BayeSSC - The executable for Serial Simcoal which can be found here: [site](http://www.stanford.edu/group/hadlylab/ssc/)   
 *  python >= 2.4   
 *  hBayeSSC.py
 *  msReject (Optional) [Install Instructions](#msreject-module)

## Input files
 *  A table of observed summary statistics for each taxon in the community. ([details](#observation-summary-statistics))
 *  An input par file for serial simcoal. ([details](#par-file))

### Observation summary statistics
[Sample file](/example_data/example_obs)  

The table of observed summary statistics consists of columns with the following header names.  hBayeSSC replaces the appropriate line in the par file with these values:

| Column name | Description |
| ----- | ----------- |
| species | name of taxa |
| nsam | number of samples to be simulated |
| nsites | number of base pairs |
| tstv  | % transitions |
| gamma |  gamma shape parameter |
| gen | numbers of years per generation |
| locuslow | low estimate of the locus mutation rate per generation |
| locushigh | high estimate of the locus mutation rate per generation |
| Nelow | low estimate for effective population size |
| Nehigh | high estimate for effective population size |
| SegSites | segregating sites |
| nucdiv | nucleotide diversity |
| Haptypes | number of haplotypes |
| HapDiver | haplotypic diversity |
| TajimasD | Tajima's D |
| F* | Fu's F |

A more complete description of these values can be found on the [BayeSSC website](http://www.stanford.edu/group/hadlylab/ssc/)

### par file
[Sample file](/example_data/example.par)  

The par file contains one prior which is not individually replaced, such as expansion magnitude (under historical events) and will apply to all populations.

## Usage

An example command to run 200 iterations of each model from 0/32 to 32/32 species with a prior on the maximum expansion time of 1000 to 500000 years would be: 
```
python hBayeSSC.py --mode initial -p example.par -i example_obs -r 200 -u full -b ./BayeSSC -t 1000:500000
```

This command will create the observed hyperstats file for the rejection analysis: 
```
python hBayeSSC.py --mode initial -p example.par -i example_obs -r 200 -u full -b ./BayeSSC -t 1000:500000 --obs_hyperstats
```

## Options
The hBayeSSC has several command line options, which can be found using the -h option when executing the script.  

```  
#> python hBayeSSC.py -h

Usage: hBayeSSC.py [options]

Options:
  -h, --help            show this help message and exit
  --mode=MODE           program operation mode [ 'initial', 'posterior' ]
                        [required]
  -p FILE, --par=FILE   par file template [required]
  -i FILE, --obs=FILE   Observation file [required]
  -r NUM, --repeat=NUM  Number of times to try a given congruent group size
                        [required]
  -u UID, --uid=UID     Unique ID to prefix generated indices [required]
  -b PATH, --bayepath=PATH
                        Path to BayeSSC application [default: Located on user
                        PATH]
  --only_hyperstats     When set, will only generate the hyperstats file
  --print_headers       When set will generate a headers.txt and exit
  -o PATH, --outdir=PATH
                        Directory to generate final outputs in (will create
                        missing folders) [default: <working directory> ]

  Regular Run:
    Options to be applied during mode 'initial'

    -m MODEL, --model=MODEL
                        Run a single model (0 to total entries in observation
                        file) [default: run all models]
    -l TYPE, --LPType=TYPE
                        Loci Rate Priori Type
    -t RANGE, --timerange=RANGE
                        The range of values to select the time from
                        (Integers). Example: 1000:20000  [required]
    --obs_stats         When set, will generate a statistics output for the
                        observation data

  Posterior Run:
    Options to be applied during mode 'posterior'

    --uid_list=FILE     Speccifies a list of UIDs to filter on for Posterior
                        processing [required]
    --run_data=FILE     run data which contains the --uid_list UIDs.  It is
                        used for the Posterior processing [required]
```  

------------------------------------------------------------------------------------


##msReject module
We use this command with a reference table of 200,000 iterations per model to do the initial acceptance of 10,000 using msReject. 
```
msReject hyperstats_observations.txt reference_table.txt 0.0015151515151515151515 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 > Posterior
```
Instructions for installing the msReject module of msBayes can be found here:
[Instructions](https://docs.google.com/document/d/1enMQaogxOs0RppAmE8KcGU3nNjzotuiycAl6I1s0KYg/edit)

Then we use the following R-script and abc.R to do the final 1,000 acceptance and parameter estimation using local linear regression.

```
library(VGAM)
library(locfit)
library(abc)

#load the posterior file of 10,000 acceptances from msreject
Prior32 <- read.table("Posterior", sep="\t", header = F)
POSTVEC32<- Prior32[,c(17:32)]

#load summary statistics for observed data (produced by hBayeSSC.py)
MZ_OBS<- read.table("hyperstats_observations.txt", sep="\t", header = F)
OBS<- MZ_OBS[,c(17:32)] 

#########back transform function####################

#replace 32 with the number of species
backtrans_z <- function(val){
	    val <- val*32
	    if(val <= 0.5) {out <- 0}
	    else if((32-0.5) < val) {out <- 32}
	    else {out <- round(val,digit=0)}
	    out <- out/32
	    return(out)
	    }

##########start estimating z-value #################

Z1true1<- MZ_OBS[,4]
ZVALUE32 <- Prior32[,4]
RANVALUE1 <- ranprior[,4]

RANVALUE <- sample(RANVALUE1, 1000)

ZVALUE32_LL <- abc(OBS, ZVALUE32, POSTVEC32, tol=0.1,method="loclinear")

ZVALUE32_BT <- sapply(ZVALUE32_LL$adj.values,backtrans_z)
ZVALUE32_BT_mode <- loc1stats(ZVALUE32_BT, prob=0.95)[1]
ZVALUE32_BT_mode

summary(ZVALUE32_BT)
```
