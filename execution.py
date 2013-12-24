#!/usr/bin/python

import sys
import os
import copy
import random
import time
from math import sqrt
from optparse import OptionParser

"""
Assumptions: 
- Assume the population count is always 1
- Assume only 1 iteration in BayeSSC per run
"""


BAYESSC_PATH = ""
RETRIES = 10

class DataMap:
    """
    -- mapping for bayssc output to DataMap
    species.txt, obs[0]
    nsam.txt, obs[1]
    nsites.txt, obs[2]
    GROUP 0,
    Haptypes, haps
    PrivHaps,
    SegSites, seg
    PairDiffs, pair
    HapDiver, hapdiv
    NucltdDiv, nucdiv
    TajimasD, tajd
    F*, fusf
    MismatDist,
    COMBINED,
    Haptypes,
    PrivHaps,
    SegSites,
    PairDiffs,
    HapDiver,
    NucltdDiv,
    TajimasD,
    F*,
    MismatDist,
    MRCA,
    PRIORS,
    Deme Size 0, Ne
    Event Size 0, expan
    Mutation Rate 0, mu
    number.txt(randint) time
    """

    def __init__(self):
        self.species = ""
        self.nsam = 0.0
        self.nsites = 0.0
        self.haps = 0.0
        self.seg = 0.0
        self.pair = 0.0
        self.hapdiv = 0.0
        self.nucdiv = 0.0
        self.tajd = 0.0
        self.fusf = 0.0
        self.ne = 0.0
        self.expan = 0.0
        self.mu = 0.0
        self.time = 0.0

    def fillFromBayeSSC(self, obsData, time, statsData):
        """
        - obsData is a single line in the Observation file
        - statsdata is the data generated in the *_stat.csv from BayeSSC
        - time is a random int, between 2 user defined values.
        """
        self.species = obsData[0]
        self.nsam = obsData[1]
        self.nsites = obsData[2]
        self.haps = statsData['haptypes']
        self.seg = statsData['segsites']
        self.pair = statsData['pairdiffs']
        self.hapdiv = statsData['hapdiver']
        self.nucdiv = statsData['nucltddiv']
        self.tajd = statsData['tajimasd']
        self.fusf =  statsData['f*']
        self.ne = statsData['deme size']
        self.expan = statsData['event size']
        self.mu = statsData['mutation rate']
        self.time = time

    def __str__(self):
        return ",".join([self.species,
                         self.nsam,    self.nsites, 
                         self.haps,    self.seg, 
                         self.pair,    self.hapdiv, 
                         self.nucdiv,  self.tajd, 
                         self.fusf,    self.ne, 
                         self.expan,   self.mu, 
                         self.time])
            

class RunningStat(object):
    #http://www.johndcook.com/skewness_kurtosis.html
    # original class was written in c++.  This is a conversion to python
    def __init__(self):
        self.Clear()

    def Clear(self):
        self.n = 0
        self.M1 = 0.0
        self.M2 = 0.0
        self.M3 = 0.0
        self.M4 = 0.0

    def Push(self, x):
        x = float(x)
        delta = 0.0
        delta_n = 0.0
        delta_n2 = 0.0
        term1 = 0.0
        n1 = float(self.n)
        self.n += 1
        delta = x - self.M1
        delta_n = delta / float(self.n)
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        self.M1 += delta_n;
        self.M4 += term1 * delta_n2 * (self.n*self.n - 3*self.n + 3) + 6 * delta_n2 * self.M2 - 4 * delta_n * self.M3;
        self.M3 += term1 * delta_n * (self.n - 2) - 3 * delta_n * self.M2;
        self.M2 += term1;

    def NumDataValues(self):
        return self.n

    def Mean(self):
        return self.M1

    def Variance(self):
        return self.M2 / (self.n - 1.0)

    def StandardDeviation(self):
        return sqrt( self.Variance() )

    def Skewness(self):
        return sqrt(float(self.n)) * self.M3/ pow(self.M2, 1.5)

    def Kurtosis(self):
        return float(self.n) * self.M4 / (self.M2 * self.M2) - 3.0;
    
    def collectStats(self):
	# division by zero can occur when dealing with times (model 0 and model 1)
	# force the division by zero to result in all stats being 0. for our sanity!
        values = [0,0,0,0]
        try:
	   values[0] = self.Mean()
	except ZeroDivisionError:
            pass
        try:
	    values[1] = self.Variance()
	except ZeroDivisionError:
            pass
        try:
            values[2] =  self.Skewness()
	except ZeroDivisionError:
            pass
	try:
	    values[3] =  self.Kurtosis()
	except ZeroDivisionError:
            pass
        return values
   
def mergeRunningStats(a, b):
    """
    Could be part of the RunningStats class, but works just as well outside the class.
    Takes 2 RunningStat objects and combines their stored values into a 3rd RunningStat.
    """
    combined = RunningStat()
    combined.n += a.n + b.n
    delta = b.M1 - a.M1
    delta2 = delta * delta
    delta3 = delta * delta2
    delta4 = delta2 * delta2
    combined.M1 = (a.n * a.M1 + b.n * b.M1) / combined.n
    combined.M2 = a.M2 + b.M2 + delta2 * a.n * b.n / combined.n
    combined.M3 = a.M3 + b.M3 + delta3 * a.n * b.n * (a.n - b.n)/(combined.n * combined.n)
    combined.M3 += 3.0 * delta * (a.n * b.M2 - b.n * a.M2) / combined.n
    combined.M4 = a.M4 + b.M4 + delta4 * a.n * b.n * (a.n * a.n - a.n * b.n + b.n * b.n) / (combined.n * combined.n * combined.n)
    combined.M4 += 6.0 * delta2 * (a.n * a.n * b.M2 + b.n * b.n * a.M2) / (combined.n * combined.n) + 4.0 * delta * (a.n * b.M3 - b.n * a.M3) / combined.n
    return combined


class BadBayesOutput(Exception):
    """
    custom exception class to signify we had a problem with the BayeSSC output/execution
    """
    def __init__(self, val):
	self.val = val
    def __str__(self):
	return repr(self.val)
	

class ParFile(object):
    """
    A somewhat simplified parser of the par file format used by BayeSSC.
    For this particular script, we need to only pull out and modify specific fields in the par file
    In order to access these particular fields, we assume a fixed format for the par file
    and count the number of data lines in order to try and correctly parse the format.
    This results in a messy parser, but it should work.
    """
    def __init__(self, inName):
	infile = open(inName)
	values = []
	line = 0
	for l in infile:
	    l = l.strip()
	    if (not l) or l.startswith("//"):
		continue   
	    if line == 0:
		self.popcnt = l
		self.popsize = []
		if int(self.popcnt.split()[0]) == 0 :
		    line += 1
	    elif line == 1:
		c = int(self.popcnt.split()[0])
		self.popsize = [l] + [ infile.next().strip() for x in xrange(int(c)-1)]
	    elif line == 2:
		self.sSize = l
	    elif line == 3:
		self.growth = l
	    elif line == 4:
		self.mCnt = l
		self.migrations = []
		if int(self.mCnt) == 0:
		    line += 1
	    elif line == 5:
		self.migrations = [l] + [infile.next().strip() for x in xrange(int(l) - 1)]
	    elif line == 6: 
		self.eCnt = l
		self.events = []
		if int(l.split()[0]) == 0:
		    line += 1
	    elif line == 7:
		c = int(self.eCnt.split()[0])
		self.events = [self.bayesSplit(l)] + [ self.bayesSplit(infile.next().strip()) for x in xrange(int(c) - 1)]
	    elif line == 8:
		self.rate = l
	    elif line == 9:
		self.loci = l
	    elif line == 10:
		self.type = l
	    elif line == 11:
		self.gamma = l
	    line += 1
	infile.close()
	self.error =  line != 12


    def bayesSplit(self, line):
        """
        certain fields in the par file use {} to demark a single value.
        We utilize a similar idea found in the BayeSSC code that only looks
        for paired {} and considers whatever that range is, to be a single value
        """
	ele = []
	i = 0
	buf = ""

	while i < len(line):
	    c = line[i]
	    if c.isspace():
		if buf:
		    ele.append(buf)
		    buf = ""
		i += 1
		continue
	    if c == '{':
		pos = line[i:].find('}')
                if(pos == 1):
                    raise RuntimeError("Unbalanced { } in the par file")
		line = line[i:]
		ele.append( "".join( line[:pos+1].split() ) )
		line = line[pos+1:]
		i = 0
	    else:
		buf += c
		i += 1
	if(buf):
	    ele.append(buf)
	return " ".join([e.replace(" ", "") for e in ele])

    def setPopulation(self, distro, low, high):
	self.popsize = ["{%s:%s,%s}"%(distro, low, high) for r in self.popsize  ]
    
    def setTime(self, time):
	npop = []
	for r in self.events:
	    ele = r.split()
	    ele[0] = time
	    npop.append(" ".join(ele))
	self.events = npop        

    def setLociRate(self, distro, low, high):
	self.rate = "{%s:%s,%s}"%(distro, low, high)

    def setLoci(self, lociCnt):
	self.loci = lociCnt

    def setSampleSize(self, size):
	self.sSize = size

    def setTSTV(self, tstv):
	v = self.type.split()
	v[1] = tstv
	self.type = " ".join(v)
    
    def setgamma(self, gamma):
	v = self.gamma.split()
	v[0] = gamma
	self.gamma = " ".join(v)
    
    def matrixStr(self):
	if(self.migrations):
	    return "%s\n%s"%(self.mCnt, "\n".join(self.migrations))
	else:
	    return self.mCnt
    
    def eventStr(self):
	if(not self.events):
	    return ""
	return "\n".join(self.events)

    def __str__(self):       
	return """//Number of population samples - tmp\n%s\n//Population sizes\n%s\n//Sample sizes\n%s\nGrowth rates\n%s\n//Number of migration matrices : If 0 : No migration between demes\n%s\n//Historical event format:\n%s\n%s\n//mutation rate\n%s\n//Number of independent loci\n%s\n//Data type, num loci, rec.rate, mut rate, gamma, shape\n%s\n//\n%s\n"""%(
                self.popcnt, "\n".join(self.popsize), self.sSize, self.growth, self.matrixStr(), self.eCnt, self.eventStr(), self.rate, self.loci, self.type, self.gamma)


class ObservationSplitter:
    """
    a simple class that depending on which splitType is selected, can be used to split the obvervations up
    """
    def __init__(self, splitType = 'uniform'):
        # select the splitter to use, based on user provided split type
        splitters ={'identity':self.__splitIndentity, 'uniform':self.__splitUniform}
        self.split = splitters.get(splitType, None)
        if(not self.split):
            raise KeyError("Unknown splitter. Valid choices: %s"%(splitters.keys()))
        
    def __splitIndentity(self, observations, con_species):
        conSpecs = observations[:con_species]
        randSpecs = observations[con_species:]
        return conSpecs, randSpecs, observations
    
    def __splitUniform(self, observations, con_species, shuffleAction = random.random):
        """
        Take the array of observations, and shuffle them using the random() function in python.
        random should be uniform to start, so as a result, shuffle should be uniform as well.
        """
        random.shuffle(observations, random = shuffleAction)	
        return self.__splitIndentity(observations, con_species)


def is_exe(fpath):
    # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def which(program):
    # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    Validate taht the user provided path, does infact exist for BayeSSC
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
    

def generatedUniformTime(timerange):
    """
    generated a time value between the high and low value provided
    """
    return random.randint(min(timerange), max(timerange))


def getStatsPath(parPath):
    """
    convert the par file path into the bayessc stats file path
    """
    return os.path.splitext(parPath)[0] + "_stat.csv"
   

def parseBayeSSCOut(filePath):
    """
    takes the output fro BayeSSC and parses it so that we 
    """
    # with our assumption of 1 pop, combined and group 0 will be the same, so filter out the dup and store in a dict with the 0 removed
    statsf = open(filePath)
    try:
	hdr = [ c.replace(" 0", "").strip().lower() for c in statsf.next().strip().split(",")]
    except:
	raise  BadBayesOutput("Header Not Found")
    datadict = {}	
    d = ""
    try:
	d = statsf.next().strip().split(",")
    except:
	raise BadBayesOutput("Data Not Found")		
    if len(d) != len(hdr):
	raise BadBayesOutput("Data row is only partial")
    used = {}
    for h, v in zip(hdr, d):
	if h not in used and v:
	    datadict[h] = v
	    used[h] = None	
    statsf.close()
    return datadict	


def parseObs(obs):
    """
    The observation file is a tab delimited file. This can easily be
    split and stored in memory.
    """
    obsf = open(obs)
    obsl = [l.strip().split() for l in obsf]
    obsf.close()
    return obsl


def runBayeSSC(obs, time, par, workdir = ".", parname = "tmp.par"):
    """
    Execute BayeSSC and then parse the data generated by the run.
    """
    global BAYESSC_PATH
    fpath = os.path.join(workdir, parname)
    o = open(fpath, "w")
    o.write(par.__str__())
    o.close()
    #print "%s -f %s 1 &> /dev/null"%(BAYESSC_PATH, fpath)   
    os.system("%s -f %s 1 &>/dev/null"%(BAYESSC_PATH, fpath))
    data = parseBayeSSCOut(getStatsPath(fpath))
    dmap = DataMap()
    dmap.fillFromBayeSSC(obs, time, data)
    return dmap


def prepareBayeSSC(obs, parData, time, LPType, PopType):
    """
    populate the par object with the correct values.  Also modify the timestamp base on data from obs file
    """
    chngtime = str( int( time / float(obs[5])) )
    par = copy.copy(parData)
    par.setPopulation(PopType,obs[8],obs[9])
    par.setTime(chngtime)
    par.setLociRate(LPType, obs[6], obs[7])
    par.setgamma(obs[4])
    par.setSampleSize(obs[1])
    par.setLoci(obs[2])
    par.setTSTV(obs[3])
    return chngtime, par


def exceuteBateSSCWIthRetry(obs, chngtime, par, workdir):
    global RETRIES
    attempts = 0
    run = True
    while run:
        try:
            row = runBayeSSC(obs, chngtime, par, workdir)
            run = False
            return row
        except BadBayesOutput:
            attempts += 1
            if attempts >= RETRIES:
                raise BadBayesOutput("Attempted to run BayeSSC %s times, each run resulted in an output error.", RETRIES)
            print >> sys.stderr, "Error running bayeSSC.  Trying again"    
    return None


def generateCONSpecs(origParName, parData, observations, counts, timerange, workdir = ".", LPType = "U", PopType = "U"):
    """
    Iterate all observations in the Congruent group and execute BayeSSC using that particular observation data
    """
    #congruent species
    rows = []
    time = float(generatedUniformTime(timerange) )   
    for obs in observations:
	counts[obs[0]] += 1
        chngtime, par = prepareBayeSSC(obs, parData, time, LPType, PopType)
        row = exceuteBateSSCWIthRetry(obs, chngtime, par, workdir)
        if row:
            rows.append(row)
    if len(rows) != len(observations):
        raise BadBayesOutput("Did not generate an output for each observation")            
    return rows, counts


def generateRANDSpecs(origParName, parData, observations, timerange, workdir = ".", LPType = "U", PopType = "U"):
    """
    Iterate all observations in the Random group and execute BayeSSC using that particular observation data.
    A new timestamp is generated for each observations. 
    """
    rows = []
    for obs in observations:
        # for each observation, get a new time
        # moved outside the while True, since it makes code consolidation easier
        time = float(generatedUniformTime(timerange) )
        chngtime, par = prepareBayeSSC(obs, parData, time, LPType, PopType)
        row = exceuteBateSSCWIthRetry(obs, chngtime, par, workdir)
        if row:
            rows.append(row)
    if len(rows) != len(observations):
        raise BadBayesOutput("Did not generate an output for each observation")
    return rows


def computeStats(congruentCnt, total, conspecData, randomData):
    """
    Using the collect data from multipl repeats, compute some statistics on particular data columns
    """
    #time = RunningStat()	
    expan = RunningStat()
    mu = RunningStat()
    Ne = RunningStat()	
    contime = RunningStat()	
    rndtime = RunningStat()	
    haps = RunningStat()	
    hapdiv = RunningStat()	
    nucdiv = RunningStat()	
    pair  = RunningStat()	
    tajd = RunningStat()	
    fusf = RunningStat()	

    for row in conspecData:
	expan.Push(row.expan)
	mu.Push(row.mu)
	Ne.Push(row.ne)
	#time.Push(row.time)
	contime.Push(row.time)
	haps.Push(row.haps)
	hapdiv.Push(row.hapdiv)
	nucdiv.Push(row.nucdiv)
	pair.Push(row.pair)
	tajd.Push(row.tajd)
	fusf.Push(row.fusf)

    for row in randomData:
	expan.Push(row.expan)
	mu.Push(row.mu)
	Ne.Push(row.ne)
	#time.Push(row.time)
	rndtime.Push(row.time)
	haps.Push(row.haps)
	hapdiv.Push(row.hapdiv)
	nucdiv.Push(row.nucdiv)
	pair.Push(row.pair)
	tajd.Push(row.tajd)
	fusf.Push(row.fusf)

    overalltime = mergeRunningStats(rndtime, contime)
    stats = [congruentCnt, total, float(congruentCnt) / float(total) , overalltime.Variance() / overalltime.Mean()] 

    stats.extend(contime.collectStats())
    stats.extend(rndtime.collectStats())
    stats.extend(overalltime.collectStats())
    #stats.extend(time.collectStats())
    stats.extend(Ne.collectStats())
    stats.extend(expan.collectStats())
    stats.extend(mu.collectStats())
    stats.extend(haps.collectStats())
    stats.extend(hapdiv.collectStats())
    stats.extend(nucdiv.collectStats())
    stats.extend(tajd.collectStats())
    stats.extend(fusf.collectStats())
    stats.extend(pair.collectStats())
    return map(str, stats)
    

def performSingleModel(splitter, timerange, uid, outdir, parName, hyperstats, observations, repeats, LPType, par, con_species, obsCnt, counts):
    """
    a model describes how many observations make up the congruent group.  for instance, model0, means we have no congruent observations.
    This method is meant to contain all actions required to execute this script on a single model.  It will take that model, and
    repeat the experiment multiple times, each time splitting the observations again and again.
    """
    print >> sys.stderr , ".",
    o = open(os.path.join(outdir, "con%s_total%s_-_%s_iterations.csv"%(con_species, obsCnt, repeats)), "w")
    for trial in xrange(int(repeats)):
	conSpecs, randSpecs, observations = splitter.split(observations, con_species)
	outstr = []
	conspecData= []
	randomData= []
	index = "%s_%s_%s_%s_%s"%(uid, obsCnt, con_species, trial, "_".join([str(random.random()), str(time.time())]).replace(".","_"))
	if conSpecs:
	    rows, counts = generateCONSpecs(parName, par, conSpecs, counts, timerange, workdir = outdir, LPType = LPType)
	    conspecData.extend(rows)
	    outstr.append( ",".join( ( str(r) for r in rows) ) )
	if randSpecs:
	    rows = generateRANDSpecs(parName, par, randSpecs, timerange, workdir = outdir, LPType = LPType)
	    randomData.extend(rows)
	    outstr.append( ",".join( ( str(r) for r in rows) ) )
	hyperstats.write( "%s,%s\n"%(index, ",".join( computeStats(len(conSpecs), obsCnt, conspecData, randomData) )) )
	o.write("%s,%s\n"%(index, ",".join(outstr)))
    o.close()
    o = open(os.path.join(outdir, "con%s_total%s_-_%s_iterations_congruent_appearance.csv"%(con_species, obsCnt, repeats)), "w")
    o.write("\n".join( [",".join(map(str, itm)) for itm in counts.iteritems() ] ) )
    o.close()
	

	
def commandlineArgs():
    """
    Command line arguments for this script.  It also provides some validation on the arguments passed in, to make sure we have some type of chance to actually succeede.
    """
    global BAYESSC_PATH
    parser = OptionParser("%prog [options]")
    parser.add_option("-p", "--par", dest = "par", help = "par file template [required]", action = "store", type = "string", metavar = "FILE")
    parser.add_option("-i", "--obs", dest = "obs", help = "Observation file [required]", action = "store", type = "string", metavar = "FILE")
    parser.add_option("-r", "--repeat", dest = "repeats", help = "Number of times to try a given congruent group size [required]", action = "store", type = "int", metavar = "NUM")
    parser.add_option("-m", "--model", dest = "model", help = "Run a single model (0 to total entries in observation file) [default: run all models] ", action = "store", type = "int", metavar = "MODEL", default = None)
    parser.add_option("-u", "--uid", dest = "uid", help = "Unique ID to prefix generated indices [required]", action = "store", type = "string", metavar = "UID")
    parser.add_option("-l", "--LRType", dest = "LRType", help = "Loci Rate Priori Type", action = "store", type = "choice", choices = ["U"], default = "U", metavar = "TYPE")
    parser.add_option("-o", "--outdir", dest = "outdir", help = "Directory to generate output in (will create missing folders) [default: %default]", action = "store", type = "string", metavar = "PATH", default = os.getcwd())
    parser.add_option("-t", "--timerange", dest= "trange", help = "The range of values to select the time from (Integers). Example: 1000:20000  [required]", action = "store", type = "string", metavar ="RANGE")
    parser.add_option("-b", "--bayepath", dest = "bayesPath", help = "Path to BayeSSC application [default: Located on user PATH]", action = "store", type = "string", metavar = "PATH", default = "BayeSSC")
    
    (options, args) = parser.parse_args()    
    
    if not options.par:
	parser.print_help()
	parser.error("par file is required")
    if not options.obs:
	parser.print_help()
	parser.error("observation file is required")
    if not options.repeats:
	parser.print_help()
	parser.error("Number of repeats is required")
    if not options.uid:
	parser.print_help()
	parser.error("A Unique ID is required")
    if os.path.exists(options.outdir) and  not os.path.isdir(options.outdir):
	parser.print_help()
	parser.error("Output path exists, but is not a directory")
    if not os.path.exists(options.outdir):
	try:
	    os.makedirs(options.outdir)
	except OSError:
	    parser.print_help()
	    parser.error("Output path cannot be created")
    options.uid = options.uid.replace(",","_").replace(" ","")
    BAYESSC_PATH = which(options.bayesPath)
    if not BAYESSC_PATH:
	parser.print_help()
	parser.error("BayeSSC application not found at supplied path: '%s'" %(options.bayesPath))
    if not options.trange:
	parser.print_help()
	parser.error("Time range is required")
    if options.trange.find(".") != -1:
	parser.print_help()
	parser.error("Time range must consist of only integers")
    options.trange = options.trange.split(":")
    if len(options.trange) != 2:
	parser.print_help()
	parser.error("Time range  must be provided in the following format:  <lowerbounds>:<upperbounds> Example: 1000:20000")
    try:
	options.trange = map(int, options.trange)
    except:
	parser.print_help()
	parser.error("Time range does not consist of valid integers")
    return options


def main():
    """
    Main loop of the appliocation
    drives how the program executes (only 1 model, or multiple models).
    """
    options = commandlineArgs()
    par = ParFile(options.par)
    observations = parseObs(options.obs)
    obsCnt = len(observations)
    fname = "hyperstats_-_%s.txt"
    if options.model == None:
	fname = fname%("models_0-%s_-_%s_iterations"%(obsCnt, options.repeats))
    else:
	fname = fname%("model%s_of_%s"%(options.model, obsCnt))
    hyperstats = open(os.path.join(options.outdir, fname), "w")
    counts = dict([ (obs[0], 0) for obs in observations])
    
    splitter = ObservationSplitter("uniform")
    if options.model == None:
	for con_species in xrange(obsCnt + 1):
	    for k in counts:
		counts[k] = 0
	    performSingleModel(splitter, options.trange, options.uid, options.outdir, options.par, hyperstats, observations, options.repeats, options.LRType, par, con_species, obsCnt, counts)
    else:
	performSingleModel(splitter, options.trange, options.uid, options.outdir, options.par, hyperstats, observations, options.repeats, options.LRType, par, options.model, obsCnt, counts)


if __name__ == "__main__":
    main()
