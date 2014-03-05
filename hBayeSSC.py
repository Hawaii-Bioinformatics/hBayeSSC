#!/usr/bin/python

from itertools import izip, chain
import sys
import os
import copy
import random
import time
import string
from math import sqrt
try:
    from math import isnan
except:
    def isnan(val):
        return not (float('-inf') <= val <= float('inf'))
try:
    from math import isinf
except:
    def isinf(val):
        return val == float('-inf') or val == float('inf')

from optparse import OptionParser


"""
--note that we are not changing ==anything== in the runs.  we just want to run a specific set of params.--
epic888, 3/3/2014 10:23:04 PM:
It looks like you write out the time, mutation rate, expansion..

David S, 3/3/2014 10:41:51 PM:
1. modify the hbayessc so that it can take a UID list  (or  posterior file) and a data_run file

David S, 10:42:10 PM:
2. generate new par files using the exact parameters define din the data_run file per given UID

David S, 10:42:36 PM:
3. run bayessc for 1 iteration and make a new hyperstats (allow N iterations)




# script to count appearances from hyperstats or UID list
# script to pull specific models by UID

Author: David Schanzenbach
Original script by: Yvonne Chan

Original Version: 20131224
- Should work with python 2.4 and any newer python 2.x release

Assumptions: 
- Assume the population count is always 1
- Assume only 1 iteration in BayeSSC per run
"""

def skipNanInf(value):
    """
    The default filtering action used by RunningStat when using Push()
    """
    try:
	value = float(value)
	return isnan(value) or isinf(value)
    except ValueError:
	return True


class BadBayesOutput(Exception):
    """
    custom exception class to signify we had a problem with the BayeSSC output/execution
    """
    def __init__(self, val):
	self.val = val
    def __str__(self):
	return repr(self.val)


class CommonData(object):
 
    def __init__(self):
        self.label = ""
        self.nsam = 0.0
        self.nsites = 0.0
        self.haps = 0.0
        self.seg = 0.0
        self.pair = 0.0
        self.hapdiv = 0.0
        self.nucdiv = 0.0
        self.tajd = 0.0
        self.fusf = 0.0

	
    def fill(self, label, nsam, nsites, data, nucdiv):
	self.label = label
	self.nsam = nsam
	self.nsites = nsites
	self.haps = data.get('haptypes', float("NaN") )
	self.seg = data.get('segsites', float("NaN") )
	self.pair = data.get('pairdiffs', float("NaN") )
	self.nucdiv = nucdiv #data.get('nucdiv', data.get('nucltddiv',float("NaN")) )
	self.hapdiv = data.get('hapdiver', float("NaN") )
	self.tajd = data.get('tajimasd', float("NaN") )
	self.fusf = data.get('f*', float("NaN") )

    def addStats(self, statsdict):
	statsdict['haps'].Push(self.haps)
	statsdict['hapdiv'].Push(self.hapdiv)
	statsdict['nucdiv'].Push(self.nucdiv)
	statsdict['pair'].Push(self.pair)
	statsdict['tajd'].Push(self.tajd)
	statsdict['fusf'].Push(self.fusf)
	return statsdict

 
class ObservationData(CommonData):
    """ Represents a row in the Observeration file """
    # default column order
    columns = ['species','nsam','nsites','tstv','gamma','gen','locuslow','locushigh','nelow','nehigh','segsites','nucdiv','haptypes','hapdiver','pairdiffs','tajimasd','f*','exphet']
    def __init__(self, data = None):
	super(ObservationData, self).__init__()
	self.gamma = 0.0
	self.gen = 0.0
	self.locuslow = 0.0
	self.locushigh = 0.0
	self.neLow = 0.0
	self.neHigh = 0.0
	self.exphet = 0.0
	self.tstv = 0.0
        if data != None:
            self.fill(data)

    def fill(self, data):
        """
        - obsData is a single line in the Observation file
        - statsdata is the data generated in the *_stat.csv from BayeSSC
        - time is a random int, between 2 user defined values.
        """
	super(ObservationData, self).fill(data['species'], data['nsam'], data['nsites'], data, data.get('nucdiv', float("NaN") )) 	           
        self.tstv = data['tstv']
        self.gamma = data['gamma']
        self.gen = data['gen']
        self.locuslow = data['locuslow']
        self.locushigh = data['locushigh']
        self.neLow = data['nelow']
        self.neHigh = data['nehigh']
        self.exphet = data.get('exphet', float("NaN") )


    def addStats(self, statsdict):
	return super(ObservationData, self).addStats(statsdict)

    def getPopRange(self):
        return (float(self.neLow), float(self.neHigh))
        
    def getMutationRange(self):
        return (float(self.locuslow), float(self.locushigh))

    def __str__(self):
        return "\t".join(map(str, [self.label, self.nsam, self.nsites,
                                  self.tstv, self.gamma, self.gen, self.locushigh, 
                                  self.locuslow, self.neLow, self.neHigh, self.seg,
                                  self.nucdiv, self.haps, self.hapdiv, self.pair, 
                                  self.tajd, self.fusf, self.exphet]) )


class BayeSSCData(CommonData):
    """Represents a row in the BayeSSC stats file """
    def __init__(self, obs = None, time = None, data = None):
	super(BayeSSCData, self).__init__()
	self.ne = 0.0
	self.expan = 0.0
	self.mu = 0.0
	self.time = 0.0
        if obs != None and time != None and data != None:
	    self.fill(obs, time, data)
	
    def fill(self, obsData, time, statsData):
        """
        - obsData is a single line in the Observation file
        - statsdata is the data generated in the *_stat.csv from BayeSSC
        - time is a random int, between 2 user defined values.
        """
	super(BayeSSCData, self).fill(obsData.label, obsData.nsam, obsData.nsites, statsData, statsData.get('nucltddiv', float("NaN") ))
	self.ne = statsData['deme size']
	self.expan = statsData['event size']
	self.mu = statsData['mutation rate']
	self.time = time

    def addStats(self, statsdict):
	statsdict = super(BayeSSCData, self).addStats(statsdict)

	statsdict['mu'].Push(self.mu)
	statsdict['ne'].Push(self.ne)
	statsdict['expan'].Push(self.expan)
        #print self.ne
	return statsdict
    
    def header(self):
	return ['species', 'nsam','nsites', 'haptype', 'segsites', 'pairdiffs', 'hapdiv', 'nucdiv', 'tajimasd', 'fusf','ne', 'expan', 'mu', 'time']

    def __str__(self):
        return "\t".join(map(str, [self.label,
                         self.nsam,    self.nsites, 
                         self.haps,    self.seg, 
                         self.pair,    self.hapdiv, 
                         self.nucdiv,  self.tajd, 
                         self.fusf,    self.ne, 
                         self.expan,   self.mu, 
                         self.time]))
            

class RunningStat(object):
    #http://www.johndcook.com/skewness_kurtosis.html
    # original class was written in c++.  This is a conversion to python
    # follows the online_kurtosis function describe at wikipedia:
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
    def __init__(self):
        self.Clear()

    def Clear(self):
        self.n = 0
        self.M1 = 0.0
        self.M2 = 0.0
        self.M3 = 0.0
        self.M4 = 0.0

    def Push(self, x, filterFunction = skipNanInf):
        x = float(x)
        if filterFunction and filterFunction(x):
            return
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
        self.M4 += term1 * delta_n2 * (self.n * self.n - 3 * self.n + 3) + 6 * delta_n2 * self.M2 - 4 * delta_n * self.M3;
        self.M3 += term1 * delta_n * (self.n - 2) - 3 * delta_n * self.M2;
        self.M2 += term1;

    def NumDataValues(self):
        return self.n

    def Mean(self):
        if self.NumDataValues() == 0:
            return float('NaN')
        return self.M1

    def Variance(self):
        return self.M2 / (self.n - 1.0)

    def StandardDeviation(self):
        return sqrt( self.Variance() )

    def Skewness(self):
        return sqrt(float(self.n)) * self.M3/ pow(self.M2, 1.5)

    def Kurtosis(self):
        return float(self.n) * self.M4 / (self.M2 * self.M2) - 3.0;

    def Dispersion(self):
        try:
            return self.Variance() / self.Mean()
        except ZeroDivisionError:
            return float('Nan')
        
    def collectMeanAndVariance(self, values = [float('NaN'), float('NaN')] ):
	if self.n == 0:
	    return values
        try:
            values[0] = "%.15f"%(self.Mean())
        except ZeroDivisionError:
            pass
        try:
            values[1] = "%.15f"%(self.Variance())
        except ZeroDivisionError:
            pass
        return values

    def collectStats(self):
        # division by zero can occur when dealing with times (model 0 and model 1)
        # force the division by zero to result in all stats being NaN. for our sanity!
        values = [float('NaN'), float('NaN'), float('NaN'), float('NaN')]
	if self.NumDataValues() == 0:
	    return values	
	values = self.collectMeanAndVariance(values)
        try:
            values[2] = "%.15f"%(self.Skewness())
        except ZeroDivisionError:
            pass
        try:
            values[3] = "%.15f"%(self.Kurtosis())
        except ZeroDivisionError:
            pass
        return values


class ParFile(object):
    """
    A somewhat simplified parser of the par file format used by BayeSSC.
    For this particular script, we need to only pull out and modify specific fields in the par file
    In order to access these particular fields, we assume a fixed format for the par file
    and count the number of data lines in order to try and correctly parse the format.
    This results in a messy parser, but it should work.
    """
    def __init__(self, inName):
	infile = open(inName, "rU")
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
		self.events = [self.splitPrior(l)] + [ self.splitPrior(infile.next().strip()) for x in xrange(int(c) - 1)]
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

    def splitPrior(self, line):
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

    def setPopulation(self, distro, poprng):
	self.popsize = ["{%s:%.15f,%.15f}"%(distro, poprng[0], poprng[1]) for r in self.popsize]
    
    def setTime(self, time):
	npop = []
	for r in self.events:
	    ele = r.split()
	    ele[0] = time
	    npop.append(" ".join(ele))
	self.events = npop        

    def setLociRate(self, distro, locrng):
	self.rate = "{%s:%.15f,%.15f}"%(distro, locrng[0], locrng[1])

    def setLoci(self, lociCnt):
	self.loci = lociCnt

    def setSampleSize(self, size):
	self.sSize = size

    def setTSTV(self, tstv):
	v = self.type.split()
	v[1] = "%.15f"%(float(tstv))
	self.type = " ".join(v)
    
    def setgamma(self, gamma):
	v = self.gamma.split()
	v[0] = "%.15f"%(float(gamma))
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
	return """//Number of population samples - tmp\n%s\n//Population sizes\n%s\n//Sample sizes\n%s\n//Growth rates\n%s\n//Number of migration matrices : If 0 : No migration between demes\n%s\n//Historical event format:\n%s\n%s\n//mutation rate\n%s\n//Number of independent loci\n%s\n//Data type, num loci, rec.rate, mut rate, gamma, shape\n%s\n//\n%s\n"""%(
                self.popcnt, "\n".join(self.popsize), self.sSize, self.growth, self.matrixStr(), self.eCnt, self.eventStr(), self.rate, self.loci, self.type, self.gamma)


class ObservationSplitter(object):
    """ a simple class that depending on which splitType is selected, can be used to split the obvervations up """
    def __init__(self, splitType = 'uniform'):
        # select the splitter to use, based on user provided split type
        splitters ={'identity':self.__splitIndentity, 'uniform':self.__splitUniform}
        self.split = splitters.get(splitType, None)
        if(not self.split):
            raise KeyError("Unknown splitter. Valid choices: %s"%(splitters.keys()))
        
    def __splitIndentity(self, observations, con_species):
        conSpecs = observations[:con_species]
        randSpecs = observations[con_species:]
        return conSpecs, randSpecs #, observations
    
    def __splitUniform(self, observations, con_species, shuffleAction = random.random):
        """
        Take the array of observations, and shuffle them using the random() function in python.
        random should be uniform to start, so as a result, shuffle should be uniform as well.
        """
        random.shuffle(observations, random = shuffleAction)	
        return self.__splitIndentity(observations, con_species)


class TimeGenerator(object):
    """ a simple class that depending on which genType is selected, can be used to generate a time value.  """
    def __init__(self, genType = 'uniform'):
        timeGens ={'uniform':self.__generatedUniformTime}
        self.generate = timeGens.get(genType, None)
        if(not self.generate):
            raise KeyError("Unknown Time generator. Valid choices: %s"%(timeGens.keys()))

    def __generatedUniformTime(self, timerange):
        """ generated a time value between the high and low value provided """
        return random.randint(min(timerange), max(timerange))

	   
def mergeRunningStats(a, b):
    """ Takes 2 RunningStat objects and combines their stored values into a 3rd RunningStat. """
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


def statsHeader():
    return ['congruent_group_size', 'total_observations', 'model_pct',
     'congruent_time_Mean', 'congruent_time_Dispersion',
     'random_time_Mean', 'random_time_Dispersion',
     'overall_time_Mean', 'overall_time_Dispersion',
     'ne_Mean', 'ne_Variance',
     'expan_Mean', 'expan_Variance',
     'mu_Mean', 'mu_Variance',
     'haptypes_Mean', 'haptypes_Variance', 'haptypes_Skewness', 'haptypes_Kurtosis',
     'hapdiv_Mean', 'hapdiv_Variance', 'hapdiv_Skewness', 'hapdiv_Kurtosis',
     'nucdiv_Mean', 'nucdiv_Variance', 'nucdiv_Skewness', 'nucdiv_Kurtosis',
     'tajimasd_Mean', 'tajimasd_Variance', 'tajimasd_Skewness', 'tajimasd_Kurtosis',
     'fusf_Mean', 'fusf_Variance', 'fusf_Skewness', 'fusf_Kurtosis',	
     'pairDiffs_Mean', 'pairDiffs_Variance', 'pairDiffs_Skewness', 'pairDiffs_Kurtosis']

def computeStats(congruentCnt, total, conspecData = None, randomData = None, obsData = None):
    """ Using the collect data from multipl repeats, compute some statistics on particular data columns """
    statsdict = dict(expan = RunningStat(), mu = RunningStat(), ne = RunningStat(), contime = RunningStat(), 
                 rndtime = RunningStat(), haps = RunningStat(), hapdiv = RunningStat(), nucdiv = RunningStat(), 
                 pair = RunningStat(), tajd = RunningStat(), fusf = RunningStat(), overalltime = RunningStat())
    if conspecData or randomData:
        for row in conspecData:
	    statsdict = row.addStats(statsdict)
            statsdict['contime'].Push(row.time)
        for row in randomData:
	    statsdict = row.addStats(statsdict)
            statsdict['rndtime'].Push(row.time)
        statsdict['overalltime'] = mergeRunningStats( statsdict['rndtime'], statsdict['contime'])
        stats = [congruentCnt, total, "%.15f"%(float(congruentCnt) / float(total)) ]             
    elif obsData:
        # compute the stats for the observation file.  Computer only those that are identical to the hyperstats output
        for row in obsData:
	    statsdict = row.addStats(statsdict)
        stats = [float('NaN'), float('NaN'), float('NaN')]
    else:
        raise BadBayesOutput("No observation data, congruent data or random data found to compute stats on")

    stats.extend(chain( *[ [statsdict[k].Mean(), statsdict[k].Dispersion() ] for k in ['contime', 'rndtime', 'overalltime']] ))
    tmp = []
    for k in  ['ne', 'expan', 'mu']:
        v = [float('nan'), float('nan')]
        v  = statsdict[k].collectMeanAndVariance(v)
        tmp.append(v)
    stats.extend(chain( *tmp))
    stats.extend(chain( *[statsdict[k].collectStats() for k in ['haps', 'hapdiv', 'nucdiv', 'tajd', 'fusf', 'pair']] ))
    return map(str, stats)
    

class BayeSSC(object):
    """ A class that represents the execution and parsing of the BayeSSC application """

    def __init__(self, execpath, retries = 10):
        self.execpath = execpath
        self.retires = retries
        
    def __getStatsPath(self, parPath):
        """ convert the par file path into the bayessc stats file path. """
        # TODO: DLS - This may break under certain conditions.  need to verify
	#	parfile = os.path.split(parPath)[-1]
	#	parfile = os.path.splitext(parfile)[0] + "_stat.csv"
	#    return os.path.join(workdir, parfile)
        return os.path.splitext(parPath)[0] + "_stat.csv"

    def __parseBayeSSCOut(self, filePath):
        """ takes the output from BayeSSC and places it into a dictionary """
        # with our assumption of 1 pop, combined and group 0 will be the same, so filter out the dup and store in a dict with the 0 removed
        statsf = open(filePath, "rU")
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
        for h, v in izip(hdr, d):
            if h not in used and v:
                datadict[h] = v
                used[h] = None	
        statsf.close()
        return datadict	

    def runBayeSSC(self, obs, ctime, par, outdir = ".", parname = "tmp.par"):
        """ Execute BayeSSC and then parse the data generated by the run. """
        fpath = os.path.join(outdir, parname)
        o = open(fpath, "w")
        o.write(par.__str__())
        o.close()
        os.system("%s -f %s 1 2>/dev/null >/dev/null"%(self.execpath, fpath))
        data = self.__parseBayeSSCOut(self.__getStatsPath(fpath))
        return BayeSSCData(obs, ctime, data)

    def exceuteBateSSCWithRetry(self, obs, chngtime, par, outdir):
        for x in xrange(self.retires):
            try:
                bayeData = self.runBayeSSC(obs, chngtime, par, outdir)
                return bayeData
            except BadBayesOutput:
                print >> sys.stderr, "Error running bayeSSC.  Trying again"    
        raise BadBayesOutput("Attempted to run BayeSSC %s times, each run resulted in an output error." %(self.retires))


class Model(object):
    FIELD_DELIM = "\t"

    def __init__(self, options, par, observations, totalObservations, splitter, timegen, bayessc):
        self.splitter = splitter
        self.par = par
        self.observations = observations
        self.obsCnt = totalObservations
        self.options  = options	
        self.indx =  "%s_%s_%%s_%%%%s_%%s"%(self.options.uid, self.obsCnt)
        self.bayessc = bayessc
        self.timeGenerator = timegen
		
    def execute(self, modelNumber, hyperstatsOut = None, runDatOut = None):
        """
        a model describes how many observations make up the congruent group.  for instance, model0, means we have no congruent observations.
        This method is meant to contain all actions required to execute this script on a single model.  It will take that model, and
        repeat the experiment multiple times, each time splitting the observations again and again.
        """
        print >> sys.stderr , ".",
        indx_raw = self.indx%(modelNumber, "_".join([str(random.random()), str(time.time())]).replace(".","_"))
        for trial in xrange(int(self.options.repeats)):
            conSpecs, randSpecs = self.splitter.split(self.observations, modelNumber)
            outstr = []
            conspecData= []
            randomData= []
            indx = indx_raw%(trial)
            if conSpecs:
                rows = self.__generateCONSpecs(self.options.par, self.par, conSpecs, self.options.trange, self.options.outdir, LPType = self.options.LPType)
                conspecData.extend(rows)
                outstr.append( Model.FIELD_DELIM.join( map(str, rows) ) )
            if randSpecs:
                rows = self.__generateRANDSpecs(self.options.par, self.par, randSpecs, self.options.trange, self.options.outdir, LPType = self.options.LPType)
                randomData.extend(rows)
                outstr.append( Model.FIELD_DELIM.join( map(str, rows) ) )

	    print >> hyperstatsOut, Model.FIELD_DELIM.join( [indx] + computeStats(len(conSpecs), self.obsCnt, conspecData, randomData) )
            if runDatOut:
                print >> runDatOut, Model.FIELD_DELIM.join( [indx] + outstr)

    def __commonExec(self, obs, parData, time, LPType, PopType, outdir, rows, modifyTime = True):
        chngtime, par = prepareNewParFile(obs, parData, time, LPType, PopType, modifyTime)
        row = self.bayessc.exceuteBateSSCWithRetry(obs, chngtime, par, outdir)
        if row:
            rows.append(row)
        return rows

    def __generateCONSpecs(self, origParName, parData, observations, timerange, outdir, LPType = "U", PopType = "U", timeGenerator =  TimeGenerator("uniform")):
        """ Iterate all observations in the Congruent group and execute BayeSSC using that particular observation data """
        rows = []
        time = float(self.timeGenerator.generate(timerange) )   
        for obs in observations:
            rows = self.__commonExec(obs, parData, time, LPType, PopType, outdir, rows)
        if len(rows) != len(observations):
            raise BadBayesOutput("Did not generate an output for each observation")            
        return rows

    def __generateRANDSpecs(self, origParName, parData, observations, timerange, outdir, LPType = "U", PopType = "U", timeGenerator =  TimeGenerator("uniform")):
        """
        Iterate all observations in the Random group and execute BayeSSC using that particular observation data.
        A new timestamp is generated for each observations. 
        """
        rows = []
        for obs in observations:
            time = float(timeGenerator.generate(timerange) )
            rows = self.__commonExec(obs, parData, time, LPType, PopType, outdir, rows)
        if len(rows) != len(observations):
            raise BadBayesOutput("Did not generate an output for each observation")
        return rows


class PostModel(Model):
    
    def __init__(self,options, par, conSpecs, conSpecsTimes, randSpecs, randSpecsTimes, bayessc):	
	super(PostModel, self).__init__(options, par, None, len(conSpecs) + len(randSpecs), None, None, bayessc)
	self.conSpecs = conSpecs
	self.conSpecsTimes = conSpecsTimes
	self.randSpecs = randSpecs
	self.randSpecsTimes= randSpecsTimes
	
   def __generatePOST(self, origParName, parData, observations, times, outdir, LPType = "U", PopType = "U"):
	""" Iterate all observations in the Congruent group and execute BayeSSC using that particular observation data """
	rows = []
	for obs, time in izip(observations, times):
	    rows = self.__commonExec(obs, parData, time, LPType, PopType, outdir, rows, False)
	if len(rows) != len(observations):
	    raise BadBayesOutput("Did not generate an output for each observation")            
	return rows	
	
    def execute(self, modelNumber, hyperstatsOut = None, runDatOut = None):
        """
        a model describes how many observations make up the congruent group.  for instance, model0, means we have no congruent observations.
        This method is meant to contain all actions required to execute this script on a single model.  It will take that model, and
        repeat the experiment multiple times, each time splitting the observations again and again.
        """
        print >> sys.stderr , ".",
        indx_raw = self.indx%(modelNumber, "_".join([str(random.random()), str(time.time())]).replace(".","_"))
        for trial in xrange(int(self.options.repeats)):
            outstr = []
            conspecData= []
            randomData= []
            indx = indx_raw%(trial)
            if conSpecs:
                rows = self.__generatePOST(self.options.par, self.par, self.conSpecs, self.conSpecsTimes, self.options.outdir)
                conspecData.extend(rows)
                outstr.append( Model.FIELD_DELIM.join( map(str, rows) ) )
            if randSpecs:
                rows = self.__generatePOST(self.options.par, self.par, self.randSpecs, self.randSpecsTimes, self.options.outdir)
                randomData.extend(rows)
                outstr.append( Model.FIELD_DELIM.join( map(str, rows) ) )

	    print >> hyperstatsOut, Model.FIELD_DELIM.join( [indx] + computeStats(len(conSpecs), self.obsCnt, conspecData, randomData) )
            if runDatOut:
                print >> runDatOut, Model.FIELD_DELIM.join( [indx] + outstr)
	
	

def prepareNewParFile(obs, parData, time, LPType, PopType, modifyTime = True):
    """ populate the par object with the correct values.  Also modify the timestamp base on data from obs file """
    if modifyTime:
	chngtime = str( int( time / float(obs.gen)) )
    else:
	chngtime = str(time)
    par = copy.copy(parData)
    par.setPopulation(PopType, obs.getPopRange())
    par.setTime(chngtime)
    par.setLociRate(LPType, obs.getMutationRange())
    par.setgamma(obs.gamma)
    par.setSampleSize(obs.nsam)
    par.setLoci(obs.nsites)
    par.setTSTV(obs.tstv)
    return chngtime, par


def parseObs(obs):
    """
    The observation file is a tab delimited file. This can easily be
    split and stored in memory.
    default ordering
    >species,nsam,nsites,tstv,gamma,gen,locuslow,locushigh,Nelow,Nehigh,SegSites,nucdiv,Haptypes,HapDiver,PairDiffs,TajimasD,F*,ExpHet   
    """

    obsf = open(obs, "rU")
    ObservationData.columns = map(string.strip, obsf.next().strip().lower().split("\t"))
    obsl = [ObservationData( dict( izip(ObservationData.columns, l.strip().split("\t")) ) ) for l in obsf]
    obsf.close()
    return obsl


def is_exe(fpath):
    # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """ Validate that the user provided path, does infact exist for BayeSSC """
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


def main_init(options, par):
    """
    main loop specific to the initial mode of the program
    """
    observations = parseObs(options.obs)
    obsCnt = len(observations)
    if options.makestats:
        obsStats = open(os.path.join(options.outdir,"hyperstats_observations.txt"), "w")
        index = "%s_%s_%s_%s_%s"%(options.uid, -1, -1, -1, "_".join([str(random.random()), str(time.time())]).replace(".","_"))
        print >> obsStats, Model.FIELD_DELIM.join( [index] + computeStats(0, obsCnt, obsData = observations) )
	obsStats.close()
    
    hyperstats = open(os.path.join(options.outdir, "hyperstats_iterations_%s.txt"%(options.repeats)), "w")
    runData = None
    if not options.onlyHyperstats:
	    runData = open(os.path.join(options.outdir, "run_data_iterations_%s.csv"%(options.repeats)), "w")
    processor = Model(options, par, observations, obsCnt, ObservationSplitter("uniform"), TimeGenerator("uniform"), BayeSSC(options.bayesPath) )       
    if options.model == None:
	for modelNum in xrange(obsCnt + 1):	    
	    processor.execute(modelNum, hyperstats, runData)
	    if hyperstats:
		hyperstats.flush()
	    if runData:
		runData.flush()
    else:
	processor.execute(options.model, hyperstats, runData)

    hyperstats.close()
    if runData:
	runData.close()	


def selectRuns(uidlst, run_dat):
    uids = dict([(l.strip().split()[0], None,) for l in open(uidlst, "rU")])
    rundat = open(run_dat, "rU")
    parsedData = []
    for l in rundat:
	line = l.strip().split()
	if line[0] not in uids:
	    continue
	uid = line[0]
	line = line[1:]
	start = 0
	end = 14
	cnt = len(line) / 14
	step = 14
	obs = []
	for x in xrange(cnt):
	    obs.append(line[start:end])
	    start = end
	    end +=  step
	model = uid.split("_")[2]
	parsedData.append([ model, obs[:model], obs[model:] ] )
    return parsedData


def main_post(options, par):
    """
    main loop specific to the posterior mode of the program
    """
    observations = parseObs(options.obs)
    #TODO: parse the run_data and the UID list to select what to process
    conspecs, randspecs = selectRuns(options.uidlst, options.run_dat)
    
    
    obsCnt = len(observations)    
    hyperstats = open(os.path.join(options.outdir, "hyperstats_iterations_%s.txt"%(options.repeats)), "w")
    runData = None
    if not options.onlyHyperstats:
	    runData = open(os.path.join(options.outdir, "run_data_iterations_%s.csv"%(options.repeats)), "w")
    #index = "%s_%s_%s_%s_%s"%(options.uid, -1, -1, -1, "_".join([str(random.random()), str(time.time())]).replace(".","_"))
    processor = PostModel(options, par, conSpecs, conSpecsTimes, randSpecs, randSpecsTimes, BayeSSC(options.bayesPath) )          
    processor.execute(model, hyperstats, runData)

    hyperstats.close()
    if runData:
	runData.close()	


def main():
    """
    Main loop of the appliocation
    drives how the program executes (only 1 model, or multiple models).
    """
    options = commandlineArgs()
    par = ParFile(options.par)
    if options.mode == 'initial':
	main_init(options, par)
    elif options.mode == 'posterior':
	main_post(options, par)
    else:
	pass



def mode_init(parser, options, args):
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
    return (options, args,)
    

def mode_post(parser, options, args):
    if not options.uidlst:
	parser.print_help()
	parser.error("UID list file is required")
    if not options.run_dat:
	parser.print_help()
	parser.error("Run data file is required")  
    return (options, args,)
    


	
def commandlineArgs():
    """
    Command line arguments for this script.  It also provides some validation on the arguments passed in, to make sure we have some type of chance to actually succeede.
    """
    global BAYESSC_PATH
    parser = OptionParser("%prog [options]")

    parser.add_option("", "--mode", dest = "mode", help = "program operation mode [required]", action = "store", type = "string", choices = [ 'initial', 'posterior' ] )
    parser.add_option("-p", "--par", dest = "par", help = "par file template [required]", action = "store", type = "string", metavar = "FILE")
    parser.add_option("-i", "--obs", dest = "obs", help = "Observation file [required]", action = "store", type = "string", metavar = "FILE")
    parser.add_option("-r", "--repeat", dest = "repeats", help = "Number of times to try a given congruent group size [required]", action = "store", type = "int", metavar = "NUM")
    parser.add_option("-u", "--uid", dest = "uid", help = "Unique ID to prefix generated indices [required]", action = "store", type = "string", metavar = "UID")
    parser.add_option("-b", "--bayepath", dest = "bayesPath", help = "Path to BayeSSC application [default: Located on user PATH]", action = "store", type = "string", metavar = "PATH", default = "BayeSSC")
    parser.add_option("", "--only_hyperstats", action="store_true", dest="onlyHyperstats", default=False, help="When set, will only generate the hyperstats file")
    parser.add_option("", "--print_headers", action="store_true", dest="headers", default=False, help="When set will generate a headers.txt and exit")
    parser.add_option("-o", "--outdir", dest = "outdir", help = "Directory to generate final outputs in (will create missing folders) [default: %default]", action = "store", type = "string", metavar = "PATH", default = os.getcwd())


    init_group = OptionGroup(parser, "Regular Run", "Options to be applied during mode 'initial'")
    
    init_group.add_option("-m", "--model", dest = "model", help = "Run a single model (0 to total entries in observation file) [default: run all models] ", action = "store", type = "int", metavar = "MODEL", default = None)
    init_group.add_option("-l", "--LPType", dest = "LPType", help = "Loci Rate Priori Type", action = "store", type = "choice", choices = ["U"], default = "U", metavar = "TYPE")
    init_group.add_option("-t", "--timerange", dest= "trange", help = "The range of values to select the time from (Integers). Example: 1000:20000  [required]", action = "store", type = "string", metavar ="RANGE")
    init_group.add_option("", "--obs_stats", action="store_true", dest="makestats", default=False, help="When set, will generate a statistics output for the observation data")

    parser.add_option_group(init_group)    

    post_group = OptionGroup(parser, "Posterior Run", "Options to be applied during mode 'posterior'")   

    post_group.add_option("", "--uid_list", action="store", dest="uidlst", default="", type = "string", metavbar = "FILE", help="Speccifies a list of UIDs to filter on for Posterior processing [required]")
    post_group.add_option("", "--run_data", action="store", dest="run_dat", default="", type = "string", metavar = "FILE", help="run data which contains the --uid_list UIDs.  It is used for the Posterior processing [required]")

    parser.add_option_group(post_group)    

    (options, args) = parser.parse_args()    

    if options.headers:
	 print "Generating header file"
	 bayshdr = ['index'] + BayeSSCData().header()
	 hyperhdr = ['index'] + statsHeader()
	 o = open("headers.txt", "w")
	 print >> o , "Iteration file header"
	 print >> o, Model.FIELD_DELIM.join(bayshdr)
	 print >> o, "\nHyperstats file header"
	 print >> o, Model.FIELD_DELIM.join(hyperhdr)
	 o.close()
	 sys.exit()


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

    if options.mode == 'initial':
	options, args = mode_init(parser, options, args)
    elif options.mode == 'posterior':
	options, args = mode_post(parser, options, args)
    else:
	parser.print_help()
	parser.error("Mode must be either 'initial' or 'posterior'")
	
    return options, 


if __name__ == "__main__":
    main()
