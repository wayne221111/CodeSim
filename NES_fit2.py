import sys
from scipy import eye, multiply, ones, dot, array, outer, rand, zeros, diag, randn, exp
from scipy.linalg import cholesky, inv, det

import math
import pybrain
from pybrain.optimization import ExactNES
import os
os.nice(20)
scale=1
offset=0
useNearest = True
globel_Dsize = 0
if(len(sys.argv) == 1):
	print "Usage: NES_fit.py filename"
	sys.exit()

def parseInt(a):
	b=[]
	for i in range(len(a)):
		b.append(int(a[i]))
	return b

def parseFloat(a):
	b=[]
	for i in range(len(a)):
		b.append( float(a[i]) )
	return b

def toStringList(a):
	b=[]
	for i in range(len(a)):
		b.append( `a[i]` )
	return b

def add(x,y):
	return x+y

def normalize(p):
	from scipy import sum
	sum_p=float(sum(p))
	p.append(1-sum_p)

def mygetline(file):
	while(True):
		s=file.readline()
		if(s == ''):
			return s
		if(s[0] != '#'):
			return s

def log10(x):
	return math.log(x, 10)

def abs(x):
	return math.fabs(x)

def outOfRange(x):
	if(x>1):
		return x-1
	if(x<0):
		return -x
	return 0

def validDist(p):
	# check if p is a valid probability distribution
	sum = 0.0
	for i in p[:lt.Dsize-1]:
		if i < 0:
			return False
		sum += i
	if sum > 1:
		return False
	
	# check mean and low_degree_probability
	degree=lt.Tags[:3]+map(int, p[lt.Dsize-1:]+0.5)
	probability=map(float, p[:lt.Dsize-1])
	normalize(probability)
	
	mean=0.0
	low_degree_probability=0.0
	for i in range(lt.Dsize):
		mean += degree[i]*probability[i]
		if degree[i] <= 10:
			low_degree_probability += probability[i]
	
	
	if mean > 25:
		return False
	if low_degree_probability < 0.6:
		return False
	return True

def nearestPoint(p):
	sum = 0.0
	shift = 0.0
	positiveDim=0
	for i in p:
		if(i>0):
			sum += i
			positiveDim = positiveDim + 1
	if(sum>1):
		shift = (sum-1) / positiveDim
	def trans(x):
		if(x>0):
			return x-shift
		else:
			return 0.0
	if(sum>1):
		return nearestPoint(map(trans, p))
	else:
		return map(trans, p)

def degree_bound(x):
	if x<4:
		return 4.0
	if x>lt.K:
		return float(lt.K)
	return x

class LT_exp:
	def __init__(self, filename):
		self.filename = filename
		#read from config file
		f=open(filename,'r')
		tmp=mygetline(f).split()
		self.K = int(tmp[0])
		self.Run = int(tmp[1])
		self.Dsize = int(mygetline(f))
		globel_Dsize = self.Dsize
		self.Tags = parseInt(mygetline(f).split())
		if(self.Dsize != len(self.Tags)):
			print "Error: Dsize and Tags doesn't match"
			sys.exit()
		self.D = parseFloat(mygetline(f).split())
		if(self.Dsize != len(self.D)):
			print "Error: Dsize and initial distribution doesn't match"
			sys.exit()
		tmp=mygetline(f).split()
		self.STEPS = int(tmp[0])
		self.Delta = float(tmp[1])
		self.targetRho = float(mygetline(f))
		self.targetFailureRate = float(mygetline(f))
		self.targetEpsilon = float(mygetline(f))
		self.optimParameter = int(mygetline(f))
		self.maxs = [0,0,0,self.STEPS*self.Delta]
		# prepare LT_BER process
		LT_BER_format = "{K} {Run}\n{Dsize}\n{tag_list}\n{STEPS} {Delta}\n{targetRho}\n{targetFailureRate}\n{targetEpsilon}\n{optimParameter}\n"

		from subprocess import Popen
		from subprocess import PIPE
		self.p = Popen('./LT_BER2.out', stdin=PIPE, stdout=PIPE)
		input = LT_BER_format.format(K=self.K, Run=self.Run, Dsize=self.Dsize, tag_list=' '.join(toStringList(self.Tags)), STEPS=self.STEPS, Delta=self.Delta, targetRho = self.targetRho, targetFailureRate= self.targetFailureRate, targetEpsilon= self.targetEpsilon, optimParameter = self.optimParameter)
		self.p.stdin.write(input)
		
		#prepare listener for ExactNES to write out results
		self.file_result = open('%s_result.txt' % filename, 'w')
		self.fitness_log = open('%s_fitness_log.txt' % filename, 'w')
		self.result_line = 0
		self.times_of_eval = 0

	def fitness_LT(self, dist):
		self.times_of_eval += 1
		
		d=map(float,dist[:self.Dsize-1]) # distribution 
		normalize(d) # calculate the probability of dependent degree
		
		t=self.Tags[:3]+map(int, dist[self.Dsize-1:]+0.5) #tabs of degree

		

		print 'O',
		if(d[0]==0):
			return self.maxs[self.optimParameter]
		input=' '.join(toStringList(t))+'\n'+' '.join(toStringList(d))
		#K='`self.K`', Run='100', Dsize='`self.Dsize`', tag_list=' '.join(self.Tags), distribution='`self.D`', STEPS='16', epsilon_list= "")
		#print input
		self.p.stdin.write('%s\n' % input)
		fit = float(self.p.stdout.readline())

		
		return fit
	def printResult(self, parameter_list, eval):
		self.result_line += 1
		p_list = map(float,parameter_list[:self.Dsize-1])
		normalize(p_list)
		self.file_result.write(`self.result_line`+'\t'+ `self.times_of_eval`+'\teval\t'+`-eval`+'\tpara\t' + '\t'.join(toStringList(p_list))+'\tdegree\t'+'\t'.join(toStringList(self.Tags[:3]+map(int,parameter_list[self.Dsize-1:]+0.5)))+'\n')
		self.file_result.flush()
		self.fitness_log.write(`-eval`+'\t')
		self.fitness_log.flush()
		dist = open('%s_dist.txt' % filename, 'w')
		dist.write(`self.K` + '\n' + `self.Dsize` + '\n' + '\t'.join(toStringList(self.Tags[:3]+map(int,parameter_list[self.Dsize-1:]+0.5))) + '\n' + '\t'.join(toStringList(p_list))+'\n')
		dist.close()

class ExactNESforLT(ExactNES):
    def _additionalInit(self):
        xdim = self.numParameters
        assert not self.diagonalOnly, 'Diagonal-only not yet supported'
        self.numDistrParams = xdim + xdim * (xdim + 1) / 2
                
        if self.momentum != None:
            self.momentumVector = zeros(self.numDistrParams)
        if self.learningRateSigma == None:
            self.learningRateSigma = self.learningRate
        
        if self.rangemins == None:
            self.rangemins = -ones(xdim)
        if self.rangemaxs == None:
            self.rangemaxs = ones(xdim)
        if self.initCovariances == None:
            if self.diagonalOnly:
                self.initCovariances = ones(xdim)
            else:
                self.initCovariances = eye(xdim)

        #self.x = rand(xdim) * (self.rangemaxs - self.rangemins) + self.rangemins
        self.x = self._initEvaluable
        self.sigma = dot(eye(xdim)*.000625, self.initCovariances)
        for i in range(lt.Dsize-1, 2*lt.Dsize-1-3):
        	self.sigma[i][i]=100
        self.factorSigma = cholesky(self.sigma)
        
        # keeping track of history
        self.allSamples = []
        self.allFitnesses = []
        self.allPs = []
        
        self.allGenerated = [0]
        
        self.allCenters = [self.x.copy()]
        self.allFactorSigmas = [self.factorSigma.copy()]
        
        # for baseline computation
        self.phiSquareWindow = zeros((self.batchSize, self.numDistrParams))
        
        
    def _produceNewSample(self, z=None, p=None):
        if z == None:
            while True:
                p = randn(self.numParameters)
                z = dot(self.factorSigma.T, p) + self.x
                z = array(map(float,z[:lt.Dsize-1])+map(degree_bound,z[lt.Dsize-1:]))
                if useNearest:
                    z = array(nearestPoint(z[:lt.Dsize-1])+map(float,z[lt.Dsize-1:]))
                if validDist(z):
                    break
        if p == None:
            p = dot(inv(self.factorSigma).T, (z - self.x))            
        self.allPs.append(p)
        self.allSamples.append(z)
        fit = self._oneEvaluation(z)
        self.allFitnesses.append(fit) 
        return z, fit
        
    def _produceSamples(self):
        """ Append batchsize new samples and evaluate them. """
        if self.numLearningSteps == 0 or not self.importanceMixing:
            for _ in range(self.batchSize):
                self._produceNewSample()
            self.allGenerated.append(self.batchSize + self.allGenerated[-1])
        else:
            olds = len(self.allSamples)
            oldDetFactorSigma = det(self.allFactorSigmas[-2])
            newDetFactorSigma = det(self.factorSigma)
            invA = inv(self.factorSigma)
    
            # All pdfs computed here are off by a coefficient of 1/power(2.0*pi, self.numDistrParams/2.)
            # but as only their relative values matter, we ignore it.
            
            # stochastically reuse old samples, according to the change in distribution
            for s in range(olds - self.batchSize, olds):
                oldPdf = exp(-0.5 * dot(self.allPs[s], self.allPs[s])) / oldDetFactorSigma
                sample = self.allSamples[s]
                newPs = dot(invA.T, (sample - self.x))
                newPdf = exp(-0.5 * dot(newPs, newPs)) / newDetFactorSigma
                r = rand()
                if r < (1 - self.forcedRefresh) * newPdf / oldPdf:
                    self.allSamples.append(sample)
                    self.allFitnesses.append(self.allFitnesses[s])
                    self.allPs.append(newPs)
                # never use only old samples
                if (olds + self.batchSize) - len(self.allSamples) < self.batchSize * self.forcedRefresh:
                    break
            self.allGenerated.append(self.batchSize - (len(self.allSamples) - olds) + self.allGenerated[-1])

            # add the remaining ones
            oldInvA = inv(self.allFactorSigmas[-2])
            while  len(self.allSamples) < olds + self.batchSize:
                r = rand()
                if r < self.forcedRefresh:
                    self._produceNewSample()
                else:
                    while True:
                        p = randn(self.numParameters)
                        newPdf = exp(-0.5 * dot(p, p)) / newDetFactorSigma
                        sample = dot(self.factorSigma.T, p) + self.x
                        sample = array(map(float,sample[:lt.Dsize-1])+map(degree_bound,sample[lt.Dsize-1:]))
                        if useNearest:
                            sample = array(nearestPoint(sample[:lt.Dsize-1])+map(float,sample[lt.Dsize-1:]))
                        if validDist(sample):
                            break
                    oldPs = dot(oldInvA.T, (sample - self.allCenters[-2]))
                    oldPdf = exp(-0.5 * dot(oldPs, oldPs)) / oldDetFactorSigma
                    if r < 1 - oldPdf / newPdf:
                        self._produceNewSample(sample, p)


"""main()"""
filename = sys.argv[1]
lt = LT_exp(filename)
# out = dup(sys.stdout, open('result_%s' % filename, 'w'))
# sys.stdout = out
# writer = ResultWriter('dist_%s' % filename)
nes = ExactNESforLT(lt.fitness_LT, array( lt.D[:lt.Dsize-1]+map(float,lt.Tags[3:]) ), minimize=True, maxEvaluations=10000, verbose=True, listener=lt.printResult)
nes_result = nes.learn()
final = map(float, nes_result[0])
normalize(final)
final = toStringList(final)
print ' '.join(final)
# dist = open('dist_%s' % filename, 'w')
# dist.write(' '.join(final) + '\n')
# dist.close()
# p = pow(2).power
# from pybrain.optimization import ExactNES
# print ExactNES(pow(32.60).power, [0], maxEvaluations = 10000, minimize=True).learn()[0][0]
