import math,random
import fastatomfm3 as atomfm
import myfeedforward4 as myfeedforward


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def read_fei_file(filename):
   """Lazy function (generator) to read a fei_file 
      frame by frame."""
   #tobohr = 1.0/0.529177
   tobohr = 1.0  #fei file is in atomic units
   with open(filename,'r') as ff: 
      xyzdata = ff.readline()
   nion = int(xyzdata.split('\n')[0])
   nlines = nion+5
   #print("NION=",nion)

   with open(filename,'r') as fp:
      while True:
         data = ''.join([fp.readline() for _ in range(nlines)])
         if not data: 
            break
         lines = data.strip().split('\n')
         energy  = float(lines[1])
         symbols = [ln.split()[1] for ln in lines[5:]]
         rions   = [tobohr*float(ln.split()[i]) for ln in lines[5:] for i in range(3,6)]
         fions   = [tobohr*float(ln.split()[i]) for ln in lines[5:] for i in range(6,9)]

         yield (symbols,rions,fions,energy)


class mlpotential_e:
   """implementation of simple ML potential
   """

   def __init__(self,nsweeps,nlayers,parameterfilename,feidatafile,energygradient0):
      self.energygradient0 = energygradient0

      parameters = []
      with open(parameterfilename,'r') as ff:
         data = ff.read()
         for line in data.split('\n'):
            ss = line.split()
            if (len(ss)>1):
               p = []
               for s in ss:
                  if s.isdigit():
                     p += [int(s)] 
                  elif isFloat(s): 
                     p += [float(s)]
               ok = p[0] in range(0,6)
               if (p[0]==0): ok = ok and (len(p)==2)
               if (p[0]==1): ok = ok and (len(p)==2)
               if (p[0]==2): ok = ok and (len(p)==4)
               if (p[0]==3): ok = ok and (len(p)==3)
               if (p[0]==4): ok = ok and (len(p)==5)
               if (p[0]==5): ok = ok and (len(p)==5)
               if ok: parameters.append(p)

      #### define Atomic Feature Mapping ####
      nparameters = len(parameters)
      #self.nparameters = nparameters
      self.afm = atomfm.AtomsFM(parameters)

      #### define NN machine ####
      self.nlayers = nlayers
      #sigmoid      = lambda x: 1.0/(1.0+math.exp(-x))
      #sigmoidp     = lambda x: math.exp(-x)/(math.exp(-x)+1.0)**2
      #sigmoidpp    = lambda x: 2*math.exp(-2*x)/(math.exp(-x)+1.0)**3 - math.exp(-x)/(math.exp(-x)+1.0)**2

      sigmoid      = lambda x: math.tanh(x)
      sigmoidp     = lambda x: (1.0/math.cosh(x))**2
      sigmoidpp    = lambda x: -2.0*math.tanh(x)*(1.0/math.sech(x))**2
      #anparameters = [nparameters]
      #asigmoid   = [lambda x: x]
      #asigmoidp  = [lambda x: 1]
      #asigmoidpp = [lambda x: 0]
      anparameters = []
      asigmoid   = []
      asigmoidp  = []
      asigmoidpp = []
      for i in range(nlayers):
         anparameters.append(nparameters)
         asigmoid.append(sigmoid)
         asigmoidp.append(sigmoidp)
         asigmoidpp.append(sigmoidpp)
      anparameters.append(1)
      asigmoid.append(lambda x: x)
      asigmoidp.append(lambda x: 1)
      asigmoidpp.append(lambda x: 0)
      self.nn_machine = myfeedforward.MyFeedForward(anparameters,asigmoid,asigmoidp,asigmoidpp)
 
      #### read in number of atoms ####
      with open(feidatafile,'r') as ff:
         feidata = ff.readline()
      nion = int(feidata)
      self.nion = nion
      print("nion=",nion)

      #### define NN machine weights for all atoms ####
      self.nn_weights  = []
      for ii in range(nion):
         self.nn_weights += self.nn_machine.initial_w()
      self.nw = len(self.nn_weights)/nion

      alpha = 0.05
      for (symbols,rions,fions,energy) in read_fei_file(feidatafile):
         #aalpha = alpha*random.random()
         aalpha = alpha*random.random()
         etmp = []
         dedw = []
         for ii in range(nion):
            fm = self.afm(rions,ii)
            print "fm=",fm
            eee = self.nn_machine.dyoutdw_gradient(fm,self.nn_weights[ii*self.nw:(ii+1)*self.nw])
            etmp += eee[0]
            dedw += eee[1]
            #print "etmp=",ii,eee[0],energy

         error = math.sqrt((sum(etmp) - energy)**2)
         derror1detmp    = 2.0*(sum(etmp) - energy)

         for i in range(len(self.nn_weights)):
            self.nn_weights[i] -= aalpha*derror1detmp*dedw[i]

      self.nn_machine.print_w(self.nn_weights[0])

      print("Checking Energies and Forces")
      #nion3 = 3*nion
      frame = 1
      sumerror = 0.0
      maxerror = 0.0
      for (symbols,rions,fions,energy) in read_fei_file(feidatafile):
         etmp0 = []
         #force3 = [0.0]*nion3
         for ii in range(nion):
            fm   = self.afm(rions,ii)
            ee0 = self.nn_machine.evaluate(fm,self.nn_weights[ii*self.nw:(ii+1)*self.nw])
            etmp0 += ee0

            #fafm = self.afm.Egradients(rions,ii)
            #eee = self.nn_machine.evaluate(fafm[0],self.nn_weights[ii*self.nw:(ii+1)*self.nw])
            #fff = self.nn_machine.gradients_evaluate(fafm[0],self.nn_weights[ii*self.nw:(ii+1)*self.nw])

            #esum  += eee[0]
            #for jj in range(nion3):
            #   for k in range(nparameters):
            #      force3[jj] -= fafm[1][jj + k*nion3]*fff[k]

         error = math.sqrt((sum(etmp0) - energy)**2)
         print frame,sum(etmp0),energy,error
         sumerror += error
         if (error>maxerror): maxerror = error
         frame += 1

      print(" - average error=", sumerror/frame, " maxerror=",maxerror)



