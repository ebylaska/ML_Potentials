
class AtomFunctions:
   """AtomFunctions - defines functions for Atom Feature mapping

   The following functions are defined
      1) f1(rion,i) = Sum(j!=i) fcut(rij)
      2) f2(rion,i) = Sum(j!=i) exp(-eta*(rij-rs)**2)*fcut(rij)
      3) f3(rion,i) = Sum(j!=i) cos(kappa*rij)*fcut(rij)
      4) f4(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik + rjk*rjk))*fcut(rij)*fcut(rik)*fcut(rjk)
      5) f5(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik))*fcut(rij)*fcut(rik)

   Defining functions  
      f1 = AtomFunctions([1,rcut])
      f2 = AtomFunctions([2,rcut,eta,rs])
      f3 = AtomFunctions([3,rcut,kappa])
      f4 = AtomFunctions([4,rcut,eta,lmbda,xi])
      f5 = AtomFunctions([4,rcut,eta,lmbda,xi])

   Using functions
      value1 = f1(rion,i)
      value2 = f2(rion,i)
      value3 = f3(rion,i)
      value4 = f4(rion,i)
      value5 = f5(rion,i)
 
   """

   def __fcut(self,rcut,r):
      import math
      ff = 0.0
      if r < rcut:
         ff = 0.5*(math.cos(math.pi*r/rcut)+1.0)
      return ff

   def __f2exp(self,eta,rs,r):
      import math
      return math.exp(-eta*(r-rs)*(r-rs))

   def __f2cos(self,kappa,r):
      import math
      return math.cos(kappa*r)

   def __dist(self,rion,i,j):
      import math
      x = rion[3*j]   - rion[3*i]
      y = rion[3*j+1] - rion[3*i+1]
      z = rion[3*j+2] - rion[3*i+2]
      return math.sqrt(x*x + y*y + z*z)

   def __theta(self,rion,i,j,k):
      import math
      if (j==k):
         vv = 0.0
      else:
         xij = rion[3*i]   - rion[3*j]
         yij = rion[3*i+1] - rion[3*j+1]
         zij = rion[3*i+2] - rion[3*j+2]
         xik = rion[3*i]   - rion[3*k]
         yik = rion[3*i+1] - rion[3*k+1]
         zik = rion[3*i+2] - rion[3*k+2]
         rij = math.sqrt(xij*xij+yij*yij+zij*zij)
         rik = math.sqrt(xik*xik+yik*yik+zik*zik)
         ff =(xij*xik + yij*yik + zij*zik)/(rij*rik)
         if (ff>1.0):    ff = 1.0
         if (ff<(-1.0)): ff = -1.0
         vv = math.acos(ff)
      return vv

   def __init__(self,parameters):
      self.__ftype = parameters[0]
      self.__parameters = parameters[1:]
      return

   def __call__(self,rion,i):
      import math
      nion = len(rion)/3
      ff = 0.0
      if self.__ftype==1:
         rcut = self.__parameters[0]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__fcut(rcut,r)
      elif self.__ftype==2:
         rcut = self.__parameters[0]
         eta  = self.__parameters[1]
         rs   = self.__parameters[2]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__f2exp(eta,rs,r)*self.__fcut(rcut,r)
      elif self.__ftype==3:
         rcut  = self.__parameters[0]
         kappa = self.__parameters[1]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__f2cos(kappa,r)*self.__fcut(rcut,r)
      elif self.__ftype==4:
         rcut  = self.__parameters[0]
         eta   = self.__parameters[1]
         lmbda = self.__parameters[2]
         xi    = self.__parameters[3]
         for j in range(nion):
            if (i!=j):
               rij  = self.__dist(rion,i,j)
               rij2 = rij*rij
               for k in range(nion):
                  if (i!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     theta = self.__theta(rion,i,j,k)
                     tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp2 = math.exp(-eta*(rij2 + rik2 + rjk2))
                     ff += tmp1*tmp2*self.__fcut(rcut,rij)*self.__fcut(rcut,rik)*self.__fcut(rcut,rjk)
         ff *= 2.0**(1-xi)
      elif self.__ftype==5:
         rcut  = self.__parameters[0]
         eta   = self.__parameters[1]
         lmbda = self.__parameters[2]
         xi    = self.__parameters[3]
         for j in range(nion):
            if (i!=j):
               rij  = self.__dist(rion,i,j)
               rij2 = rij*rij
               for k in range(nion):
                  if (i!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     theta = self.__theta(rion,i,j,k)
                     tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp2 = math.exp(-eta*(rij2 + rik2))
                     ff += tmp1*tmp2*self.__fcut(rcut,rij)*self.__fcut(rcut,rik)
         ff *= 2.0**(1-xi)
      else:
         print("no function defined, ftype=",self.__ftype)
     
      return ff

class AtomsFM:

   def __init__(self,parameterslist):
      self.atomfunctions = []
      for parameters in parameterslist:
         self.atomfunctions.append(AtomFunctions(parameters))
      return

   def __call__(self,rion,i):
      fm = []
      for f in self.atomfunctions:
         fm.append(f(rion,i))
      return fm



#### examples of defining and using AtomFunctions ####
f1 = AtomFunctions([1,6.0/0.529177])
f2 = AtomFunctions([2,6.0/0.529177,0.01,0.00])
f3 = AtomFunctions([3,6.0/0.529177,0.21])
f4 = AtomFunctions([4,6.0/0.529177,0.01,1.00,1.0])
f5 = AtomFunctions([5,6.0/0.529177,0.01,1.00,1.0])
print "f1=",f1
print "f2=",f2
print "f3=",f3
print "f4=",f4
print "f5=",f5

rion = [0.0, 0.0, 0.0]
rion += [0.0, 0.0, 2.0]
rion += [0.0, 0.0, -2.0]
rion += [0.0, 2.0, 2.0]
rion += [-2.0, 2.0, 2.0]

print "f1(rion,1)=",f1(rion,1)
print "f2(rion,1)=",f2(rion,1)
print "f3(rion,1)=",f3(rion,1)
print "f4(rion,1)=",f4(rion,1)
print "f5(rion,1)=",f5(rion,1)

#### example of defining and using AtomFM (atom feature mapping) ####
allparameters = [[1,6.0/0.529177],
                 [2,6.0/0.529177,0.01,0.00],
                 [3,6.0/0.529177,0.21],
                 [4,6.0/0.529177,0.01,1.00,1.0],
                 [5,6.0/0.529177,0.01,1.00,1.0]]

afm = AtomsFM(allparameters)
print "afm(rion,1)=",afm(rion,1)


#### Larger example of defining and using AtomFM (atom feature mapping) ####
#### Using the parameter values from Nongnuch Artrith and Jorg Behler,PHYSICAL REVIEW B 85, 045439 (2012) ####
allparameters2 = []
allparameters2 += [[2,11.338,0.001,0.0]]
allparameters2 += [[2,11.338,0.010,0.0]]
allparameters2 += [[2,11.338,0.020,0.0]]
allparameters2 += [[2,11.338,0.035,0.0]]
allparameters2 += [[2,11.338,0.060,0.0]]
allparameters2 += [[2,11.338,0.100,0.0]]
allparameters2 += [[2,11.338,0.200,0.0]]
allparameters2 += [[2,11.338,0.400,0.0]]

allparameters2 += [[4,11.338,0.0001,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0001, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0001,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0001, 1.0,2.0]]

allparameters2 += [[4,11.338,0.0030,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0030, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0030,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0030, 1.0,2.0]]

allparameters2 += [[4,11.338,0.0080,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0080, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0080,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0080, 1.0,2.0]]

allparameters2 += [[4,11.338,0.0150,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0150, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0150,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0150, 1.0,2.0]]
allparameters2 += [[4,11.338,0.0150,-1.0,4.0]]
allparameters2 += [[4,11.338,0.0150, 1.0,4.0]]
allparameters2 += [[4,11.338,0.0150,-1.0,16.0]]
allparameters2 += [[4,11.338,0.0150, 1.0,16.0]]

allparameters2 += [[4,11.338,0.0250,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0250, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0250,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0250, 1.0,2.0]]
allparameters2 += [[4,11.338,0.0250,-1.0,4.0]]
allparameters2 += [[4,11.338,0.0250, 1.0,4.0]]
allparameters2 += [[4,11.338,0.0250,-1.0,16.0]]
allparameters2 += [[4,11.338,0.0250, 1.0,16.0]]

allparameters2 += [[4,11.338,0.0450,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0450, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0450,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0450, 1.0,2.0]]
allparameters2 += [[4,11.338,0.0450,-1.0,4.0]]
allparameters2 += [[4,11.338,0.0450, 1.0,4.0]]
allparameters2 += [[4,11.338,0.0450,-1.0,16.0]]
allparameters2 += [[4,11.338,0.0450, 1.0,16.0]]

allparameters2 += [[4,11.338,0.0800,-1.0,1.0]]
allparameters2 += [[4,11.338,0.0800, 1.0,1.0]]
allparameters2 += [[4,11.338,0.0800,-1.0,2.0]]
allparameters2 += [[4,11.338,0.0800, 1.0,2.0]]
allparameters2 += [[4,11.338,0.0800,-1.0,4.0]]
allparameters2 += [[4,11.338,0.0800, 1.0,4.0]]
allparameters2 += [[4,11.338,0.0800,-1.0,16.0]]
allparameters2 += [[4,11.338,0.0800, 1.0,16.0]]

afm2 = AtomsFM(allparameters2)
print "afm2(rion,0)=",afm2(rion,0),len(afm2(rion,0))
print "afm2(rion,1)=",afm2(rion,1),len(afm2(rion,1))
print "afm2(rion,2)=",afm2(rion,2),len(afm2(rion,2))
print "afm2(rion,3)=",afm2(rion,3),len(afm2(rion,3))
print "afm2(rion,4)=",afm2(rion,4),len(afm2(rion,4))
