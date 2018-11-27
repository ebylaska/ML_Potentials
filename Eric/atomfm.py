import math

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
      ff = 0.0
      if r < rcut:
         ff = 0.5*(math.cos(math.pi*r/rcut)+1.0)
      return ff

   def __f2exp(self,eta,rs,r):
      return math.exp(-eta*(r-rs)*(r-rs))

   def __f2cos(self,kappa,r):
      return math.cos(kappa*r)

   def __dfcutdr(self,rcut,r):
      ff = 0.0
      if r < rcut:
         ff = -0.5*math.pi/5.0*math.sin(math.pi*x/5.0)
      return ff

   def __df2expdr(self,eta,rs,r):
      return -2.0*eta*(r-rs)*math.exp(-eta*(r-rs)*(r-rs)) 

   def __df2cosdr(self,kappa,r):
      return -kappa*math.sin(kappa*r)

   def __dist(self,rion,i,j):
      x = rion[3*j]   - rion[3*i]
      y = rion[3*j+1] - rion[3*i+1]
      z = rion[3*j+2] - rion[3*i+2]
      return math.sqrt(x*x + y*y + z*z)

   def __theta(self,rion,i,j,k):
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

   def __ctheta(self,rion,i,j,k):
      if (j==k):
         vv = 1.0
      else:
         xij = rion[3*i]   - rion[3*j]
         yij = rion[3*i+1] - rion[3*j+1]
         zij = rion[3*i+2] - rion[3*j+2]
         xik = rion[3*i]   - rion[3*k]
         yik = rion[3*i+1] - rion[3*k+1]
         zik = rion[3*i+2] - rion[3*k+2]
         rij = math.sqrt(xij*xij+yij*yij+zij*zij)
         rik = math.sqrt(xik*xik+yik*yik+zik*zik)
         vv =(xij*xik + yij*yik + zij*zik)/(rij*rik)
         if (vv>1.0):    vv = 1.0
         if (vv<(-1.0)): vv = -1.0
      return vv

   def __init__(self,parameters):
      self.__ftype = parameters[0]
      self.__parameters = parameters[1:]
      return

   def __call__(self,rion,i):
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
                  if (i!=k) and (j!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     ctheta = self.__ctheta(rion,i,j,k)
                     #tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp1 = (1.0-lmbda*ctheta)**xi
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
                  if (i!=k) and (j!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     ctheta = self.__ctheta(rion,i,j,k)
                     #tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp1 = (1.0-lmbda*ctheta)**xi
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


