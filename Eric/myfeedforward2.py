
import numpy as np9

#def matmul(m,n,k,a,b):
#   aa = []
#   bb = []
#   for i in range(m): aa.append([a[i+j*m] for j in range(k)])
#   for i in range(k): bb.append([b[i+j*k] for j in range(n)])
#   cc = np9.matrix(aa)*np9.matrix(bb)
#   return  [cc[i,j] for j in range(n) for i in range(m)]

class MyFeedForward2(object):
   """A feed forward network class
    following properties:

    Attributes:
       
   """
   def __init__(self,nodes,fs,fps,fpps,biases=[]):
      self.nlayers = len(nodes)
      self.nodes   = nodes
      self.f       = fs
      self.fp      = fps
      self.fpp     = fpps
      self.x = []
      self.y = []
      self.yp = []
      self.ypp = []
      self.biases = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         self.x.append([0.0]*n)
         self.y.append([0.0]*n)
         self.yp.append([0.0]*n)
         self.ypp.append([0.0]*n)
         self.biases.append([0.0]*n)
      for ib in range(len(biases)): self.biases[ib] = biases[ib]
      self.matsize = 0
      self.wshift = []
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         self.wshift.append(self.matsize)
         self.matsize += m*n
         #ww = np9.random.normal(0.0,pow(m,-0.5),(m,n))
         #self.w += [ww[i,j] for j in range(n) for i in range(m)]

   #
   # generate a vector of initial weights to be optimized
   #
   def initial_w(self):
      w = []
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         ww = np9.random.normal(0.0,pow(m,-0.5),(m,n))
         w += [ww[i,j] for j in range(n) for i in range(m)]
      return w

   #
   # local matrix multiply function
   #
   def __matmul(self,m,n,k,a,b):
      aa = []
      bb = []
      for i in range(m): aa.append([a[i+j*m] for j in range(k)])
      for i in range(k): bb.append([b[i+j*k] for j in range(n)])
      cc = np9.matrix(aa)*np9.matrix(bb)
      return  [cc[i,j] for j in range(n) for i in range(m)]


   #
   # Computes A*diag(fp(x))
   #
   def __matfpmullayer(self,m,n,a,layer):
      c = [0.0]*m*n
      for j in range(n):
         #tfp = self.fp[layer](self.x[layer][j])
         tfp = self.yp[layer][j]
         for i in range(m):
            c[i+j*m] = a[i+j*m]*tfp
      return c

   #
   # evaluate the feed forward neural network
   #
   def evaluate(self,xin,w):
      for i in range(self.nodes[0]): 
         self.x[0][i] = xin[i]
         self.y[0][i] = self.f[0](xin[i]+self.biases[0][i])
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         self.x[layer+1] = self.__matmul(m,1,n,w[shift:shift+m*n],self.y[layer])
         for i in range(m): 
            self.y[layer+1][i] = self.f[layer+1](self.x[layer+1][i]+self.biases[layer+1][i])

      return self.y[layer+1]

   #
   # evaluate the gradient wrt x, dyout/dx
   #
   def gradients_evaluate(self,xin,w):
      #...initialize values for x
      yout = self.evaluate(xin,w)

      # define yp
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         for k in range(n):
            self.yp[layer][k] = self.fp[layer](self.x[layer][k]+self.biases[layer][k])

      #...calculate dydx
      nlayers = self.nlayers
      nmax = self.nodes[nlayers-1]
      dydx = [0.0]*nmax*nmax
      for i in range(nmax):
         #dydx[i+i*nmax] = self.fp[nlayers-1](self.x[nlayers-1][i])
         dydx[i+i*nmax] = self.yp[nlayers-1][i]
      for layer in range(nlayers-2,-1,-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         dydx = self.__matmul(nmax,n,m,dydx,w[shift:shift+m*n])
         dydx = self.__matfpmullayer(nmax,n,dydx,layer)

      return dydx


   #
   # calculate the error and gradient of w using energy training
   #
   def w_energy_gradient(self,xtrain,ytrain,w):
      Error     = 0.0
      dErrordw  = [0.0]*self.matsize
      yout = self.evaluate(xtrain,w)
      for i in range(len(ytrain)):
         Error += (yout[i]-ytrain[i])**2
         de =  2.0*(yout[i]-ytrain[i])

         layer = self.nlayers - 1
         tmp0 = de*self.fp[layer](self.x[layer][i]+self.biases[layer][i])
         for layer in range(self.nlayers-2,-1,-1):
            n = self.nodes[layer]
            m = self.nodes[layer+1]
            shift = self.wshift[layer]
            ww = w[shift:shift+m*n]
            for k in range(n):
               for j in range(m):
                  ww[j+k*m] *=  self.fp[layer](self.x[layer][k]+self.biases[layer][k])

            if (layer==self.nlayers-2):
               for k in range(n):
                  dErrordw[shift+i+k*m] += tmp0*self.y[layer][k]
               tmp1 = []
               for j in range(n):
                  xx = tmp0*ww[i+j*m]
                  tmp1.append(xx)
            else:
               for k in range(n):
                  for j in range(m):
                     dErrordw[shift+j+k*m] += tmp1[j]*self.y[layer][k]
               tmp1 = self.__matmul(1,n,m,tmp1,ww)
            
      return (Error,dErrordw)

   #
   # calculate the error and gradient of w using derivative training
   #
   def w_energy_gradient1(self,xtrain,dytraindx,w):
      Error     = 0.0
      dErrordw  = [0.0]*self.matsize
      dyoutdx = self.gradients_evaluate(xtrain,w)
      for i in range(len(dytraindx)):
         Error += (youtdx[i]-dytraindx[i])**2

      ### define ypp and dxdw template
      dxdw = []
      dypdw = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         dxdw.append([0.0]*n)
         dypdw.append([0.0]*n)
         for k in range(n):
            self.ypp[layer][k] =  self.fpp[layer](self.x[layer][k]+self.biases[layer][k])
      
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         for j in range(n):
            for i in range(m):
               # dw/dwij
               wwp = [0.0]*m*n
               wwp[i+j*m] = 1.0

               # define needed dxdw's and dfp/dw's
               dxdw[layer+1] = self.__matmul(m,1,n,wwp,self.y[layer])
               for k in range(m): dypdw[layer+1][k] = self.ypp[layer+1][k]*dxdw[layer+1][k]
               for q in range(layer+2,self.nlayers):
                  nqm = self.nodes[q-1]
                  mqm = self.nodes[q]
                  shiftqm = self.wshift[q-1]
                  wwq = w[shiftq:shiftqm+mqm*nqm]
                  tmp     = self.__matfpmullayer(mqm,nqm,wwq,q-1)
                  dxdw[q] = self.__matmul(mqm,1,nqm,tmp,dxdw[q-1])
                  for k in range(mqm): dypdw[q][k] = self.ypp[q][k]*dxdw[q][k]
                
               # calculate dydxij
               nlayers = self.nlayers
               nmax    = self.nodes[nlayers-1]
               dydxij  = [0.0]*nmax*nmax
               for i in range(nmax):
                  dydxij[k+k*nmax] = self.yp[nlayers-1][k]
               for q in range(nlayers-2,layer,-1):
                  nq = self.nodes[q]
                  mq = self.nodes[q+1]
                  shiftq = self.wshift[q]
                  dydxij = self.__matmul(nmax,nq,mq,dydxij,w[shiftq:shiftq+mq*nq])
                  dydijx = self.__matfpmullayer(nmax,nq,dydxij,q)
               dydxij = self.__matmul(nmax,n,m,dydxij,wwp)
               dydxij = self.__matfpmullayer(nmax,n,dydxij,layer)
               for q in range(layer-1,-1,-1):
                  nq = self.nodes[q]
                  mq = self.nodes[q+1]
                  shiftq = self.wshift[q]
                  dydxij = self.__matmul(nmax,nq,mq,dydxij,w[shiftq:shiftq+mq*nq])
                  dydxij = self.__matfpmullayer(nmax,nq,dydxij,q)

               for rr in range(nlayers-1,layer,-1):
                  dydxrr  = [0.0]*nmax*nmax
                  if (rr==(nlayers-1)):
                     for k in range(nmax):
                        dydxrr[k+k*nmax] = self.dypdw[nlayers-1][k]
                  else:
                     for k in range(nmax):
                        dydxrr[k+k*nmax] = self.yp[nlayers-1][k]
                     for q in range(nlayers-2,rr,-1):
                        nq = self.nodes[q]
                        mq = self.nodes[q+1]
                        shiftq = self.wshift[q]
                        dydxrr = self.__matmul(nmax,nq,mq,dydxrr,w[shiftq:shiftq+mq*nq])
                        dydxrr = self.__matfpmullayer(nmax,nq,dydxrr,q)
                     nrr = self.nodes[rr]
                     mrr = self.nodes[rr+1]
                     shiftrr = self.wshift[rr]
                     dydxrr = self.__matmul(nmax,nrr,mrr,dydxrr,w[shiftrr:shiftrr+mrr*nrr])
                     for l in range(nrr):
                        for k in range(mrr):
                           dydxrr[k+l*mrr] *= self.dypdw[rr][l]
                  for q in range(rr-1,-1,-1):
                     nq = self.nodes[q]
                     mq = self.nodes[q+1]
                     shiftq = self.wshift[q]
                     dydxrr = self.__matmul(nmax,nq,mq,dydxrr,w[shiftq:shiftq+mq*nq])
                     dydxrr = self.__matfpmullayer(nmax,nq,dydxrr,q)
                  for k in range(len(dydxij)):
                    dydxij[k] += dydxrr[k]
               
               for k in range(len(dydxij)):
                  de = 2.0*(dyoutdx[k]-dytraindx[k])
                  dErrordw[i+j*m] += de*dydxij[k]

      return (Error,dErrordw)
      

#octave:16> a
#a =
#
#   0.455001   0.850642   0.835266   0.484402
#   0.619349   0.412941   0.582336   0.060414
#   0.200913   0.380119   0.644372   0.293560
#
#octave:17> b
#b =
#
#   0.712170   0.811203   0.904325   0.112777   0.214010
#   0.944679   0.458811   0.101055   0.764077   0.208326
#   0.336804   0.408634   0.989100   0.829954   0.218230
#   0.541491   0.851742   0.542150   0.099963   0.169209
#
#octave:18> c
#c =
#
#   1.67124   1.51329   1.58621   1.44292   0.53883
#   1.06003   0.98130   1.21056   0.87472   0.35588
#   0.87816   0.85073   1.01661   0.87724   0.31248

# a = [0.455001, 0.619349, 0.200913, 0.850642, 0.412941, 0.380119, 0.835266, 0.582336, 0.644372, 0.484402, 0.060414, 0.293560]
# b = [0.712170, 0.944679, 0.336804, 0.541491, 0.811203, 0.458811, 0.408634, 0.851742, 0.904325, 0.101055, 0.989100, 0.542150, 0.112777, 0.764077, 0.829954, 0.099963, 0.214010, 0.208326, 0.218230, 0.169209]
# c = [1.67124, 1.06003, 0.87816, 1.51329, 0.98130, 0.85073, 1.58621, 1.21056, 1.01661, 1.44292, 0.87472, 0.87724, 0.53883, 0.35588, 0.31248]

# v = [   0.533104, 0.067522, 0.422639, 0.228982]

