
import math
import numpy as np9

#def matmul(m,n,k,a,b):
#   aa = []
#   bb = []
#   for i in range(m): aa.append([a[i+j*m] for j in range(k)])
#   for i in range(k): bb.append([b[i+j*k] for j in range(n)])
#   cc = np9.matrix(aa)*np9.matrix(bb)
#   return  [cc[i,j] for j in range(n) for i in range(m)]

class MyFeedForward(object):
   """A feed forward network class
    following properties:

    Attributes:
       
   """
   def __init__(self,nodes,fs,fps,fpps):
      self.nlayers = len(nodes)
      self.nodes   = nodes
      self.f       = fs
      self.fp      = fps
      self.fpp     = fpps
      self.x = []
      self.y = []
      self.yp = []
      self.ypp = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         self.x.append([0.0]*n)
         self.y.append([0.0]*n)
         self.yp.append([0.0]*n)
         self.ypp.append([0.0]*n)
            
      self.matsize = 0
      self.bshift = []
      self.wshift = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         self.bshift.append(self.matsize)
         self.matsize += n
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
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         w += [0.0000]*n
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         ww = np9.random.normal(0.0,pow(m,-0.5),(m,n))
         w += [ww[i,j] for j in range(n) for i in range(m)]
      return w

   def print_w(self,w):
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         bshift = self.bshift[layer]
         print("b layer,n=",layer,n)
         str = ""
         for i in range(n):
            str += "%f " % w[bshift+i]
         str += "\n"
         print(str)
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         print("w layer,n,m=",layer,n,m)
         str = ""
         for i in range(m):
            for j in range(n):
               str += "%f " % w[shift+i+j*m]
            str += "\n"
         print(str)

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
         #self.y[0][i] = self.f[0](xin[i]+self.biases[0][i])
         self.y[0][i] = self.f[0](xin[i]+w[i])
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift  = self.wshift[layer]
         bshift = self.bshift[layer+1]
         self.x[layer+1] = self.__matmul(m,1,n,w[shift:shift+m*n],self.y[layer])
         for i in range(m): 
            self.y[layer+1][i] = self.f[layer+1](self.x[layer+1][i]+w[bshift+i])

      return self.y[self.nlayers-1][:]

   #
   # evaluate the gradient wrt x, dyout/dx
   #
   def gradients_evaluate(self,xin,w):
      #...initialize values for x
      yout = self.evaluate(xin,w)

      # define yp
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         bshift = self.bshift[layer]
         for k in range(n):
            self.yp[layer][k] = self.fp[layer](self.x[layer][k]+w[bshift+k])

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
         bshift = self.bshift[layer]
         tmp0 = de*self.fp[layer](self.x[layer][i]+w[bshift+i])
         dErrordw[bshift+i] += tmp0
         for layer in range(self.nlayers-2,-1,-1):
            n = self.nodes[layer]
            m = self.nodes[layer+1]
            shift  = self.wshift[layer]
            bshift = self.bshift[layer]
            ww = w[shift:shift+m*n]
            for k in range(n):
               for j in range(m):
                  ww[j+k*m] *=  self.fp[layer](self.x[layer][k]+w[bshift+k])

            if (layer==self.nlayers-2):
               for k in range(n):
                  dErrordw[shift+i+k*m] += tmp0*self.y[layer][k]
               tmp1 = []
               for j in range(n):
                  xx = tmp0*ww[i+j*m]
                  tmp1.append(xx)
                  dErrordw[bshift+j] += xx
            else:
               for k in range(n):
                  for j in range(m):
                     dErrordw[shift+j+k*m] += tmp1[j]*self.y[layer][k]
               tmp1 = self.__matmul(1,n,m,tmp1,ww)
               #print("size tmp1=",len(tmp1),layer,n,m)
               for j in range(n):
                  xx = tmp1[j]
                  dErrordw[bshift+j] += xx
            
      return (Error,dErrordw)

   #
   # calculate the error and gradient of w using derivative training
   #
   def w_energy_gradient1(self,xtrain,dytraindx,w):
      Error     = 0.0
      dErrordw  = [0.0]*self.matsize
      dyoutdx = self.gradients_evaluate(xtrain,w)
      for i in range(len(dytraindx)):
         Error += (dyoutdx[i]-dytraindx[i])**2

      ### define ypp and dxdw template
      dxdw = []
      dypdw = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         bshift = self.bshift[layer]
         dxdw.append([0.0]*n)
         dypdw.append([0.0]*n)
         for k in range(n):
            self.ypp[layer][k] =  self.fpp[layer](self.x[layer][k]+w[bshift+k])

      #print("dyp2/dw1=",self.ypp[2][0]*self.y[1][0])
      #print("dyp1/dw1=",0.0)
      #print("dyp2/dw0=",self.ypp[2][0]*w[1]*self.yp[1][0]*self.y[0][0])
      #print("dyp1/dw0=",self.ypp[1][0]*self.y[0][0])

      #print("dError/dw0=",2.0*(dyoutdx[0]-dytraindx[0])*(self.ypp[2][0]*w[1]*self.yp[1][0]*self.y[0][0]*w[1]*self.yp[1][0]*w[0]*self.yp[0][0] + self.yp[2][0]*w[1]*self.ypp[1][0]*self.y[0][0]*w[0]*self.yp[0][0]+self.yp[2][0]*w[1]*self.yp[1][0]*self.yp[0][0]))
      #print("dError/dw1=",2.0*(dyoutdx[0]-dytraindx[0])*(self.ypp[2][0]*self.y[1][0]*w[1]*self.yp[1][0]*w[0]*self.yp[0][0] + self.yp[2][0]*self.yp[1][0]*w[0]*self.yp[0][0]))

      #print("dg2/dw0=",self.ypp[2][0]*w[1]*self.yp[1][0]*self.y[0][0]*w[1]*self.yp[1][0]*w[0]*self.yp[0][0] )
      #print("dg1/dw0=",self.yp[2][0]*w[1]*self.ypp[1][0]*self.y[0][0]*w[0]*self.yp[0][0])
      #print("dg0/dw0=",self.yp[2][0]*w[1]*self.yp[1][0]*self.yp[0][0])

      #print("dg1/dw1=",self.ypp[2][0]*w[1]*self.yp[1][0]*self.y[0][0]*w[1]*self.yp[1][0]*w[0]*self.yp[0][0])
      #print("dg0/dw1=",self.yp[2][0]*self.yp[1][0]*w[0]*self.yp[0][0])
      #prin()
      
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         for j in range(n):
            for i in range(m):
               # dw/dwij
               wwp = [0.0]*m*n
               wwp[i+j*m] = 1.0
               
               for q in range(self.nlayers):
                  dxdw[q]  = [0.0]*n
                  dypdw[q] = [0.0]*n

               # define needed dxdw's and dfp/dw's
               dxdw[layer+1] = self.__matmul(m,1,n,wwp,self.y[layer])
               for k in range(m): dypdw[layer+1][k] = self.ypp[layer+1][k]*dxdw[layer+1][k]
               for q in range(layer+2,self.nlayers):
                  nqm = self.nodes[q-1]
                  mqm = self.nodes[q]
                  shiftqm = self.wshift[q-1]
                  wwq = w[shiftqm:shiftqm+mqm*nqm]
                  tmp     = self.__matfpmullayer(mqm,nqm,wwq,q-1)
                  dxdw[q] = self.__matmul(mqm,1,nqm,tmp,dxdw[q-1])
                  for k in range(mqm): dypdw[q][k] = self.ypp[q][k]*dxdw[q][k]
                
               # calculate dydxij
               nlayers = self.nlayers
               nmax    = self.nodes[nlayers-1]
               dydxij  = [0.0]*nmax*nmax
               for k in range(nmax):
                  dydxij[k+k*nmax] = self.yp[nlayers-1][k]
               for q in range(nlayers-2,layer,-1):
                  nq = self.nodes[q]
                  mq = self.nodes[q+1]
                  shiftq = self.wshift[q]
                  dydxij = self.__matmul(nmax,nq,mq,dydxij,w[shiftq:shiftq+mq*nq])
                  dydxij = self.__matfpmullayer(nmax,nq,dydxij,q)
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
                        dydxrr[k+k*nmax] = dypdw[nlayers-1][k]
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
                        for k in range(nmax):
                           dydxrr[k+l*nmax] *= dypdw[rr][l]
                  
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
                  dErrordw[shift+i+j*m] += de*dydxij[k]
               #print("layer,shift,i,j,shift+i+j*m,dErrordw=",layer,shift,i,j,shift+i+j*m,dErrordw)
               

      return (Error,dErrordw)
      
   #
   # calculate the gradient of yout wrt w 
   #
   def dyoutdw_gradient(self,xtrain,w):
      dyoutdw  = [0.0]*self.matsize*self.nodes[-1]
      yout = self.evaluate(xtrain,w)
      for i in range(self.nodes[-1]):
         shifti = self.matsize*i
         de =  1.0
         layer = self.nlayers - 1
         bshift = self.bshift[layer]
         tmp0 = de*self.fp[layer](self.x[layer][i]+w[bshift+i])
         for layer in range(self.nlayers-2,-1,-1):
            n = self.nodes[layer]
            m = self.nodes[layer+1]
            shift = self.wshift[layer]
            bshift = self.bshift[layer]
            ww = w[shift:shift+m*n]
            for k in range(n):
               for j in range(m):
                  ww[j+k*m] *=  self.fp[layer](self.x[layer][k]+w[bshift+k])

            if (layer==self.nlayers-2):
               for k in range(n):
                  dyoutdw[shifti+shift+i+k*m] += tmp0*self.y[layer][k]
               tmp1 = []
               for j in range(n):
                  xx = tmp0*ww[i+j*m]
                  tmp1.append(xx)
            else:
               for k in range(n):
                  for j in range(m):
                     dyoutdw[shifti+shift+j+k*m] += tmp1[j]*self.y[layer][k]
               tmp1 = self.__matmul(1,n,m,tmp1,ww)
            
      return (yout,dyoutdw)


   #
   # calculate the gradient dyout/dx and the second derivative  of ddyout/dxdw 
   #
   def ddyoutdxdw_gradient(self,xtrain,w):
   
      print("msize=",self.matsize*self.nodes[0]*self.nodes[-1],self.nlayers-1)
      ddyoutdxdw  = [0.0]*self.matsize*self.nodes[0]*self.nodes[-1]
      dyoutdx     = self.gradients_evaluate(xtrain,w)

      ### define ypp and dxdw template
      dxdw = []
      dypdw = []
      for layer in range(self.nlayers):
         n = self.nodes[layer]
         bshift = self.bshift[layer]
         dxdw.append([0.0]*n)
         dypdw.append([0.0]*n)
         for k in range(n):
            self.ypp[layer][k] =  self.fpp[layer](self.x[layer][k]+w[bshift+k])

      
      for layer in range(self.nlayers-1):
         n = self.nodes[layer]
         m = self.nodes[layer+1]
         shift = self.wshift[layer]
         for j in range(n):
            for i in range(m):
               # dw/dwij
               wwp = [0.0]*m*n
               wwp[i+j*m] = 1.0
               
               for q in range(self.nlayers):
                  dxdw[q]  = [0.0]*n
                  dypdw[q] = [0.0]*n

               # define needed dxdw's and dfp/dw's
               dxdw[layer+1] = self.__matmul(m,1,n,wwp,self.y[layer])
               for k in range(m): dypdw[layer+1][k] = self.ypp[layer+1][k]*dxdw[layer+1][k]
               for q in range(layer+2,self.nlayers):
                  nqm = self.nodes[q-1]
                  mqm = self.nodes[q]
                  shiftqm = self.wshift[q-1]
                  wwq = w[shiftqm:shiftqm+mqm*nqm]
                  tmp     = self.__matfpmullayer(mqm,nqm,wwq,q-1)
                  dxdw[q] = self.__matmul(mqm,1,nqm,tmp,dxdw[q-1])
                  for k in range(mqm): dypdw[q][k] = self.ypp[q][k]*dxdw[q][k]
                
               # calculate dydxij
               nlayers = self.nlayers
               nmax    = self.nodes[nlayers-1]
               dydxij  = [0.0]*nmax*nmax
               for k in range(nmax):
                  dydxij[k+k*nmax] = self.yp[nlayers-1][k]
               for q in range(nlayers-2,layer,-1):
                  nq = self.nodes[q]
                  mq = self.nodes[q+1]
                  shiftq = self.wshift[q]
                  dydxij = self.__matmul(nmax,nq,mq,dydxij,w[shiftq:shiftq+mq*nq])
                  dydxij = self.__matfpmullayer(nmax,nq,dydxij,q)
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
                        dydxrr[k+k*nmax] = dypdw[nlayers-1][k]
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
                        for k in range(nmax):
                           dydxrr[k+l*nmax] *= dypdw[rr][l]
                  
                  for q in range(rr-1,-1,-1):
                     nq = self.nodes[q]
                     mq = self.nodes[q+1]
                     shiftq = self.wshift[q]
                     dydxrr = self.__matmul(nmax,nq,mq,dydxrr,w[shiftq:shiftq+mq*nq])
                     dydxrr = self.__matfpmullayer(nmax,nq,dydxrr,q)
                  for k in range(len(dydxij)):
                    dydxij[k] += dydxrr[k]
                  
               
               for k in range(len(dydxij)):
                  shiftk = k*self.matsize
                  de     = 1.0
                  ddyoutdxdw[shiftk+shift+i+j*m] += de*dydxij[k]
               

      return (dyoutdx,ddyoutdxdw)
      


def main():
#
   alpha0 = 0.01
   alpha = 0.0001
   beta1 = 0.9
   beta2 = 0.999
   eps   = 1e-8

   beta = 2.0
   #sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
   #sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
   #sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
   ap = 1.0
   xp = 4.5
   bp = 3.0
   penalty  = lambda x: ap*(0.5*(math.tanh(bp*(x-xp)) - math.tanh(bp*(x+xp))) + 1.0)
   penaltyp = lambda x: ap*0.5*bp*( (1/math.cosh(bp*(x-xp)))**2 - (1.0/math.cosh(bp*(x+xp)))**2)

   sigmoid      = lambda x: 0.5*(math.tanh(beta*x)+1.0)
   sigmoidp     = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
   sigmoidpp    = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
   xmoid1   = lambda x: x
   xmoidp1  = lambda x: 1.0
   xmoidpp1 = lambda x: 0.0

   #bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
   machine = MyFeedForward([3,2,3,2],[xmoid1,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,xmoidpp1])

   weights = machine.initial_w()
   nw = len(weights)

   xs = -0.567
   ys =  0.67
   zs =  -0.17
   es = 0.73839
   ws = 0.03839
   gg = machine.w_energy_gradient([xs,ys,zs],[es,ws],weights)
   error = gg[0]
   g1    = gg[1]
   print("error=",error)

   es1 = machine.evaluate([xs,ys,zs],weights)
   print("es1=",es1)
   print("xs,ys,zs,es,es1=",xs,ys,zs,es,ws,es1,(es-es1[0])**2 + (ws-es1[1])**2)
   print("len(g1)=",len(g1), 7+2+6+3)
   print("gb=",g1[:7])
   print("gw=",g1[7:])

   gg0 = machine.w_energy_gradient([xs,ys,zs],[es,ws],weights)

   delta = 0.00001
   for ii in range(len(g1)):
      weights[ii] += delta
      gg1 = machine.w_energy_gradient([xs,ys,zs],[es,ws],weights)
      weights[ii] -= 2*delta
      gg2 = machine.w_energy_gradient([xs,ys,zs],[es,ws],weights)

      print("gradient=",ii,weights[ii],(gg1[0]-gg2[0])/(2*delta),gg0[1][ii])

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


if __name__ == "__main__":
   main()
