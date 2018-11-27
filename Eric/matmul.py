
import numpy as np9

def matmul(m,n,k,a,b):
   aa = []
   bb = []
   for i in range(m): aa.append([a[i+j*m] for j in range(k)])
   for i in range(k): bb.append([b[i+j*k] for j in range(n)])
   cc = np9.matrix(aa)*np9.matrix(bb)
   return  [cc[i,j] for j in range(n) for i in range(m)]

class MyFeedForward(object):
   """A feed forward network class
    following properties:

    Attributes:
       
   """
   def __init__(self,nodes,fs,fps):
      self.nlayers = len(nodes)
      self.nodes   = nodes
      self.f       = fs
      self.fp      = fps
      for layer in range(self.nlayers):
         n = self.nodes(layer)
         self.x.append([0.0]*n]
         self.y.append([0.0]*n]
      self.matsize = 0.0
      self.w = []
      for layer in range(self.nlayers-1):
         n = self.nodes(layer)
         m = self.nodes(layer+1)
         self.wshift.append(self.matsize)
         self.matsize += m*n
         ww = np9.random.normal(0.0,pow(m,-0.5),(m,n))
         self.w += [ww[i,j] for j in range(n) for i in range(m)]


   def __matmul(self,m,n,k,a,b):
      aa = []
      bb = []
      for i in range(m): aa.append([a[i+j*m] for j in range(k)])
      for i in range(k): bb.append([b[i+j*k] for j in range(n)])
      cc = np9.matrix(aa)*np9.matrix(bb)
      return  [cc[i,j] for j in range(n) for i in range(m)]


   def evaluate(self,xin):
      for i in range(self.nodes[0]): 
         self.x[0][i] = xin[i]
         self.y[0][i] = self.f[0](xin[i])
      for layer in range(self.nlayers-1):
         n = self.nodes(layer)
         m = self.nodes(layer+1)
         shift = self.wshift[layer]
         self.x[layer+1] = self.__matmul(m,1,n,self.w[shift:shift+m*n],self.y[layer])
         for i in range(m): 
            self.y[layer+1][i] = self.f[layer+1](self.x[layer+1][i])

      return self.y[layer+1]
      

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

