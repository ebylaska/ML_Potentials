import myfeedforward5 as myfeedforward
import math,random

beta = 1.5
sigmoid      = lambda x: math.tanh(beta*x)
sigmoidp     = lambda x: beta*(1.0/math.cosh(beta*x))**2
sigmoidpp    = lambda x: -2.0*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
xmoid = lambda x: x*x
xmoidp = lambda x: 2*x
xmoidpp = lambda x: 2
machine = myfeedforward.MyFeedForward([1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp])
weights = machine.initial_w()


alpha = 0.05
A = 0.1
for i in range(1000000):
   x = 2*random.random()-1
   e = A*x*x
   e1 = machine.evaluate([x],weights)[0]
   gg = machine.w_energy_gradient([x],[e],weights)
   error = gg[0]
   g1 = gg[1]
   for j in range(len(g1)):
      weights[j] -= alpha*g1[j]
   #print "x,e,e1,error,w=",x,e,e1,error,weights
   print x,e,e1,(e1-e)**2,error,weights[0],weights[1],weights[2],weights[3]


