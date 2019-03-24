import myfeedforward5 as myfeedforward
import math,random

beta = 1.0
sigmoid      = lambda x: math.tanh(beta*x)
sigmoidp     = lambda x: beta*(1.0/math.cosh(beta*x))**2
sigmoidpp    = lambda x: -2.0*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
xmoid = lambda x: x
xmoidp = lambda x: 1
xmoidpp = lambda x: 0
machine = myfeedforward.MyFeedForward([1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp])
weights = machine.initial_w()


d = 3
alpha = 0.05
A = 0.1
dx = 0.02
dw = 0.01
x  = -0.1
xp = x+dx
xm = x-dx
e = A*x*x

gg = machine.w_energy_gradient([x],[e],weights)

weightsp = weights[:]
weightsm = weights[:]
weightsp[d] += dw
weightsm[d] -= dw
errorp = (machine.evaluate([x],weightsp)[0]-e)**2
errorm = (machine.evaluate([x],weightsm)[0]-e)**2

gd= (errorp-errorm)/(2.0*dw)

print "gg=",gg
print "d,gd,gg[1][d]=",d,gd,gg[1][d]



