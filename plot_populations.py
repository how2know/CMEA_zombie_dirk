from numpy import *
from pylab import *

S = loadtxt("S.txt")
Z = loadtxt("Z.txt")
R = loadtxt("R.txt")
t = loadtxt('time.txt')
plot(t,S,label="S")
plot(t,Z,label="Z")
plot(t,R,label="R")

xlabel(r'$Time$')
ylabel(r'$Individuals$')

legend(loc = "best")

#savefig("../../figs/SZR_fe_" + str(t.size-1) + ".png")
#savefig("../../figs/SZR_dirk_" + str(t.size-1) + ".png")

show()

