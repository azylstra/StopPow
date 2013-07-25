from StopPow import *

models = [] # variety of stop pow models

# SRIM:
# solid aluminum
models.append( StopPow_SRIM("data/Hydrogen in Aluminum.txt") )

# Li-Petrasso:
# H plasma at 1e23 ions/cc and 1keV temperature
mf = FloatVector(2)
mf[0] = 1
mf[1] = 1/1800
Zf = FloatVector(2)
Zf[0] = 1
Zf[1] = -1
Tf = FloatVector(2)
Tf[0] = 1
Tf[1] = 1
nf = FloatVector(2)
nf[0] = 1e23
nf[1] = 1e23
models.append( StopPow_LP(1,1,mf,Zf,Tf,nf) )

# Bethe-Bloch
# Solid density C (~diamond)
mf2 = FloatVector(1)
mf2[0] = 12
Zf2 = FloatVector(1)
Zf2[0] = 6
nf2 = FloatVector(1)
nf2[0] = 1.76e23
models.append( StopPow_BetheBloch(1,1,mf2,Zf2,nf2) )

# perform variety of tests for each model:
for x in models:
	x.set_mode(x.MODE_LENGTH)
	print("dEdx(10 MeV) = " , x.dEdx(10.0) , " MeV/um")
	x.set_mode(x.MODE_RHOR)
	print("dEdx(10 MeV) = " , x.dEdx(10.0) , " MeV/(mg/cm2)")


	x.set_mode(x.MODE_LENGTH)
	print("Eout(10 MeV, 100um) = " , x.Eout(10,100) )
	print("Ein(10 MeV, 100um) = " , x.Ein(10,100) )
	print("Thickness(10 MeV, 9 MeV) = " , x.Thickness(10,9) )
	print("----------------")

from datetime import *
t0 = datetime.now()
for i in range(1000):
	s = StopPow_LP(1,1,mf,Zf,Tf,nf)
t1 = datetime.now()
print('{:.1f}'.format((t1 - t0).total_seconds()*1e3) + "us per L-P creation")
from datetime import *
t0 = datetime.now()
for i in range(1000):
	models[1].dEdx(10.)
t1 = datetime.now()
print('{:.1f}'.format((t1 - t0).total_seconds()*1e3) + "us per L-P dE/dx call")
from datetime import *
t0 = datetime.now()
for i in range(1000):
	models[1].Eout(10.,100)
t1 = datetime.now()
print('{:.1f}'.format((t1 - t0).total_seconds()*1e3) + "us per L-P Eout call")