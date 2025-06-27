from __future__ import division, print_function, absolute_import
from math import pi,log,exp
from scipy.optimize import brentq #solver to find roots (zero points) of functions
from scipy.interpolate import interp1d
import numpy as np

import CoolProp as CP
from ACHP.MicroChannelEvaporator import MicroChannelEvaporatorClass

from ACHP.Correlations import f_h_1phase_MicroTube,ShahEvaporation_Average,LMPressureGradientAvg,AccelPressureDrop,TwoPhaseDensity,Bertsch_MC_Average,KM_Evap_Average,KandlikarEvaporation_average
from ACHP.MicroFinCorrelations import MultiLouveredMicroFins, MicroFinInputs, IsFinsClass
from ACHP.DryWetSegment import DWSVals, DryWetSegment
from ACHP.ACHPTools import ValidateFields

#Example usage for a parametric study
from CoolProp.CoolProp import PropsSI
# import pylab
import matplotlib.pyplot as pylab

num_points= 51
T_dews= np.linspace(275,292,num_points)
TT= np.empty(num_points)
Q_2p= np.empty(num_points)
w_2p= np.empty(num_points)
w_sh= np.empty(num_points)
Q_tot= np.empty(num_points)
h_2p= np.empty(num_points)
h_sh= np.empty(num_points)

FinsTubes=MicroFinInputs()

FinsTubes.Tubes.NTubes=61.354           #Number of tubes (per bank for now!)
FinsTubes.Tubes.Nbank=1                 #Number of banks (set to 1 for now!)
FinsTubes.Tubes.Npass=3                 #Number of passes (per bank-averaged)
FinsTubes.Tubes.Nports=1                #Number of rectangular ports
FinsTubes.Tubes.Ltube=0.30213           #length of a single tube
FinsTubes.Tubes.Td=0.0333               #Tube outside width (depth)
FinsTubes.Tubes.Ht= 0.002               #Tube outside height (major diameter)
FinsTubes.Tubes.b=0.00635               #Tube spacing
FinsTubes.Tubes.tw=0.0003               #Tube wall thickness
FinsTubes.Tubes.twp=0.0003              #Port (channel) wall thickness
FinsTubes.Tubes.beta=1                  #Port (channel) aspect ratio (=width/height)
FinsTubes.Tubes.kw=117                  #wall thermal conductivity

FinsTubes.Fins.FPI=11.0998              #Fin per inch
FinsTubes.Fins.Lf=0.0333                #Fin length
FinsTubes.Fins.t=0.000152               #Fin thickness
FinsTubes.Fins.k_fin=117                #Fin thermal conductivity

FinsTubes.Air.Vdot_ha=0.5663            #Air volume flow rate in m^3/s
FinsTubes.Air.Tmean=299.9
FinsTubes.Air.Tdb=299.9                   #Air inlet temperature, K
FinsTubes.Air.p=101325                  #Air pressure in Pa
FinsTubes.Air.RHmean=0.51
FinsTubes.Air.RH=0.51                   #Air inlet relative humidity
FinsTubes.Air.FanPower=438              #Fan power, Watts

FinsTubes.Louvers.Lalpha=20             #Louver angle, in degree
FinsTubes.Louvers.lp=0.001              #Louver pitch
FinsTubes.Louvers.Llouv=0.005737        #Louver cut length

#Abstract State
Ref = 'R744'
Backend = 'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
AS = CP.AbstractState(Backend, Ref)

kwargs={'AS': AS,
        'Ref':'R744',
        'mdot_r':  0.0708,
        'psat_r':  PropsSI('P','T',T_dews[0],'Q',1.0,Ref),
        'Fins': FinsTubes,
        'FinsType': 'MultiLouveredMicroFins',
        'hin_r': PropsSI('H','P',PropsSI('P','T',282,'Q',1.0,Ref),'Q',0.15,Ref),
        'Verbosity': 0,
        'h_a_tuning':1,
        'h_tp_tuning':1,
        'DP_tuning':1,
    }

MicroEvap=MicroChannelEvaporatorClass(**kwargs) #generate new micro-channel evaporator instance and update kwargs

for i in range(0, len(T_dews)):
    kwargs={'psat_r':  PropsSI('P','T',T_dews[i],'Q',1.0,Ref)}
    MicroEvap.Update(**kwargs)
    MicroEvap.Calculate()
    Q_tot[i] = MicroEvap.Q
    Q_2p[i]= MicroEvap.Q_2phase
    w_2p[i]= MicroEvap.w_2phase
    w_sh[i]= MicroEvap.w_superheat
    h_2p[i]= MicroEvap.h_r_2phase
    h_sh[i]= MicroEvap.h_r_superheat

print ("Demonstrate output list")
#print (Evap.OutputList())
for id, unit, value in MicroEvap.OutputList():
    print (str(id) + ' = ' + str(value) + ' ' + str(unit))

pylab.plot(T_dews,Q_2p,T_dews,Q_tot)
pylab.title('Parametric Study With Fixed flowrates - Capacity')
pylab.legend(['two-phase','total'],loc='best')
pylab.title('Parametric Study With Fixed flowrates - Capacity')
pylab.xlabel('Evaporation Dew Temperature in Kelvin')
pylab.ylabel('Capacity in Watt')
#pylab.savefig('Evaporator_py_capacity.pdf')
pylab.show()
pylab.plot(T_dews,h_2p,T_dews, h_sh)
pylab.title('Parametric Study with fixed Flowrates - Heat Transfer Coefficients')
pylab.legend(['two-phase','superheat'],loc='best')
pylab.xlabel('Evaporation Dew Temperature in Kelvin')
pylab.ylabel('Heat Transfer Coefficient in W/m2-K')
#pylab.savefig('Evaporator_py_HTC.pdf')
pylab.show()
pylab.plot(T_dews,w_2p, T_dews, w_sh)
pylab.title('Parametric Study with fixed Flowrates - Area Fraction')
pylab.legend(['two-phase', 'superheat'],loc='best')
pylab.xlabel('Evaporation Dew Temperature in Kelvin')
pylab.ylabel('Two-phase Wetted Area Fraction')
pylab.ylim(-0.01,1.01)
#pylab.savefig('Evaporator_py_wetted_area.pdf')
pylab.show()