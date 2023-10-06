import numpy as np
## From Compressor Design
dDelta3 = 0.1697689384124046 #cm copy from compressor outlet boundary layer thickness
rh3 = 1.25 # cm hub radius cm copy from compressor 
RPM = 60000 #Copy from compressor
## From cycle
AFlow1 = 24.185 # cm2
AFlow3 = 35.389 # cm2
mdot = 2.70E-01 #kg/sec
Mx1 = 2.99E-01
Mx3 = 2.90E-01
Vx1 = 165.742 # m/sec
Vx3 = 155.510 # m/sec
Pt1 = 1.63E+05 #Pa
Pt3 = 1.11E+05 #Pa
Tt1 = 802.588 #K
Cp1 = 1125.828 #J/(kg*K)
Cp3 = 1112.555 #J/(kg*K)
Gam1 = 1.344
Gam3 = 1.350
Rconst = 288.469 #J/(kg*K)
ht2 = 8.39E+05 # J/kg
ht3 = 7.81E+05 # J/kg

####Design parameters
NB3 = 10 #number of turbine blades
hB3 = 0.2 #cm Turbine blade thickness
NB1 = 24
hB1 = 0.2 #cm nozzle blade thickness
preOverPredictFrac = 0.05 #5% over prediction
rm2 = 4.52 #cm Pick rm2 radius
M2Abs = 1.1 #Assume outlet of nozzle is choked (absolute mach number is unity)
####Design parameters

#Compute station 3 radius assume boundary layer thickness from compressor outlet (No outlet swirl)
a = np.pi
b = -2*np.pi*dDelta3 -(2*dDelta3+hB3)*NB3
c = -(AFlow3 + np.pi*rh3**2 + 2*np.pi*rh3*dDelta3 - rh3*(2*dDelta3+hB3)*NB3)
rt3 = (-b+np.sqrt(b**2 - 4*a*c))/(2*a) #cm
print('Station 3 tip radius is :' + str(rt3) + " cm")
omega = RPM*(2*np.pi/60) # rad/sec
beta3h = np.arctan((-omega*rh3*0.01)/Vx3)*(180/np.pi) # deg
beta3t = np.arctan((-omega*rt3*0.01)/Vx3)*(180/np.pi) # deg
print("Station 3 beta angles are (Hub and Tip) :" + str(beta3h)+"deg & "+str(beta3t)+" deg")

#Compute tangential velocity from assumed radius (Assumed no outlet swirl at design point)
dDelta2 = dDelta3/2 #Assume station 2 boundary layer thickness is half of that of station 3
Gam2 = 0.5*(Gam1+Gam3)
Cp2 = 0.5*(Cp1+Cp3) #J/(kg*K)
Ut2 = omega*(rm2*0.01) #m/sec blade tip speed at station 2
ht3OP = ht3*(1-preOverPredictFrac) #make turbine over expand to avoid risk
Vtheta2 = (ht2-ht3OP)/Ut2 #m/sec

#Compute station 2 nozzle angle
T2 = Tt1*(1+0.5*(Gam2-1)*(M2Abs**2))**(-1) #K
P2 = Pt1*((T2/Tt1)**(Gam2/(Gam2-1))) #Pa
SOS2 = np.sqrt(Gam2*Rconst*T2) #speed of sound m/sec at station 2
Mtheta2Abs = Vtheta2/SOS2 #Tangential absolute mach number
Mx2Abs = np.sqrt(M2Abs**2 - Mtheta2Abs**2) #absolute station 2 axial mach number
Vx2 = Mx2Abs*SOS2 #m/sec station 2 absolute axial velocity
alpha2 = np.arctan(Vtheta2/Vx2)*(180/np.pi) # deg
print("Alpha 2 angle is : "+str(alpha2)+" deg")
print("Absolute station 2 mach number is " + str(np.sqrt(Mx2Abs**2 + Mtheta2Abs**2)))
Wtheta2 = Vtheta2 - Ut2 #m/sec
beta2 = np.arctan(Wtheta2/Vx2)*(180/np.pi) # deg
W2 = np.sqrt(Wtheta2**2 + Vx2**2) #m/sec
M2Rel = W2/SOS2
print("Beta 2 angle is : "+str(beta2) + " deg")
print("Relative station 2 mach number is "+ str(M2Rel))

#Compute station 2 and 1 flow area
rho2 = P2/(Rconst*T2) #kg/m3
AFlow2 = (mdot/rho2/Vx2)*10000 #cm2
Hc2 = (AFlow2 + 2*np.pi*rm2*(2*dDelta2))/(2*np.pi*rm2 - NB1*(2*dDelta2 + hB1)) #cm
print("Station 2 Channel height is "+str(Hc2)+" cm")
rm1 = (AFlow1 + Hc2*(2*dDelta2+hB1)*NB1)/(2*np.pi*Hc2 - 2*np.pi*(2*dDelta2)) #cm
print('Station 1 radius is :' + str(rm1) + " cm")

#Compute reaction of the turbine
P1 = Pt1*(1+0.5*(Gam1-1)*Mx1**2)**(Gam1/(1-Gam1))
P3 = Pt3*(1+0.5*(Gam3-1)*Mx3**2)**(Gam3/(1-Gam3))
print("Station 1 static pressure is:" + str(P1) + " Pa")
print("Station 2 static pressure is:" + str(P2) + " Pa")
print("Station 3 static pressure is:" + str(P3) + " Pa")
Reaction = (P2-P3)/(P1-P3)
print("Reaction of the stage is: " + str(Reaction))

#Compute blade speed ratio
V2 = np.sqrt(Vx2**2 + Vtheta2**2)
BSR = Ut2/V2
print("Blade speed ratio is "+str(BSR))
