'''
To calculate the damage energy (ED) using the Lindhard method for a given recoil energy of a primary knock-on atom (PKA).
ref: 
1. Gary Was, Fundamentals of Radiation Materials Science, p90-91
2. M.J. Norgett et al, Method of calculatitig displacement dose rates
3. M.T. Robinson, THE ENERGY DEPENDENCE OF NEUTRON RADIATION DAMAGE IN SOLIDS
'''
import math
import matplotlib.pyplot as plt
# Parameters:
# T: PKA energy (recoil energy)
# A1, Z1: mass number, atomic number of the projectile
# A2, Z2: ~ of the target
def ED(T,A1,Z1,A2,Z2):
    return T/(1+k(A1,Z1)*g(e(A1,Z1,A2,Z2,T)))
def g(e):
    return 3.4008*pow(e,1/6)+0.40244*pow(e,3/4)+e 
def k(A1,Z1):
    return 0.1337*pow(Z1,1/6)*pow(Z1/A1,1/2)
def e(A1,Z1,A2,Z2,T):
    return (A2*T)/(A1+A2)*a(Z1,Z2)/(Z1*Z2*pow(1.60e-19,1))*(4*3.14*8.85e-12)
def a(Z1,Z2):
    return pow(9*3.14*3.14/128,1/3)*0.529e-10*pow(pow(Z1,2/3)+pow(Z2,2/3),-1/2)

T_list = [ pow(10,i/10) for i in range(10,61,1)]
# For identical atoms:
# 1. Be 
A1=A2=9.01
Z1=Z2=4
Deffi_Be_list = []
for T in T_list:
    Ed = ED(T,A1,Z1,A2,Z2)
    Deffi_Be_list.append(Ed/T)
plt.semilogx(T_list, Deffi_Be_list, color="green", label="Be")

#2. U 
A1=A2=238
Z1=Z2=92
Deffi_U_list = []
for T in T_list:
    Ed = ED(T,A1,Z1,A2,Z2)
    Deffi_U_list.append(Ed/T)
plt.semilogx(T_list, Deffi_U_list, color="red", label="U")

#3. W
A1=A2=183.84
Z1=Z2=74
Deffi_W_list = []
for T in T_list:
    Ed = ED(T,A1,Z1,A2,Z2)
    Deffi_W_list.append(Ed/T)
plt.semilogx(T_list, Deffi_W_list, color=(0.1,0.2,0.5), label="W")

#For different atoms:
#4. projectile: He; target: W
A1=4
Z1=2
A2=183.84
Z2=74
Deffi_HeW_list = []
for T in T_list:
    Ed = ED(T,A1,Z1,A2,Z2)
    Deffi_HeW_list.append(Ed/T)
plt.semilogx(T_list, Deffi_HeW_list, color=(0.8,0.5,0.1), label="HeW")

plt.legend()
plt.xlabel("Recoil energy (eV)")
plt.ylabel("Damage efficiency")
plt.xlim((10,1000000))
plt.ylim((0,1))
plt.show()

