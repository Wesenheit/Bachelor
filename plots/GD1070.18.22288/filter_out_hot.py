import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyphot
import pystellibs
import extinction
from FAOV import Period
from scipy.optimize import curve_fit
from pyphot import (unit, Filter)
kpc=3.08567758*10**21
plt.rcParams["font.family"]="serif"
data_V=pd.read_csv('ASASSN-V J170801.81-410255.6.csv')
P=45.1467
data_OGLE=np.loadtxt("GD1070.18.22288.dat",dtype=np.float64)
time_OGLE=data_OGLE[:,0]# - 2400000.5 

data_g=pd.read_csv('light_curve_3407d752-9b3c-4ae0-8ae5-7da73ea2266b.csv')
data_g=data_g[data_g["mag"]<13.4]
data_g=data_g[data_g["mag_err"]<0.02]

T=np.log10(3882)
logL1=1.334
logL2=0.814
EBV=0.290
d=0.97
lib = pyphot.get_library()
def harmonic(x,*p):
    return (p[0]+p[1]*np.sin(2*np.pi*x/P)+p[2]*np.cos(2*np.pi*x/P)+
            p[3]*np.sin(2*2*np.pi*x/P)+p[4]*np.cos(2*2*np.pi*x/P)+
            p[5]*np.sin(3*2*np.pi*x/P)+p[6]*np.cos(4*2*np.pi*x/P)+
            p[7]*np.sin(3*2*np.pi*x/P)+p[8]*np.cos(4*2*np.pi*x/P))

library=pystellibs.BaSeL()
law=extinction.ccm89(np.array(library.wavelength),EBV*3.1,3.1)
#deredden_ogle=extinction.remove(law,)
ogleflux=10**(-0.4*data_OGLE[:,1])
ogleflux=ogleflux/np.mean(ogleflux)
sol=curve_fit(harmonic,time_OGLE,ogleflux,p0=np.array([1,0,0,0,0,0,0,0,0]))
params=sol[0]
log_red_lum_1=logL1+np.log10(harmonic(data_V["hjd"].values- 2400000.5 ,*params))
mag_list_V=[]
mag_list_G=[]
V_filter=lib["GROUND_JOHNSON_V"]
G_filter=lib["SDSS_g"]
for log_lum in log_red_lum_1:
    spec=library.generate_stellar_spectrum(T,2,log_lum,0.013)
    deredd=np.power(10,-0.4*law)*spec/(4*np.pi*d**2*kpc**2)
    mag_V = -2.5 * np.log10(V_filter.get_flux(np.array(library.wavelength)*unit["AA"], deredd*unit["flam"]).value) - V_filter.Vega_zero_mag
    mag_list_V.append(mag_V)


log_red_lum_2=logL1+np.log10(harmonic(data_g["HJD"]- 2400000.5 ,*params))

for log_lum in log_red_lum_2:
    mag_G = -2.5 * np.log10(G_filter.get_flux(np.array(library.wavelength)*unit["AA"], deredd*unit["flam"]).value) - G_filter.AB_zero_mag
    mag_list_G.append(mag_G)

mag_arr_V=np.array(mag_list_V)
mag_arr_G=np.array(mag_list_G)
flux_filter_V=10**(-0.4*data_V["mag"].values)-10**(-0.4*mag_arr_V)
flux_filter_G=10**(-0.4*data_g["mag"].values)-10**(-0.4*mag_arr_G)
mag_filtered_V=-2.5*np.log10(flux_filter_V)
mag_filtered_G=-2.5*np.log10(flux_filter_G)

new_period_G=Period(data_g["HJD"].values-2400000.5,mag_filtered_G,10,100)
new_period_V=Period(data_V["hjd"].values-2400000.5,mag_filtered_V,10,100)
fig,axes=plt.subplots(1,2,figsize=(12,8))
ax1=axes.flatten()[0]
ax1.errorbar(np.remainder(data_V["hjd"].values-2400000.5,2*new_period_V)/new_period_V,data_V["mag"],data_V["mag_err"],fmt="g.",label="observed")
ax1.scatter(np.remainder(data_V["hjd"].values-2400000.5,2*new_period_V)/new_period_V,mag_filtered_V,marker=".",label="prediced",color="black")
ax1.grid()
ax1.legend()
ax1.set_xlabel("Phase")
ax1.set_ylabel("V [mag]")
ax2=axes.flatten()[1]
ax2.errorbar(np.remainder(data_g["HJD"].values-2400000.5,2*new_period_G)/new_period_G,data_g["mag"],data_g["mag_err"],fmt="g.",label="observed")
ax2.scatter(np.remainder(data_g["HJD"].values-2400000.5,2*new_period_G)/new_period_G,mag_filtered_G,marker=".",label="prediced",color="black")
ax2.grid()
ax2.legend()
ax2.set_xlabel("Phase")
ax2.set_ylabel("SDSS g [mag]")
ax1.invert_yaxis()
ax2.invert_yaxis()
ax1.set_title("V, $P={:.3f}$".format(new_period_V))
ax2.set_title("SDSS g, $P={:.3f}$".format(new_period_G))
plt.savefig("filtr.png",dpi=1000)

