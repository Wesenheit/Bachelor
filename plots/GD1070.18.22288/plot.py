import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyphot
import pystellibs
import extinction
plt.rcParams["font.family"]="serif"

fig,axes=plt.subplots(2,2,figsize=(12,8))
P=45.146726862302486
bochium=np.loadtxt("bochium.csv",skiprows=1)
ax1=axes.flatten()[0]
T_boch=[55686.40097, 55687.27062, 56411.253669, 56412.268657, 56553.009201, 56555.009595, 56557.042616, 56796.35103, 56805.349178, 56806.355972, 56808.377662, 56809.375058, 56825.093866, 56828.093611, 56856.146771, 56857.149595, 56863.156238, 56867.226354, 56870.217419, 57161.294954, 57163.289352, 57164.288368, 57170.269792, 57171.226528, 57174.269838, 57175.267639, 57178.21669, 57179.210359, 57229.122269, 57230.092998, 57236.11956, 57972.057269, 57974.029387, 57978.112847, 57982.05941, 57985.05985, 57995.041134, 58004.083345, 58010.986111, 58015.978299, 58019.01772, 58022.018993, 58031.978565, 58035.993171, 58039.994757, 58044.983125]
mag_boch=[11.362845, 11.392318, 11.4020195, 11.417272, 11.555197, 11.627949, 11.619389, 11.3612, 11.60123, 11.589565, 11.491863, 11.463467, 11.590549, 11.594222, 11.372154, 11.361826, 11.333306, 11.443802, 11.562275, 11.468747, 11.54358, 11.571766, 11.51438, 11.481539, 11.409691, 11.388672, 11.320535, 11.323443, 11.483155, 11.519778, 11.570538, 11.437965, 11.445021, 11.500936, 11.501298, 11.470494, 11.431529, 11.51373, 11.426185, 11.419385, 11.463061, 11.477302, 11.440823, 11.420283, 11.457893, 11.485142]
time=np.remainder(T_boch,P)/P
ax1.scatter(time,mag_boch,marker=".",color="black")
#ax1.grid()
ax1.set_title("Bochum Disc Survey")
ax1.set_xlabel("Phase")
ax1.set_ylabel("SDSS i [mag]")
ax1.invert_yaxis()

ax2=axes.flatten()[1]
ogle=np.loadtxt("GD1070.18.22288.dat",skiprows=1)
time_ogle=np.remainder(ogle[:,0],P)/P
ax2.scatter(time_ogle,ogle[:,1],marker=".",color="black")
#ax2.grid()
ax2.set_title("OGLE")
ax2.set_ylabel("I [mag]")
ax2.set_xlabel("Phase")
ax2.invert_yaxis()

ax3=axes.flatten()[2]
data_V=pd.read_csv('ASASSN-V J170801.81-410255.6.csv')
#data_V2=pd.read_csv("AP52061023.csv")
#data_V=pd.concat((data_V1,data_V2))
phase_ASAS_V=np.remainder(data_V["hjd"]-2400000.5 ,P)/P
ax3.scatter(phase_ASAS_V,data_V["mag"],color="black",marker=".")
#ax3.grid()
ax3.set_title("ASAS-SN (V)")
ax3.set_ylabel("V [mag]")
ax3.set_xlabel("Phase")
ax3.invert_yaxis()

ax4=axes.flatten()[3]
data_g=pd.read_csv('light_curve_3407d752-9b3c-4ae0-8ae5-7da73ea2266b.csv')
data_g=data_g[data_g["mag"]<13.4]
data_g=data_g[data_g["mag_err"]<0.02]
phase_ASAS_g=np.remainder(data_g["HJD"]-2400000.5 ,P)/P
ax4.scatter(phase_ASAS_g,data_g["mag"].astype(float),marker=".",color="black")
#ax4.grid()
ax4.set_title("ASAS-SN (g)")
ax4.set_xlabel("Phase")
ax4.set_ylabel("SDSS g [mag]")
ax4.invert_yaxis()
plt.suptitle("GD1070.18.22288",fontsize=20)
plt.tight_layout()
plt.savefig("lc_comparsion.png",dpi=600)






