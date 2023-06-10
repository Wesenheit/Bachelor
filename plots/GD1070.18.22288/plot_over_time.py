import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
data=pd.read_csv("light_curve_eafe2ed7-6586-452d-a8a8-3f1b24b28dd6.csv")
data=data[data["Filter"]=="g"]
data=data[data["flux_err"]<1]
data=data[data["mag"].astype(float).values<13.50]
data=data[data["mag"].astype(float).values>12.7]
fig,ax=plt.subplots(figsize=(6,4))
ax.scatter(data["HJD"].astype(float).values-2450000,data["mag"].astype(float).values,
           marker=".",color="black",alpha=1,s=10)
ax.invert_yaxis()
ax.set_xlabel("Time [HJD-2450000]")
ax.set_ylabel("SDSS g [mag]")
plt.tight_layout()
plt.savefig("visibility_over_time.png",dpi=1000)