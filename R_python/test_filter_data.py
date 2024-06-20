# Code to test using Simon Brown's R code to filter data.
# Set up R
import numpy as np
import xarray
import CPMlib
import logging
from R_python import spatial_filter_r
import nc_time_axis


# load up the data.
path = (CPMlib.CPM_dir/"Example_CPM_data").glob("pr_rcp85_land-cpm_uk_2.2km_01_1hr_1994*.nc")
#ds = xarray.load_dataset(path)
ds = xarray.open_mfdataset(path,chunks=dict(time=24,grid_longitude=100,grid_latitude=100)).sortby('time')
logging.basicConfig(level=logging.DEBUG,force=True)
scotland = ds.pr.sel(grid_longitude=slice(357.5, 361), grid_latitude=slice(1.5, 7.5))  # Just want scotland
#scotland = ds.pr.sel(grid_longitude=slice(359.5, 360.5), grid_latitude=slice(5, 6.))#.isel(time=slice(0,48))  # Just want mini scotland
filter_scotland, ancil =spatial_filter_r.xarray_filter(scotland) # occasionally freezes.
max_rain = scotland.max(['grid_longitude','grid_latitude'])
max_rain_filtered = filter_scotland.max(['grid_longitude','grid_latitude'])
n_filter = filter_scotland.isnull().sum(['grid_longitude','grid_latitude'])
## plot the monthly max rain at each point.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#import commonLib

# plot max rain in domain
fig,axis=plt.subplots(figsize=(11,6),nrows=2,num='Domain_mx_rain',clear=True,sharex=True)
max_rain_filtered.plot(label='filtered',ax=axis[0])
max_rain.plot(label='raw',ax=axis[0])
# plot the number of filtered cells
n_filter.plot(ax=axis[1])
axis[1].set_title("Number cells set missing")
fig.legend()
fig.show()
#commonLib.saveFig(fig)


fig,axis = plt.subplots(nrows=1,ncols=3,subplot_kw=dict(projection=ccrs.OSGB()),
                        num='SB_filter_impact',clear=True,figsize=(11,8))
lev=np.arange(5,25) # levels
time_title = str(scotland.time[0].dt.strftime("%Y-%b").values)
for var,ax,title in zip([filter_scotland,scotland],axis.flatten(),['Filtered','Raw']):
    cm=var.max('time').plot(ax=ax,levels=lev,transform=CPMlib.projRot,add_colorbar=False,cmap='PuBu')
                         
    ax.set_title(f"{title} ")
    ax.set_title("Count missing ")
    ax.coastlines()
ax=axis.flatten()[-1]
filter_scotland.isnull().sum('time').plot(ax=ax, transform=CPMlib.projRot)
ax.set_title("Missing #")
ax.coastlines()

fig.suptitle(f"Max rain (mm/hr) for {time_title}")
fig.colorbar(cm,ax=axis,orientation='horizontal',fraction=0.1,aspect=40,pad=0.05)
#fig.tight_layout()
fig.show()
#commonLib.saveFig(fig)

