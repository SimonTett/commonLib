import cartopy.feature

import GSHHS_WDBII
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

GSHSS = GSHHS_WDBII.GSHHS_WDBII()
extent = (-15, 30, 20, 61)
fig, axis = plt.subplots(nrows=2, ncols=2, figsize=[12, 8],
                       subplot_kw={'projection': ccrs.PlateCarree()}, clear=True,
                       num=f'ghss_example_europe',sharex='all',sharey='all')
for resoln,ax in zip(GSHSS.allowed_scales[1:5],axis.flatten()):
    # generate the features we want
    rivers = GSHSS.rivers(scale=resoln)
    canals = GSHSS.canals(scale=resoln)
    inter_rivers = GSHSS.inter_rivers(scale=resoln)
    coastlines = GSHSS.coastlines(scale=resoln)
    lakes = GSHSS.lakes(scale=resoln)

    ax.set_extent(extent) # set the extent
    ax.stock_img()

    # plot them
    for f in [coastlines, lakes, rivers,canals,inter_rivers]:
        if f is not None:
            ax.add_feature(f)
    ax.set_title(f"GHSS at {resoln} Resolution")
    ax.gridlines(zorder=10)
    ax.add_feature(cartopy.feature.BORDERS,edgecolor='red',zorder=20)

fig.tight_layout()
fig.show()

# now make world wide plot.
fig, axis = plt.subplots(nrows=2, ncols=1, figsize=[12, 8],
                       subplot_kw={'projection': ccrs.PlateCarree()}, clear=True,
                       num=f'ghss_example_world',sharex='all',sharey='all')
for resoln,ax in zip(GSHSS.allowed_scales[3:],axis.flatten()):
    # generate the features we want
    rivers = GSHSS.rivers(scale=resoln,zorder=30)
    canals = GSHSS.canals(scale=resoln)
    inter_rivers = GSHSS.inter_rivers(scale=resoln)
    coastlines = GSHSS.coastlines(scale=resoln)
    lakes = GSHSS.lakes(scale=resoln)
    ax.set_global() # set the extent
    ax.stock_img()

    # plot them
    for f in [coastlines, lakes, rivers,canals,inter_rivers]:
        ax.add_feature(f)
    ax.set_title(f"GHSS at {resoln} Resolution")
    ax.gridlines(zorder=10)
    ax.add_feature(cartopy.feature.BORDERS,edgecolor='red',zorder=20)

fig.tight_layout()
fig.show()

