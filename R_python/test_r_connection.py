# scratch code to test if R stuff if working.
# get the r-packages installed.
import xarray
import CPMlib
from R_python import gev_r
import numpy as np
import pandas as pd
import rpy2.rinterface_lib
import rpy2.robjects as robjects
import rpy2.robjects.packages
from rpy2.robjects import numpy2ri
from rpy2.robjects import default_converter
import rpy2.robjects.pandas2ri as rpandas2ri

ext = rpy2.robjects.packages.importr('extRemes')
base = rpy2.robjects.packages.importr('base')
np_cv_rules = default_converter + numpy2ri.converter+ rpandas2ri.converter
#df=pd.DataFrame(np.array([10,20,40,78,90,100.]),columns=['x'])
radar_events = xarray.load_dataset(CPMlib.radar_dir/"radar_events_summary_5km_1hr_scotland.nc")
df=pd.DataFrame(radar_events.max_precip.sel(quantv=0.5),columns=['x'])
wts = np.ones(df.shape[0])#.squeeze()
#wts[29]=10.
wts[29]=100
wts[20]=100.
wts = robjects.vectors.FloatVector(wts) # convert to a R vector.
#robjects.globalenv['wt'] = wts
#x=df.values.squeeze()
x_symbol = robjects.r('x')
with rpy2.robjects.conversion.localconverter(np_cv_rules):
    #robjects.globalenv['df'] = df
    #r_result = robjects.r('fevd(x,data=df,weights=wt)')
    r_result = ext.fevd(x_symbol,data=df,weights=wts)
    r_sum = base.summary(r_result,silent=True)
    print(r_sum)

print('params ',r_sum.rx2('par'))