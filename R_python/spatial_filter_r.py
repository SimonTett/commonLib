# Provide spaital filtering using Simon Brown's R code.
import numpy as np
import xarray

import logging
import rpy2.situation
import typing
import  importlib.resources

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
my_logger=logging.getLogger(f"CPM_rain_analysis.{__name__}")
# get the r-packages installed.
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
# R package names
packnames = (['terra','viridis', 'fields', 'raster', 'codetools'])  # list of packages to install/make available,
# From example rpy2 install what needs to be installed.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    my_logger.warning(f"installing {names_to_install}")
    utils.install_packages(robjects.vectors.StrVector(names_to_install))
for package in packnames:
    rpackages.importr(package)  # so available

import rpy2.rinterface_lib
from rpy2.robjects import numpy2ri

from rpy2.robjects import default_converter

# Create a converter that starts with rpy2's default converter
# to which the numpy conversion rules are added.
np_cv_rules = default_converter + numpy2ri.converter  # rules for conversion

# read in and compile Simon Brown's code.
with importlib.resources.open_text('R_python',"filter_fn.R") as file: # open up the R source code.
    r_code = file.read()

filter_r = rpy2.robjects.packages.SignatureTranslatedAnonymousPackage(r_code, "filter")  # compile it!

# values for filtering

badth_values = dict(DJF=180, JJA=160)  # Seasonally dep precip for lines
badthp_values = dict(DJF=12, JJA=31)  # and for points.


def filter_data(var: xarray.Dataset,
                badth: float,
                badthp: float,
                badb: int = 2,
                bads: int = 1,
                max_precip: float = 10.0) -> (np.ndarray,bool):
    #
    """
    Apply Simon Brown's R filter function to a 2d numpy array. This filter attempts to remove artifacts due to
      non-conservation along the co-ordinate axis.
    :param var: An xarray variable
    :param badth: precipitation scaling (?) for lines.
    :param badthp: precipitation scaling (?) for points.
    :param badb: Number of pixels to set missing near bad data - big
    :param bads: Number of pixels to set missing near bad data - small
    :param max_precip: Only filter if any precip in domain is above this (mm/hr).

    returns filtered array (or original data)
    """
    if var.ndim != 2:
        raise ValueError("Must be a 2d array")

    if var.max() < max_precip:
        my_logger.debug("No filtering")
        return var.to_numpy()  # no filtering so just return the values.

    my_logger.debug(f"filtering badth={badth},badthp={badthp},badb={badb},bads={bads}")
    with rpy2.robjects.conversion.localconverter(np_cv_rules): # run the R code on the numpy array
        filtered = filter_r.filt_point_vbar_hbar_wnorm(var.to_numpy(), badth=badth, badthp=badthp, badb=badb, bads=bads)

    return filtered

def xarray_filter(data_array:xarray.DataArray,
                  badb: int = 2,
                  bads: int = 1,
                  max_precip: float = 10.0,
                  stack_dims:typing.Optional[typing.List[str]] = None,
                  ) -> (xarray.DataArray,xarray.Dataset):
    """

    :param data_array: DataArray to be filtered.
    :param badb: Number of pixels to set missing near bad data - big
    :param bads: Number of pixels to set missing near bad data - small
    :param max_precip: Only filter if any precip in domain is above this (mm/hr).
    :param stack_dims:  Dimensions to be stacked along. Default is time and ensemble_no. First dimension must be time-coord.
    :return: filtered dataArray and ancillary data as dataset.Ancillary dataset contains values for badth, badthp, badb, bads and max_precip
    """

    
    if stack_dims is None:
        stack_dims =['time','ensemble_member'] # default dimensions to stack over
    time_coord=stack_dims[0]

    result=[] # where we store the filtered result
    result_ancil=[] # where we store the ancillary data.
    for (c,g) in data_array.resample({time_coord:'QS-DEC'}): # iterate over seasons.
        seasons = set(g.time.dt.season.values) #  check only have one season.
        if len(seasons) != 1:
            raise ValueError(f"Can only cope with unique season. Have {seasons}")
        season = seasons.pop()  # get the season
        badth = badth_values.get(season, 170.0)  # Outside JJA or DJF use value 1/2 way between
        badthp = badthp_values.get(season, 22.0)  # Outside JJA or DJF use value 1/2 way between
        g.load()
        da = g.stack(stack_dim=stack_dims) # stack !
        da = da.groupby('stack_dim').\
            map(filter_data,args=(badth,badthp),badb=badb, bads=bads, max_precip=max_precip,shortcut=True) 
        # run the filtering -- one field at a time. Stacking this togetehr  afterwards.
        da = da.unstack('stack_dim') 
        # unstack. Causes some numpy RuntimeWarning's wiht modern version of pandas
        # make the ancillary data.
        timec=dict(time2=[c])
        ds=xarray.Dataset(
                dict(
                badth=xarray.DataArray(data=[badth],coords=timec),
                badthp=xarray.DataArray(data=[badthp],coords=timec),
                max_precip=xarray.DataArray(data=[max_precip],coords=timec),
                badb=xarray.DataArray(data=[badb],coords=timec),
                bads=xarray.DataArray(data=[bads],coords=timec),
                n_filter=da.isnull().sum(time_coord), # count bad.
                )
        )
        result_ancil.append(ds)
        result.append(da)
        logging.info(f"Processed data for {c}") 
    result = xarray.concat(result,dim='time')
    result_ancil = xarray.concat(result_ancil,dim='time2')
    result.attrs = data_array.attrs.copy() # explicitly copying the attributes.
    result.attrs['Processing'] = data_array.attrs.get('Processing',[])+['Filtered to remove non-conservation extreme rain']
    result_ancil.attrs = data_array.attrs.copy() #copy any attributes
    result_ancil.attrs['Processing'] = data_array.attrs.get('Processing', []) + [
        'Parameters for filtering']


    return result, result_ancil




