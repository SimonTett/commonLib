"""
SW to gev fit with R in python.
"""
import numpy
import xarray
import scipy.stats
import pandas as pd
import os
import numpy as np
import typing
import logging
import pathlib

my_logger = logging.getLogger(__name__)
use_weights_warn = True
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.pandas2ri as rpandas2ri
import rpy2.rinterface_lib
from rpy2.robjects import numpy2ri

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
# R package names
packnames = (['extRemes'])  # list of packages to install.
# From example rpy2 install what needs to be installed.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    print("installing ", names_to_install)
    utils.install_packages(robjects.vectors.StrVector(names_to_install))
for package in packnames:
    rpackages.importr(package)  # so available

base = rpackages.importr('base')  # to give us summary

type_array_or_float = typing.Union[np.ndarray, float]


def ks_fit(
        shape: float,
        location: type_array_or_float,
        scale: type_array_or_float,
        data: np.ndarray,
        dist: scipy.stats.rv_continuous = scipy.stats.genextreme
) -> float:
    """
    Do a Kolmogorov-Smirnov test on  data with varying location and scale params
    Data is normalised by subtracting location and dividing by scale. Not obvious how to normalise for varying shape!
    Answer is to use a constant shape!
    :param shape: shape parameter -- constant value
    :param location: location parameter  either a float or array of floats
    :param scale: scale parameter either a float or array of floats
    :param data: data to be tested
    :param dist: distribution to be tested against. Default is genextreme.
    :return: p-value of the test.
    TODO: extend to allow non constant data weights. See https://github.com/scipy/scipy/issues/12315.
     Issue seems to be test-stat rather than weighted calc.
     See https://doi.org/10.1017/CBO9780511977176.014 p358 for potential method.
    """
    data_norm = (data - location) / scale  # normalise to a standard distribution.
    dist = dist(shape, loc=0.0, scale=1.0)
    ks_stat = scipy.stats.kstest(data_norm, dist.cdf)  # TODO -- modify kstest to allow weights
    return ks_stat.pvalue


def gen_params(
        shape: float, location: float, scale: float,
        covariates: np.ndarray,
        d_shape: typing.Optional[np.ndarray] = None,
        d_location: typing.Optional[np.ndarray] = None,
        d_scale: typing.Optional[np.ndarray] = None,
) -> \
        (type_array_or_float, type_array_or_float, type_array_or_float):
    """
    Generate parameters for a GEV fit with co-variates.
    :param shape: shape parameter
    :param location: location parameter
    :param scale: scale parameter
    :param covariates: The co-variates.
    :param d_shape: change in shape parameter with covariates
    :param d_location: change in location parameter with covariates
    :param d_scale: change in scale parameter with covariates

    :return: shape, location, scale .
    """
    ncovariates = covariates.shape[0]
    result = dict(shape=shape, location=location, scale=scale)
    for name, array in zip(['d_shape', 'd_location', 'd_scale'], [d_shape, d_location, d_scale]):
        if array is None:
            continue  # skip processing.
        if array.size != ncovariates:
            raise ValueError(f"Expected {ncovariates} {name} values got {array.size}")
        key = name.replace('d_', '')
        result[key] += covariates.T.dot(array)

    return result["shape"], result["location"], result["scale"]


def gev_fit(
        *args: typing.List[np.ndarray],
        use_weights: bool = False,
        shapeCov: bool = False,
        minValues: int = 20,
        **kwargs
):
    """
    Do GEV fit using R and return named tuple of relevant values.
    :param x: Data to be fit
    :param cov: co-variate value (if None then not used)
    :param returnType. Type to return, allowed is:
        named -- return a named tuple (default)
        tuple -- return a tuple -- useful for apply_ufunc
        DataSet -- return a DataSet
    :param shapeCov -- If True allow the shape to vary with the co-variate.
    :return: A dataset of the parameters of the fit.
    #TODO rewrite this to directly call the R fn.
    """
    global use_weights_warn  # used to control printing out a warning message
    x = args[0]
    if use_weights:
        weights = args[-1]
        args = args[0:-1]
    else:
        weights = None  # to stop complaints from IDE about possibly undefined var
    ncov = len(args) - 1

    npts = 3 + 2 * ncov
    if shapeCov:
        npts += ncov
    # end of dealing with covariate and trying to get to the right shape.
    params = np.broadcast_to(np.nan, npts)
    se = np.broadcast_to(np.nan, npts)
    cov_params = np.broadcast_to(np.nan, [npts, npts])
    nllh = np.array([np.nan])
    aic = np.array([np.nan])
    ks = np.array([np.nan])
    L = ~np.isnan(x)
    # check we have enough data.
    sumOK = L.sum()
    if sumOK < minValues:
        my_logger.warning(f'Not enough data for fit. Have {sumOK} need {minValues}')
        return (params, se, cov_params, nllh, aic, ks)

    df_data = [x[L]]  # remove nan from data]
    cols = ['x']

    r_code = 'fevd(x=x,data=df'
    for indx, cov in enumerate(args[1:]):
        if np.isnan(cov[L]).any():
            ValueError('Missing values in covariate')
        df_data.append(cov[L])  # remove places where x was nan from cov.
        cols.append(f"cov{indx:d}")
    if len(cols) > 1:
        cov_expr = "~" + " + ~".join(cols[1:])  # expression for covariances
        r_code = r_code + ",location.fun=" + cov_expr + ",scale.fun=" + cov_expr
        if shapeCov:
            r_code += ',shape.fun=' + cov_expr
    if use_weights:  # got weights so include that in the R code.
        r_code += ',weights=wt'
        if np.isnan(weights[L]).any():
            ValueError('Missing values in weights')

        wts = robjects.vectors.FloatVector(weights[L])  # convert to a R vector.
        robjects.globalenv['wt'] = wts  # and push into the R environment.

    r_code += ')'  # add on the trailing bracket.
    df = pd.DataFrame(np.array(df_data).T, columns=cols)
    # check for missing and raise an error if any
    if df.isnull().any().any():
        ValueError('Missing values in data')
    with (robjects.default_converter + rpandas2ri.converter + numpy2ri.converter).context():
        robjects.globalenv['df'] = df  # push the dataframe with info into R

    try:
        r_fit = robjects.r(r_code)  # do the fit
        fit = base.summary(r_fit, silent=True)  # get the summary fit info
        # extract the data
        params = fit.rx2('par')  # get the parameters.
        se = fit.rx2('se.theta')  # get the std error
        cov_params = fit.rx2('cov.theta')  # get the covariance.
        if isinstance(params, rpy2.rinterface_lib.sexp.NULLType):
            # params not present (for some reason) set everything to nan
            params = np.broadcast_to(np.nan, npts)
        if isinstance(se, rpy2.rinterface_lib.sexp.NULLType):
            # std err not present (for some reason) set everything to nan
            se = np.broadcast_to(np.nan, npts)
        if isinstance(cov_params, rpy2.rinterface_lib.sexp.NULLType):
            # cov not present (for some reason) set everything to nan
            cov_params = np.broadcast_to(np.nan, [npts, npts])
        params = np.array(params)
        se = np.array(se)
        cov_params = np.array(cov_params)

        # ordering is:
        # loc, loc-cov1, loc-cov2,... scale, scale-cov1, scale-cov-2, shape, [shape-cov1 etc]
        # negate shape params as R and python conventions differ.
        start_shape = -1
        if shapeCov:
            start_shape = -ncov
        params[start_shape:] = params[start_shape:] * (-1)
        nllh = np.array(fit.rx2('nllh'))
        aic = np.array(fit.rx2('AIC'))
        #TODO extend python ks-test to cope with weights see https://doi.org/10.1017/CBO9780511977176.014 p358 for method.
        if not use_weights:  # compute the k-s fit if we are not weighting
            if ncov != 0:  # got some covariance
                if shapeCov:
                    d_shape = params[3 + 2 * ncov:]
                else:
                    d_shape = None
                d_location = params[1:1 + ncov]
                d_scale = params[2 + ncov:2 + 2 * ncov]
            else:
                d_shape = None
                d_location = None
                d_scale = None
            shape, location, scale = gen_params(float(params[start_shape]), float(params[0]), float(params[1 + ncov])
                                                , np.array(df_data[1:]),  # covariates
                                                d_shape=d_shape, d_location=d_location, d_scale=d_scale
                                                )
            ks = np.array(ks_fit(shape, location, scale, df_data[0]))
        elif use_weights_warn:
            my_logger.warning("Weights are not used in KS test")
            use_weights_warn = False  # only want the message once!

    except rpy2.rinterface_lib.embedded.RRuntimeError:
        # some problem in R with the fit. Complain

        my_logger.warning("Runtime error with R when doing GEV fit")

    return (params, se, cov_params, nllh, aic, ks)  # return the data.


import scipy.stats


def gev_fit_python(
        data: np.ndarray, nmin: int = 10,
        **kwargs
        ) -> (np.ndarray, np.ndarray):
    L = ~np.isnan(data)  # check nan
    if L.sum() <= nmin:  # Not enough data
        fit = np.array([np.nan] * 3)  # return nans
        ks = np.array([np.nan] * 2)
    else:
        d = data[L]
        fit = scipy.stats.genextreme.fit(d, **kwargs)  # in order shape,location, scale
        # do ks-test
        dist = scipy.stats.genextreme(*fit)
        ks = scipy.stats.kstest(d, dist.cdf)
        # convert to numpy arrays
        ks = np.array([ks.statistic, ks.pvalue])
        fit = np.array(fit)
    return fit, ks


def xarray_gev_python(
        ds: xarray.Dataset,
        dim: str = 'time_ensemble',
        file: typing.Optional[pathlib.Path] = None,
        recreate_fit: bool = False,
        **kwargs
        ) -> xarray.Dataset:
    """
    Fit a GEV to xarray data using scipy.stats. Less powerful than R. Note this code has not been solidly tested.
    :param name: Name of dataArray
    :param ds: dataset for which GEV is to be fit
    :param dim: The dimension over which to collapse.
    :param file -- If defined save fit to this file. If file exists then read data from it and so do not actually do fit.
    :param recreate_fit -- If True even if file exists compute fit.
    :param kwargs: any kwargs passed through to the fitting function
    :return: a dataset containing:
        Parameters -- the parameters of the fit; location, scale, shape
        ks_result -- the statistic and pvalue from the kstest.
    """
    if (file is not None) and file.exists() and (
            not recreate_fit):  # got a file specified, it exists and we are not recreating fit
        fit = xarray.load_dataset(file)  # just load the dataset and return it
        my_logger.info(f"Loaded existing data from {file}")
        return fit

    my_logger.debug('Doing fit')
    params, ks_result = xarray.apply_ufunc(
        gev_fit_python, ds,
        input_core_dims=[[dim]],
        output_core_dims=[['parameter'], ['ks']],
        vectorize=True, kwargs=kwargs
    )
    my_logger.debug('Done fit. Making dataset')
    pnames = ['shape', 'location', 'scale']
    fit_names = ['statistic', 'pvalue']

    params = params.rename("Parameters").assign_coords(parameter=pnames)
    ks_result = ks_result.rename("ks_result").assign_coords(ks=fit_names)
    fit = xarray.merge([params, ks_result])
    if file is not None:
        fit.to_netcdf(file)  # save the dataset.
        my_logger.info(f"Wrote fit information to {file}")
    return fit


def xarray_gev(
        data_array: xarray.DataArray,
        cov: typing.Optional[typing.List[xarray.DataArray]|xarray.DataArray] = None,
        shape_cov=False,
        dim: [typing.List[str], str] = 'time_ensemble',
        file=None, recreate_fit: bool = False,
        verbose: bool = False,
        name: typing.Optional[str] = None,
        weights: typing.Optional[xarray.DataArray] = None,
        extra_attrs:typing.Optional[dict] = None,
        **kwargs
):
    #
    """
    Fit a GEV to xarray data using R.

    :param data_array: dataArray for which GEV is to be fit
    :param cov: covariate (If None not used) --a list of dataarrays or None.
    :param shape_cov: If True then allow the shape to vary with the covariate.
    :param weights: Weights for each sample. If not specified, no weighting will be done.
    :param dim: The dimension(s) over which to collapse.
    :param file -- if defined save fit to this file. If file exists then read data from it and so not actually do fit.
    :param recreate_fit -- if True even if file exists compute fit.
    :param verbose -- be verbose if True
    :param name: Name of the fit. Stored in result attributes under name.
    :param kwargs: any kwargs passed through to the fitting function
    :return: a dataset containing:
        Parameters -- the parameters of the fit; location, location wrt cov, scale, scale wrt cov, shape, shape wrt cov
        StdErr -- the standard error of the fit -- same parameters as Parameters
        nll -- negative log likelihood of the fit -- measure of the quality of the fit
        AIC -- aitkin information criteria.
        ks -- KS test result
    """
    if (file is not None) and file.exists() and (
            not recreate_fit):  # got a file specified, it exists and we are not recreating fit
        data_array = xarray.load_dataset(file)  # just load the dataset and return it
        if verbose:
            print(f"Loaded existing data from {file}")
        return data_array

    kwargs['shapeCov'] = shape_cov

    if cov is None:
        cov = []
    if not isinstance(cov, list):
        cov = [cov]

    cov_names = [c.name for c in cov]

    ncov = len(cov)
    if isinstance(dim, str):
        dim = [dim]
    input_core_dims = [dim] * (1 + ncov)

    output_core_dims = [['parameter']] * 2 + [['parameter', 'parameter2'], [], [], []]
    gev_args = [data_array] + [c.broadcast_like(data_array) for c in cov] # broadcast covariates to data_array
    if weights is not None:
        gev_args += [weights]
        input_core_dims += [dim]
        kwargs.update(use_weights=True)


    params, std_err, cov_param, nll, AIC, ks = xarray.apply_ufunc(gev_fit, *gev_args,
                                                                  input_core_dims=input_core_dims,
                                                                  output_core_dims=output_core_dims,
                                                                  vectorize=True, kwargs=kwargs
                                                                  )
    pnames = []
    for n in ['location', 'scale', 'shape']:
        pnames += [n]
        if not shape_cov and n == 'shape':
            continue  # no Dshape_dX
        for cn in cov_names:
            pnames += ['D' + n + '_' + cn]

    # name variables and then combine into one dataset.

    params = params.rename("Parameters")
    std_err = std_err.rename("StdErr")
    cov_param = cov_param.rename("Cov")
    nll = nll.rename('nll').squeeze()
    AIC = AIC.rename('AIC').squeeze()
    ks = ks.rename('ks').squeeze()
    data_array = xarray.Dataset(dict(Parameters=params, StdErr=std_err, Cov=cov_param, nll=nll, AIC=AIC, KS=ks)
                                ).assign_coords(
        parameter=pnames, parameter2=pnames
    )
    if name:
        data_array.attrs.update(name=name)
    if extra_attrs:
        data_array.attrs.update(extra_attrs)
    if file is not None:
        file.parent.mkdir(exist_ok=True, parents=True)  # make directory
        data_array.to_netcdf(file)  # save the dataset.
        if verbose:
            print(f"Wrote fit information to {file}")
    return data_array


## use apply ufunc to generate distributions...

def fn_isf(c, loc, scale, p: typing.Optional[np.ndarray] = None, dist=scipy.stats.genextreme):
    shape = list(c.shape) + list(p.shape)
    fd = dist(np.broadcast_to(c, shape), loc=np.broadcast_to(loc, shape), scale=np.broadcast_to(scale, shape))
    # handle single p value.
    # if len(p) == 1:
    #    p=p.reshape(1,-1)

    # breakpoint()
    x = fd.isf(np.broadcast_to(p, shape))
    # x=fd.isf(p)
    # x = fdist.isf(p)  # values for 1-cdf.
    return x


def fn_sf(c, loc, scale, x=None, dist=scipy.stats.genextreme):
    p = dist.sf(x, c, loc=loc, scale=scale)
    # p = fdist.sf(x)  # 1-cdf for given x
    return p


def fn_interval(c, loc, scale, alpha=None, dist=scipy.stats.genextreme):
    # fdist = dist(c, loc=loc, scale=scale)
    range = dist.interval(alpha, c, loc=loc, scale=scale)  # range for dist
    return np.array([range[0], range[1]])


def xarray_sf(x, params, output_dim_name='value'):
    """
    Compute the survival value for different values based on dataframe of fit parameters.
    :param params: xarray dataarray of shape, location and scale values
    :param output_dim_name: name of output dimension. Default is "value" but set it to what ever you are using. E.g "Rx1hr"
    :param kwargs: passed to fn_sf which does the computation. Must contain x which is used for the computation.
    :return:dataset of survival function values (1-cdf for values specified)

    """
    # need to add a dummy singleton dimension to params
    params = params.assign_coords(probability=1)
    sf = xarray.apply_ufunc(fn_sf, params.sel(parameter='shape'), params.sel(parameter='location'),
                            params.sel(parameter='scale'),
                            output_core_dims=[[output_dim_name]],
                            vectorize=True, kwargs=dict(x=x)
                            )
    sf = sf.assign_coords({output_dim_name: x}).rename('sf')

    return sf


def xarray_interval(alpha, params):
    """
    Compute the interval for different values based on dataframe of fit parameters.
    :param params: xarray dataarray of shape, location and scale values
    :param alpha -- alpha value for interval fn.
    :return:dataset of intervals
    """
    interval = xarray.apply_ufunc(fn_interval, params.sel(parameter='shape'), params.sel(parameter='location'),
                                  params.sel(parameter='scale'),
                                  output_core_dims=[['interval']],
                                  vectorize=True, kwargs=dict(alpha=alpha)
                                  )
    offset = (1 - alpha) / 2
    interval = interval.Parameters.assign_coords(interval=[offset, 1 - offset]).rename('interval')

    return interval


def xarray_isf(
        p: np.ndarray,
        params: xarray.DataArray,
        output_dim_name: typing.Optional[str] = None,
        input_dim_name: typing.Optional[str] = None
) -> xarray.DataArray:
    """
    Compute the inverse survival function for specified probability values
    :param output_dim_name: name of output_dim -- default is probability
    :param params: dataset of parameter values.
    :param p: sf values at which values to be computed
    :param kwargs:Additional keyword arguments passes to fn_isf. Make sure p is set.
    :return:
    """
    if output_dim_name is None:
        output_dim_name = 'probability'
    if input_dim_name is None:
        input_dim_name = 'quantv'

    # FIXME -- this is failing with a broadcast error. Sigh. I hate apply_ufunc.
    aa = tuple([params.sel(parameter=k, drop=True) for k in ['shape', 'location', 'scale']])
    x = xarray.apply_ufunc(fn_isf, *aa, input_core_dims=[[input_dim_name]] * len(aa),
                           output_core_dims=[[output_dim_name, input_dim_name]],
                           vectorize=True, kwargs=dict(p=p)
                           )
    x = x.assign_coords({output_dim_name: p}).rename('isf')
    return x


def param_cov(params, cov):
    raise Warning("Use param_at_cov")
    p = ['location', 'scale', 'shape']
    p2 = ["D" + a.lower() for a in p]
    params_c = params.Parameters.sel(parameter=p2).assign_coords(parameter=p) * cov + params.Parameters.sel(parameter=p)
    return params_c


def param_at_cov(params, cov):
    p = ['location', 'scale', 'shape']
    p2 = ["D" + a.lower() for a in p]
    params_c = params.sel(parameter=p2).assign_coords(parameter=p) * cov + params.sel(parameter=p)
    return params_c


def xarray_gev_isf(
        params: xarray.DataArray,
        pvalues: typing.Union[np.ndarray, typing.List[float]],
        distribution: typing.Optional[scipy.stats.rv_continuous] = None, ) -> xarray.DataArray:
    """

    :param distribution: distribution to be used for fit
    :param params: dataArray of parameters with co-ords parameter and names shape, location and scale
    :param pvalues: probability values for which thresholds are computed.
    :return: dataarray
    """
    if distribution is None:
        distribution = scipy.stats.genextreme
    # convert list to np.ndarray
    if isinstance(pvalues, list):
        pvalues = np.array(pvalues)  # convert to a numpy array
    pvalues = np.unique(pvalues)  # get the unique values.
    # extract data expanding dimension and generate frozen dist.
    p_list = [params.sel(parameter=p).expand_dims(dim=dict(pvalues=1), axis=-1) for p in ['shape', 'location', 'scale']]
    fit = distribution(*p_list)
    # compute the inverse sf (pvalue -> threshold)
    result = fit.isf(np.expand_dims(pvalues, axis=list(range(0, params.ndim - 1))))
    # now to make it into a DataArray
    coords = {coord: params[coord] for coord in params.coords if coord != 'parameter'}
    dims = [dim for dim in params.dims if dim != 'parameter']
    coords['pvalues'] = pvalues
    dims.append('pvalues')
    result = xarray.DataArray(data=result, coords=coords, dims=dims, name='isf')
    # add on another co-ord -- the return_value.
    result = result.assign_coords(return_period=('pvalues', 1.0 / pvalues))
    # could do more with meta-data but that is for another time,
    return result


def xarray_gev_sf(
        params: xarray.DataArray,
        thresholds: typing.Union[np.ndarray, typing.List[float], float],
        distribution: typing.Optional[scipy.stats.rv_continuous] = None, ) -> xarray.DataArray:
    """

    :param distribution: distribution to be used for fit
    :param params: dataArray of parameters with co-ords parameter and names shape, location and scale
    :param thresholds: probability values for which thresholds are computed.
    :return: dataarray
    """
    if distribution is None:
        distribution = scipy.stats.genextreme
    # convert list to np.ndarray
    if isinstance(thresholds, (list, float)):
        thresholds = np.array(thresholds)  # convert to a numpy array
    thresholds = np.unique(thresholds)  # get the unique values.
    # extract data expanding dimension and generate frozen dist.
    p_list = [params.sel(parameter=p).expand_dims(dim=dict(pvalues=1), axis=-1) for p in ['shape', 'location', 'scale']]
    fit = distribution(*p_list)
    # compute the  sf (threshold -> sf)
    result = fit.sf(np.expand_dims(thresholds, axis=list(range(0, params.ndim - 1))))
    # now to make it into a DataArray
    coords = {coord: params[coord] for coord in params.coords if coord != 'parameter'}
    dims = [dim for dim in params.dims if dim != 'parameter']
    coords['threshold'] = thresholds
    dims.append('threshold')
    result = xarray.DataArray(data=result, coords=coords, dims=dims, name='isf')
    # could do more with meta-data but that is for another time,
    return result
