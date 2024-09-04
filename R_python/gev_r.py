"""
SW to gev fit with R in python.
"""
import logging
import pathlib
import typing

import numpy as np
import pandas as pd
import scipy.stats
import xarray

my_logger = logging.getLogger(__name__)
use_weights_warn = True
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.pandas2ri as rpandas2ri
import rpy2.rinterface_lib
from rpy2.robjects import numpy2ri
import re
import importlib.resources
import scipy.stats

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
# R package names
packnames = (['extRemes'])  # list of packages to install.
# From example rpy2 install what needs to be installed.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    my_logger.info(f"installing " + " ".join(names_to_install))
    utils.install_packages(robjects.vectors.StrVector(names_to_install))
for package in packnames:
    rpackages.importr(package)  # so available

base = rpackages.importr('base')  # to give us summary
type_array_or_float = typing.Union[np.ndarray, float]
with (importlib.resources.files('R_python') / "fevd_summ.R").open('rt') as file:  # open up the R source code.
    fevd_sum_code = file.read()
r_gev = rpy2.robjects.packages.SignatureTranslatedAnonymousPackage(fevd_sum_code, "r_gev")  # compile it!


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
        my_logger.debug(f'Not enough data for fit. Have {sumOK} need {minValues}')
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


def gev_fit_new(df:pd.DataFrame,
                shape_cov: bool = False,
                min_values: int = 20,
                verbose: bool = False,
                weights: typing.Optional[np.ndarray] = None,
                initial: typing.Optional[dict[str, typing.Union[float, int,  list]]] = None,
                **kwargs) -> tuple[pd.Series, pd.Series, pd.DataFrame, float, float, float]:
    """
        Do GEV fit using R and return values.
        :param data_arrays : Data to be fit, covariate values and optionally weights.
            data_array[0] is data to be fit.
            data_array[1:] are the covariate values.
            If weights is True then data_array[-1] is the weights.
        :param shape_cov : If True allows the shape to vary with the co-variate.
        :param min_values : Minimum number of values to do the fit.
        :param verbose : If True be verbose -- passed to the R code.
        :param weights : Weights for array
        :param initial : Initial values for the fit. Should contain initial values for fevd.
                Contains location , scale,shape -- see doc for extRemes::fevd
                Will negate any shape values as R and python conventions differ.
        :param shape_cov -- If True allow the shape to vary with the co-variate.
        :return: the parameter values -- as a pandas series,
                 the std error -- as a pandas series,
                 the covariance matrix -- as a pandas data frame,
                    the negative log likelihood, as a float
                    the AIC, as a float
                    the ks test p value as a float
        """

    def fail_result(shapes:dict[str,list], cov_names:typing.Optional[list[str]]=None, shape_cov: bool = False) -> dict:
        """
        Generate a failure result
        :return: dict of results all with nan values
        """
        result = {name: np.nan if shp == (1,) else np.broadcast_to(np.nan, shp) for name, shp in shapes.items()}
        # deal with the dataframes/series

        index = ['location', 'scale', 'shape']
        for n in cov_names:
            index += [f'Dlocation_{n}', f'Dscale_{n}']
            if shape_cov:
                index += [f'Dshape_{n}']
        result['par']: pd.Series = pd.Series(result['par'], index=index)
        result['se_theta']: pd.Series = pd.Series(result['se_theta'], index=index)
        result['cov_theta']: pd.DataFrame = pd.DataFrame(result['cov_theta'], index=index, columns=index)
        return result

    global use_weights_warn  # used to control printing out a warning message

    if 'use_phi' in kwargs:
        raise NotImplementedError('use_phi not implemented')
    ncov = len(df.columns) - 1  # number of covariates
    npts = 3 + 2 * ncov  # number of points in the result
    if shape_cov:
        npts += ncov  # add in the shape covariates
    if ncov > 0:
        cov_names = df.columns[1:].tolist()  # names of the covariates
    else:
        cov_names = []

    L = ~np.isnan(df.iloc[:,0])  # mask
    # shapes of the results. par is the parameter values, se_theta is the standard errors, cov_theta is the covariance matrix
    shapes = dict(par=npts, se_theta=npts, cov_theta=[npts, npts], nllh=[1], AIC=[1], ks=[1])
    # check we have enough data.
    sumOK = L.sum()
    if sumOK < min_values:
        my_logger.warning(f'Not enough data for fit. Have {sumOK} need {min_values}')
        result = fail_result(shapes,cov_names=cov_names,shape_cov=shape_cov)
        return tuple(result.values())
    # check for missing and raise an error if any
    bad_cols=[]
    for name,col in df.iloc[:,1:].items():
        if np.isnan(col[L]).any():
            bad_cols.append(name)
    if bad_cols:
        ValueError(f'Missing values in covariate for cols {bad_cols}')
    r_args = kwargs.copy()
    # deal with things we know how to deal with which we also use in the python code or need to explicitly convert for R.
    r_args['verbose'] = verbose
    if weights is not None:
        # check weights has the right size and if not complain.
        if weights.shape != df.iloc[:,0].shape:
            raise ValueError('Weights have wrong shape')
        r_args['weights'] = robjects.vectors.FloatVector(weights[L])  # extract weights and make it an r floatVector
    if initial:
        r_initial = initial.copy()  # don't want to change the original
        # convert all values to a float
        shape = initial.get('shape')
        if shape is not None and isinstance(shape, (float, int, np.ndarray)):
            r_initial['shape'] = -shape
        if shape is not None and isinstance(shape, list):
            r_initial['shape'] = [-v for v in shape]

        r_args['initial'] = robjects.r.list(**r_initial)  # convert the initial values to an R list


    if ncov > 0:
        cov_expr = "~" + " + ~".join(df.columns[1:])  # expression for covariances
        cov_expr = robjects.Formula(cov_expr)
        r_args['location.fun'] = cov_expr
        r_args['scale.fun'] = cov_expr
        if shape_cov:
            r_args['shape.fun'] = cov_expr

    try:
        with (robjects.default_converter + rpandas2ri.converter + numpy2ri.converter).context():
            fit = r_gev.fevd_sum(x=df.columns[0], data=df[L], **r_args)  # run the R
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        my_logger.warning(f'Error in R code: {e}')
        result = fail_result(shapes,cov_names,shape_cov=shape_cov)
        return tuple(result.values())
    # annoyingly the names that R generates change if covariates are used so  need to rename them.
    par_names = [f.replace('mu', 'location', 1).replace('sigma', 'scale', 1) for f in fit['par.names']]
    # and then any xxx0 have the 0 removed
    pattern = '^(location|scale|shape)0$'
    par_names = [re.sub(pattern, r'\1', f) for f in par_names]  # const values
    # and name the covariates Dlocation_cname[0], Dlocation_cname[1] etc

    for indx, cname in enumerate(cov_names):
        pattern = fr'^(location|scale|shape){indx + 1:d}$'
        par_names = [re.sub(pattern, rf'D\1_{cname}', f) for f in par_names]  # covariate values

    # now to extract the values and convert them to series, dataframes or floats as appropriate
    result = dict()
    for k in shapes.keys():
        shp = shapes[k]
        f = fit.get(k.replace('_', '.'))
        if (f is None) or isinstance(f, rpy2.rinterface_lib.sexp.NULLType):
            # Not set so set all to nan
            f = np.broadcast_to(np.nan, shp)

        f = np.array(f)
        if shapes[k] == npts:  # make vectors a series with names
            f = pd.Series(f, index=par_names).rename(k)
        elif shapes[k] == [npts, npts]:  # make a data frame
            f = pd.DataFrame(f, index=par_names, columns=par_names)
        else:
            pass
        # shapes are opposite sign in R and python so negate them in par and cov_theta
        if k == 'par':
            f = f.mask(f.index.str.startswith('shape'), -f)
        elif k == 'cov_theta':
            L = f.index.str.startswith('shape')
            f = f.mask(np.broadcast_to(L, f.shape), -f)
            L = f.columns.str.startswith('shape')
            f = f.mask(np.broadcast_to(L, f.shape), -f)
        else:
            pass
        result[k] = f
    # and do  the ks test
    # TODO extend python ks-test to cope with weights see https://doi.org/10.1017/CBO9780511977176.014 p358 for method.
    if  weights is None:  # compute the k-s fit if we are not weighting
        params = result['par']
        d_shape = None
        d_location = None
        d_scale = None
        if ncov != 0:  # got some covariance
            # get the covariate values
            d_location = params.loc[[f'Dlocation_{d}' for d in cov_names]]
            d_scale = params.loc[[f'Dscale_{d}' for d in cov_names]]
            if shape_cov:
                d_shape = params.loc[[f'Dshape_{d}' for d in cov_names]]

        shape, location, scale = gen_params(float(params.loc['shape']), float(params.loc['location']),
                                            float(params.loc['scale'])
                                            , np.array(df.iloc[:, 1:].T),  # covariates
                                            d_shape=d_shape, d_location=d_location, d_scale=d_scale
                                            )
        result['ks'] = np.array([ks_fit(shape, location, scale, df.iloc[:, 0])])
    elif use_weights_warn:
        my_logger.warning("Weights are not used in KS test")
        use_weights_warn = False  # only want the message once!

    for k, v in result.items():
        if v.shape == (1,):
            result[k] = float(np.squeeze(v))

    return tuple(result.values())


def gev_fit_wrapper(*data_arrays:tuple[np.ndarray],
                    shape_cov: bool = False,
                    min_values: int = 20,
                    verbose: bool = False,
                    weights: bool = False,
                    cov_names: typing.Optional[typing.List[str]] = None,
                    initial: bool = False,
                    **kwargs) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
        If inital is set then data_arrays[-1] is the initial values for the fit.
        The names of the parameters are in kwargs['initial_names'] which must be present.
        Sigh -- the joys of apply_unfuc.
        Calls gev_fit_new and converts the results to numpy arrays. See gev_fit_new for other  details.
        :returns: the parameter values, the std error, the covariance matrix, the negative log likelihood, the AIC, the ks test p value and names of the parameters.
    """
    # need to convert data_arrays to a list
    data_arrays=list(data_arrays)
    if initial:
        init = [[float(d)] for d in data_arrays.pop()]
        keys = kwargs.pop('initial_names')
        initial = dict(zip(keys,init))

        # Not clean at all... Arises because apply_ufunc can only broadcast its positional arguments
        # Now to handle Dlocation_xx etc which are added to the list for location|shape|scale
        for key in keys:
            if re.match(r'D(location|shape|scale)_',key):
                param = key.split('_')[0][1:] # real hack -- keys are Dxxxx_
                initial[param] = initial[param]+initial.pop(key)
    else:
        initial = None
    if weights:
        weights = data_arrays.pop()
    else:
        weights = None

    ncov = len(data_arrays) - 1
    npts = 3 + 2 * ncov
    if shape_cov:
        npts += ncov
    if cov_names is None:
        cov_names = [f'cov{i:d}' for i in range(1, ncov + 1)]

    # now to construct the dataframe which will be passed to the R code.
    # generate the data frame
    df = [pd.Series(data_arrays[0], name='x')]
    for array, name in zip(data_arrays, ['x']+cov_names):
        df.append(pd.Series(array).rename(name))  # add in the covariates.
    df= pd.DataFrame(np.stack(data_arrays, axis=1), columns=['x'] + cov_names)
    result = gev_fit_new(df, shape_cov=shape_cov, min_values=min_values, verbose=verbose, weights=weights,
                         initial=initial,  **kwargs)
    par_names = result[0].index.values
    result = [np.array(r) for r in result]
    result += [par_names]
    return tuple(result)


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
        use_dask: bool = False,
        **kwargs
) -> xarray.Dataset:
    """
    Fit a GEV to xarray data using scipy.stats. Less powerful than R. Note this code has not been solidly tested.
    :param name: Name of dataArray
    :param ds: dataset for which GEV is to be fit
    :param dim: The dimension over which to collapse.
    :param file -- If defined save fit to this file. If file exists then read data from it and so do not actually do fit.
    :param recreate_fit -- If True even if file exists compute fit.
    :oaram use_dask: If True use dask to parallelize the computation.
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

    if use_dask:
        extra_args = dict(dask='parallelized',
                          dask_gufunc_kwargs=dict(allow_rechunk=True),
                          output_sizes=dict(parameter=3, ks=2))
    else:
        extra_args = {}
    params, ks_result = xarray.apply_ufunc(
        gev_fit_python, ds,
        input_core_dims=[[dim]],
        output_core_dims=[['parameter'], ['ks']],
        vectorize=True, kwargs=kwargs, **extra_args
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
        cov: typing.Optional[typing.Union[typing.List[xarray.DataArray],xarray.DataArray]] = None,
        shape_cov=False,
        dim: [typing.List[str], str] = 'time_ensemble',
        file=None, recreate_fit: bool = False,
        verbose: bool = False,
        name: typing.Optional[str] = None,
        weights: typing.Optional[xarray.DataArray] = None,
        initial: typing.Optional[xarray.DataArray] = None,
        extra_attrs: typing.Optional[dict] = None,
        use_dask: bool = False,
        **kwargs
):
    #
    """
    Fit a GEV to xarray data using R.

    :param data_array: dataArray for which GEV is to be fit
    :param cov: covariate (If None not used) --a list of dataarrays or None. Will be broadcast to data_array.
    :param shape_cov: If True then allow the shape to vary with the covariate.
    :param weights: Weights for each sample. If not specified, no weighting will be done.
    :param dim: The dimension(s) over which to collapse.
    :param file -- if defined save fit to this file. If file exists then read data from it and so not actually do fit.
    :param recreate_fit -- if True even if file exists compute fit.
    :param verbose -- be verbose if True
    :param name: Name of the fit. Stored in result attributes under name.
    :param use_dask: If True use dask to do the fitting.
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
        my_logger.info(f"Loaded existing data from {file}")
        return data_array

    kwargs['shape_cov'] = shape_cov
    kwargs['verbose'] = verbose
    dim_not_collapsed = set(data_array.dims) - set([dim])  # dimensions which are not collapsed
    if cov is None:
        cov = []
    if not isinstance(cov, list):
        cov = [cov]

    cov_names = [c.name for c in cov]

    ncov = len(cov)
    if isinstance(dim, str):
        dim = [dim]
    input_core_dims = [dim] * (1 + ncov)

    output_core_dims = [['parameter']] * 2 + [['parameter', 'parameter2'], [], [], [], ['parameter']]
    my_logger.debug('Setting up GEV args')
    gev_args = [data_array] + [c.broadcast_like(data_array) for c in cov]  # broadcast covariates to data_array

    if use_dask:
        extra_args = dict(dask='parallelized',
                          dask_gufunc_kwargs=dict(allow_rechunk=True,
                                                  output_sizes=dict(parameter=3 + 2 * ncov, parameter2=3 + 2 * ncov)))
        my_logger.debug('Using dask')
    else:
        extra_args = {}
    if weights is not None:
        gev_args += [weights]
        input_core_dims += [dim]
        kwargs.update(weights=True)
        my_logger.debug('Using weights')
    if initial is not None:
        gev_args += [initial.rename(parameter="initial_parameter")]
        input_core_dims += [['initial_parameter']]
        kwargs.update(initial=True,initial_names=initial.parameter.values)
        extra_args.update(on_missing_core_dim='copy')
        my_logger.debug('Using initial values')

    my_logger.debug('Doing fit')
    kwargs.update(cov_names=cov_names) # pass in the names of the covariances.
    params, std_err, cov_param, nll, AIC, ks, param_names = xarray.apply_ufunc(gev_fit_wrapper, *gev_args,
                                                                               input_core_dims=input_core_dims,
                                                                               output_core_dims=output_core_dims,
                                                                               vectorize=True, kwargs=kwargs,
                                                                               **extra_args
                                                                               )
    my_logger.debug('Done fit. Making dataset')
    # extract the names of the parameters. Mild pain as they are returned from the scaler fit function.
    param_names = param_names.stack(fdim=dim_not_collapsed).isel(fdim=0).values.tolist()

    # name variables and then combine into one dataset.
    params = params.rename("Parameters")
    std_err = std_err.rename("StdErr")
    cov_param = cov_param.rename("Cov")
    nll = nll.rename('nll').squeeze()
    AIC = AIC.rename('AIC').squeeze()
    ks = ks.rename('ks').squeeze()
    data_array = xarray.Dataset(dict(Parameters=params, StdErr=std_err, Cov=cov_param, nll=nll, AIC=AIC, KS=ks)
                                ).assign_coords(
        parameter=param_names, parameter2=param_names)
    if use_dask:
        my_logger.debug('Computing for dask')
        logger = logging.getLogger('distributed.utils_perf')
        logger.setLevel('ERROR')
        data_array = data_array.compute()
        my_logger.debug('Done dask GEV computation')

    if name:
        data_array.attrs.update(name=name)
    if extra_attrs:
        data_array.attrs.update(extra_attrs)
    if file is not None:
        file.parent.mkdir(exist_ok=True, parents=True)  # make directory
        data_array.to_netcdf(file)  # save the dataset.
        my_logger.info(f"Wrote fit information to {file}")
    return data_array


## use apply ufunc to generate distributions...

def fn_isf(c, loc, scale, p: typing.Optional[np.ndarray] = None, dist=scipy.stats.genextreme):
    shape = list(c.shape) + list(p.shape)
    fd = dist(np.broadcast_to(c, shape), loc=np.broadcast_to(loc, shape), scale=np.broadcast_to(scale, shape))
    # handle single p value.
    # if len(p) == 1:
    #    p=p.reshape(1,-1)

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
