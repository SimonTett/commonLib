import unittest

import pandas as pd
import xarray

from R_python import gev_r
from scipy import stats
import numpy as np
import pandas.testing as pd_tests
import xarray.testing as x_tests
class test_R_gev(unittest.TestCase):
    def test_gev_fit_new(self):
        data = stats.genextreme.rvs(0.1, loc=10, scale=5, size=10000)
        expected = pd.Series([0.1, 10, 5],index=['shape','location','scale']).rename('par')
        pars,se,cov,nll,aic,ks = gev_r.gev_fit_new([data], verbose=False)
        self.assertEqual(len(pars), 3)
        self.assertEqual(len(se), 3)
        self.assertEqual(cov.shape, (3,3))
        pd_tests.assert_series_equal(pars.round(1), expected, atol=1e-1,check_like=True)
        # test weights don't change the result
        wts = np.ones(len(data))
        wts[::2] = 2.0
        pars2,se2,cov2,nll2,aic2,ks2 = gev_r.gev_fit_new([data, wts],weights=True, verbose=False)
        pd_tests.assert_series_equal(pars.round(1), expected, atol=1e-1, check_like=True)
        self.assertTrue(np.any(pars2 != pars)) # values should be different
        self.assertTrue(np.isnan(ks2)) # ks should be nan
        # test with covariates. Will increase the mean with time
        time  = np.arange(0,len(data))/100.
        data2 = stats.genextreme.rvs(0.1, loc=10, scale=5, size=10000) *(1+time*0.1) # expect to get a change in location and scale
        expected = pd.Series([0.1, 10, 5,1,0.5],index=['shape','location','scale','Dlocation_cov1','Dscale_cov1']).rename('par')
        pars3,se3,cov3,nll3,aic3,ks3 = gev_r.gev_fit_new([data2,time], verbose=True,initial=dict(location=[15,0.9],scale=[6.0,0.4],shape=0.01))
        pd_tests.assert_series_equal(pars3.round(1), expected, rtol=1e-1, check_like=True)

    def test_gev_xarray(self):
        # test case with xarray. Will have a two-dimensional array with time and space
        time = np.arange(0,10000)/100.
        space = np.arange(0,2)
        data = np.stack([stats.genextreme.rvs(0.5, loc=10*(s+1), scale=5, size=10000) * (1 + time * 0.1) for s in space],axis=1)
        ds = xarray.DataArray(data, dims=['time','space'], coords={'time':time,'space':space})
        expect_data = np.stack([np.array([0.5,10.0*(s+1),5.0,s+1.,0.5]) for s in space],axis=1)
        expect = xarray.DataArray(expect_data,coords=dict(parameter=['shape','location','scale','Dlocation_time','Dscale_time'],space=space))
        cov = xarray.DataArray(time, name='time')
        fit = gev_r.xarray_gev(ds, ds.time,dim='time')
        self.assertEqual(fit.sizes, dict(parameter=5,parameter2=5,space=2))
        expect = expect.reindex_like(fit['Parameters']).T # make shape like Parameters
        x_tests.assert_allclose(fit['Parameters'],expect,rtol=1e-1)
        # try with initial choices. Same choice for everything
        fit = gev_r.xarray_gev(ds, ds.time, dim='time',initial=dict(shape=0.1,location=[10,1.0],scale=[5.0,0.5]))
        x_tests.assert_allclose(fit['Parameters'], expect, rtol=1e-1)





if __name__ == '__main__':
    unittest.main()
