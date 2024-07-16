# mini python script to test if can make gev fit work with dask.
from R_python import gev_r
import ausLib
import xarray
import multiprocessing
if __name__ == "__main__":
    multiprocessing.freeze_support()  # needed for obscure reasons I don't get!
    my_logger=ausLib.setup_log(1)
    c=ausLib.dask_client() # have dask running.
    ds=xarray.tutorial.open_dataset("air_temperature",chunks=dict(time=100))
    my_logger.info('Loading dataset')
    air = ds.air.load()
    my_logger.info('Dataset loaded')
    my_logger.info('Doing GEV fit using python')
    fitp=gev_r.xarray_gev_python(ds.air.chunk(lat=5,lon=5,time=-1),dim='time',use_dask=True) # gev fit No compu
    my_logger.info('Computing')
    fitp=fitp.compute()
    my_logger.info('Doign R GEV fit in memory')
    fitr=gev_r.xarray_gev(air,dim='time',use_dask=False) # gev fit No compu
    my_logger.info('Doing GEV fit using R')
    fitr=gev_r.xarray_gev(ds.air.chunk(lat=5,lon=5,time=-1),dim='time',use_dask=True) # gev fit No compu
    my_logger.info('Computing ')
    fitr=fitr.compute()
    my_logger.info('Done')