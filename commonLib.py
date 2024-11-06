"""
Library of functions/methods that are worth reusing and so shared across multiple
projects.
"""
import datetime
import pathlib
import cftime
import pandas as pd
import requests
import io
import platform
import shutil
import functools  # for cache
import xarray
import matplotlib.pyplot as plt
import typing
import sys
import logging

my_logger = logging.getLogger(__name__)


def init_log(
        log: logging.Logger,
        level: typing.Union[typing.Literal['DEBUG', 'INFO', 'WARNING', 'ERROR'], int] = 'WARNING',
        log_file: typing.Optional[typing.Union[pathlib.Path, str]] = None,
        mode: str = 'a'
) -> None:
    """
    Set up logging on a logger! Will clear any existing logging
    :param log: logger to be changed
    :param level: level to be set.
    :param log_file:  if provided pathlib.Path to log to file
    :param mode: mode to open log file with (a  -- append or w -- write). Only used in log_file provided.
    :return: nothing -- existing log is modified.
    """
    log.handlers.clear()
    log.setLevel(level)
    formatter = logging.Formatter('%(asctime)s %(levelname)s:  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler(sys.stderr)
    ch.setFormatter(formatter)
    log.addHandler(ch)
    # add a file handler.
    if log_file:
        if isinstance(log_file, str):
            log_file = pathlib.Path(log_file)
        log_file.parent.mkdir(exist_ok=True, parents=True)
        fh = logging.FileHandler(log_file, mode=mode + 't')  #
        fh.setLevel(level)
        fh.setFormatter(formatter)
        log.addHandler(fh)
    log.propagate = False


def high_resoln():
    """
    Set up matplotlib for high resolution screen where want to make things about 75% larger.
    I think needed because windows default dpi means 100% scaling.
    Will close all existing figures.
    :return: nada
    """
    plt.close('all')
    plt.matplotlib.rcParams['figure.dpi'] = 175  # make figures 75% bigger.


def read_cet(file=None, retrieve=False, direct='data', mean='seasonal', temp_type='mean'):
    """

    :param file: name of file to use. If None then file name will be constructed
    :param retrieve: If true (default is False) retrieve data. Otherwise, read local data.
       If local data does not exist then data will be retrieved.

    :param direct: name of directory where data to be retrieved to or read from.
    :param mean: What mean to be retrieved. (seasonal|monthly|daily)
    :param temp_type: What temp_type of CET to be retrieved (mean|max|min)
    :return: xarray
    """
    mo_cet_root = 'https://www.metoffice.gov.uk/hadobs/hadcet'
    mo_cet_root = 'https://www.metoffice.gov.uk/hadobs/hadcet/data/'
    # ALL but seasonal likely need updating.
    urls = dict(dailymean='meantemp_daily_totals.txt', monthlymean='meantemp_monthly_totals.txt',
                seasonalmean='meantemp_seasonal_totals.txt', dailymin='cetmindly1878on_urbadj4.dat',
                monthlymin='cetminmly1878on_urbadj4.dat', seasonalmin='sn_HadCET_min.txt',
                dailymax='cetmaxdly1878on_urbadj4.dat', monthlymax='cetmaxmly1878on_urbadj4.dat',
                seasonalmax='sn_HadCET_max.txt', )

    nskip = dict(monthly=4, seasonal=9, daily=0)
    month_lookups = dict(JAN=1, FEB=2, MAR=3, APR=4, MAY=5, JUN=6, JUL=7, AUG=8, SEP=9, OCT=10, NOV=11, DEC=12, DJF=1,
                         MAM=4, JJA=7, SON=10, WIN=1, SPR=4, SUM=7, AUT=10
                         )  # month

    if file is None:
        file = f"cet_{mean}_{temp_type}.nc"
    path = pathlib.Path(direct) / file
    if (not path.exists()) or retrieve:
        # retrieve data from MO.
        url = mo_cet_root + '/' + urls[mean + temp_type]
        # will trigger an error if mean or temp_type not as expected though error won't be very helpful...
        print(f"Retrieving data from {url}")
        headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:77.0) Gecko/20190101 Firefox/77.0'}  # fake we are interactive...
        r = requests.get(url, headers=headers)
        rdata = io.StringIO(r.text)
        data = pd.read_csv(rdata, skiprows=nskip.get(mean, 0), header=[0], sep=r'\s+', na_values=[-99.9])
        # need to use cftime to make time-coords... so all a bit of a pain! and will be doubly so for daily data..
        if mean != 'daily':
            dates = []
            values = []
            for c in data.columns:
                month = month_lookups.get(c.upper())
                if month is None:
                    continue
                if temp_type == 'daily':
                    raise Exception(f"Can't handle {temp_type} data")
                else:
                    dates.extend([cftime.datetime(yr, month, 1, calendar='proleptic_gregorian') for yr in data.Year])
                    values.extend(data.loc[:, c].values)
        else:
            # Function to convert string to cftime datetime. Thanks chatgpt co-pilot
            def to_cftime(date_str):
                return cftime.datetime.strptime(date_str, '%Y-%m-%d', calendar='gregorian')
            dates = data.Date.apply(to_cftime)
            values = data.Value

        ts = xarray.DataArray(values, coords=dict(time=dates)).rename(f'CET{mean}{temp_type}').sortby('time')
        pathlib.Path(direct).mkdir(parents=True, exist_ok=True)  # make (if needed directory to put the data
        ts.to_netcdf(path)  # write out the data.
    else:
        # just retrieve the cached data.
        ts = xarray.load_dataarray(path)

    return ts


def saveFig(
        fig,
        name: typing.Optional[str] = None,
        savedir: typing.Optional[pathlib.Path] = None,
        figtype: typing.Union[str, typing.List[str]] = 'png', dpi: int = 300,
        transpose: bool = False
):
    """

    :param fig -- figure to save
    :param name (optional) set to None if undefined use fig.get_label
    :param savedir (optional) directory as a pathlib. Path to save figure to. Default is figures
    :param figtype graphics type to save as. If list then types to save as
    :param dpi: resolution to save at.    else:
    :param transpose: If True then generate a transposed version of the figure before saving.
    Transposition will be reversed at end.
        This allows figures created for papers to go to landscape.
        If no name provided name will be fig.get_label+'_transpose'
    """
    if transpose:
        fig.set_size_inches(fig.get_size_inches()[::-1]) # transpose the figure #
    if dpi is None:
        dpi = 300
    # set up defaults
    # work out sub_plot_name.
    if name is None:
        name = fig.get_label()
        if transpose:
            name += '_transpose'
    if savedir is None:
        savedir = pathlib.Path('figures')  # always relative to where we are
    # possibly create the savedir.
    savedir.mkdir(parents=True, exist_ok=True)  # create the directory
    if isinstance(figtype, str):
        figtype = [figtype]
    for ftype in figtype:  # loop over fig types
        out_file_name = savedir / (name + "." + ftype)
        fig.savefig(out_file_name, dpi=dpi,bbox_inches='tight')
        my_logger.debug(f"Saved figure to {out_file_name} at {dpi} dpi")
    # and undo transpose
    if transpose:
        fig.set_size_inches(fig.get_size_inches()[::-1])


class plotLabel:
    """
    Class for plotting labels on sub-plots
    """

    def __init__(self, upper=False, roman=False, fontdict={}):
        """
        Make instance of plotLabel class
        parameters:
        :param upper -- labels in upper case if True
        :param roman -- labels use roman numbers if True
        """

        import string
        if roman:  # roman numerals
            strings = ['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii', 'ix', 'x', 'xi', 'xii']
        else:
            strings = [x for x in string.ascii_lowercase]

        if upper:  # upper case if requested
            strings = [x.upper() for x in strings]

        self.strings = strings[:]
        self.num = 0
        self.fontdict = fontdict

    def label_str(self):
        """
        Return the next label
        """
        string = self.strings[self.num] + " )"
        self.num += 1
        self.num = self.num % len(self.strings)
        return string

    def plot(self, ax=None, where=None):
        """
        Plot the label on the current axis.
        :param ax -- axis to plot on. Default is current axis (using plt.gca())
        :param where -- (x,y) tuple saying where  to plot label using axis coords. Default is (-0.03,1.03)
        """

        if ax is None:
            plt_axis = plt.gca()
        else:
            plt_axis = ax
        try:
            if plt_axis.size > 1:  # got more than one element
                for a in plt_axis.flatten():
                    self.plot(ax=a, where=where)
                return
        except AttributeError:
            pass

        # now go and do the actual work!

        text = self.label_str()
        if where is None:
            x = -0.03
            y = 1.03
        else:
            (x, y) = where

        plt_axis.text(x, y, text, transform=plt_axis.transAxes, horizontalalignment='right', verticalalignment='bottom',
                      fontdict=self.fontdict
                      )


# cache handling
cache_dirs = dict()


def setup_cache(config_dict):
    cache_dirs.update(config_dict)


def gen_cache_file(filepath, verbose=False):
    """
    Generate a cache file path from filepath
    :param filepath: filepath on slow file system
    :return:cached filepath.
    """

    # cache_dir =
    fileName = str(filepath)
    cache_file = filepath  # if nothing found then we just return the filepath
    for dir_name, cache_name in cache_dirs.items():
        if dir_name in fileName:
            cache_file = fileName.replace(dir_name, cache_name)
            cache_file = pathlib.Path(cache_file)
            if verbose:
                print(f"Replaced {dir_name} with {cache_name} for {fileName}")
            continue  # no need to process more
    return cache_file


@functools.cache
def cache_filename(filepath: pathlib.Path, verbose=False, use_cache=None):
    """
    Generate local cache. If the local cache does not exist then copy filepath to it.
    Returns the cache filename (which might be the same as filepath)
    Uses functools.cache so second (or subsequent) time ran in a session will just return the cached filepath
    :param use_cache: Logical -- If True use cache. If None then if on specified platform.node use cache.
    :param verbose: If True be verbose
    :param filepath: path to the file.
    :return: file read
    """
    if use_cache is None:
        use_cache = (platform.node() in ['geos-w-048'])  # list of platforms where want to cache
    if not use_cache:  # not using cache -- just return the input
        if verbose:
            print("Not using cache ")
        return filepath

    cache_file = gen_cache_file(filepath)

    if cache_file.exists() and (not filepath.exists()):  # cache exists and filepath doesn't.
        if verbose:
            print(f"Cache file: {cache_file} exists while {filepath} does not")
    elif not (cache_file.exists() and (cache_file.stat().st_mtime > filepath.stat().st_mtime)):
        # Want to use cache_file if it exists and its mtime is greater than the mtime of the file being cached
        if verbose:
            print(f"Copying data from {filepath} to {cache_file}")
        cache_file.parent.mkdir(parents=True, exist_ok=True)  # make the cache dir if neeed.
        shutil.copy(filepath, cache_file)
    if verbose:
        print(f"Cache file is {cache_file}")

    return cache_file  # return cached file.


def create_time(
        path: pathlib.Path,
        fmt: str = '%Y-%m-%d %H:%M.%S'
        ) -> str:
    """
    Return  creation time of file as a string
    :param path: path to the file
    :param fmt: format to return time in
    :return: formatted time of creation for file
    """
    return datetime.datetime.fromtimestamp(path.stat().st_ctime).strftime(fmt)
