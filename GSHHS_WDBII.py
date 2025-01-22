"""
Support for obtaining and reading the GSHHSS coastlines and WDBII river data.
See https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
Comes with methods to provide coastlines, lakes and rivers.
Area of features and class of rivers depends on the resolution.
GSHHS coastlines/lakelines and rivers seems slightly offset from one-another.
"""

import cartopy
import cartopy.crs
import cartopy.io
import pathlib
import cartopy.feature
import cartopy.io.shapereader
import fiona
import itertools
import functools


# use existing functionality for location of GSHHS (and WDBII) zip file


class GSHHS_WDBII(cartopy.io.Downloader):
    """
    Class for dealing with GSHHS/WBII data.
    init method will attempt  to download data (if not already there).
    """
    #GSHHS_URL = cartopy.io.shapereader.GSHHSShpDownloader._GSHHS_URL_TEMPLATE  # Note this is not a public method so could change.
    GSHHS_URL = 'https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/latest/gshhg-shp-2.3.7.zip'
    template_path = pathlib.Path(cartopy.config['data_dir']) / 'GSHHS_WDBII' / pathlib.Path(GSHHS_URL).name
    allowed_scales = ('f', 'h', 'i', 'l', 'c')
    translate = dict(full='f', high='h', intermediate='i', low='l', crude='c')  # synonyms for the scales.
    default_scale = 'h'  # if nothing specified use this scale.
    # default plot kws for different types of data.
    default_plot_kws = dict(coastline=dict(edgecolor='black', facecolor='none'),
                            river=dict(facecolor='none', edgecolor='blue'),
                            intermittent_river=dict(facecolor='none', edgecolor='blue', linestyle='dashed'),
                            canal=dict(facecolor='none', edgecolor='blue', linestyle='dotted'),
                            lake=dict(facecolor='cornflowerblue', edgecolor='cornflowerblue')
                            )

    def __init__(self,
                 url_template:str=GSHHS_URL,
                 target_path_template:str=str(template_path),
                 pre_downloaded_path_template:str=''):
        """

        :param url_template: template for URL -- really a string as independent of anything else.
        :param target_path_template: path as str to where file will be downloaded to.
        :param pre_downloaded_path_template: Not really used but see superclass for documentation
        """
        super().__init__(url_template, target_path_template,
                         pre_downloaded_path_template)

    def check_trans_scale(self, scale):
        """
        Check that scale is OK and translate scales like "full" to f. Uses allowed_scales & translate
        If scale, after translation, is not in allowed scale then ValueError will be raised.
        :param scale: scale. If None the self.default_scale will be returned.
        :return: translated scale
        """

        if scale is None:
            rtn_scale = self.default_scale
        else:
            rtn_scale = self.translate.get(scale, scale)
        if rtn_scale not in self.allowed_scales:
            raise ValueError(f"Scale {rtn_scale} not in allowed scales {self.allowed_scales}")
        return rtn_scale

    def path(self, scale='h', level=1, river=False):
        """
        Compute path for reading appropriate GSHHG or WDBII resource.
        :param scale: scale used. Should be one of:
            f: full resoln
            h: high resoln
            i: intermediate resoln
            l: low resoln
            c: crude resoln
            See   http://www.soest.hawaii.edu/pwessel/gshhg/index.html
            There is not much point using intermediate to crude resoln as natural_earth
              covers those scales in a more consistent way.
        :param level: Level wanted.
         Interpretation depends if using "shoreline" data or WDBII river data. See river parameter.
         Shoreline data:
            1: boundary between land and ocean, except Antarctica.
            2: boundary between lake and land.
            3: boundary between island-in-lake and lake.
            4: boundary between pond-in-island and island.
            5: boundary between Antarctica ice and ocean.
            6: boundary between Antarctica grounding-line and ocean.
        River data:
            1: Double-lined rivers (river-lakes).
            2: Permanent major rivers.
            3: Additional major rivers.
            4: Additional rivers.
            5: Minor rivers.
            6: Intermittent rivers - major.
            7: Intermittent rivers - additional.
            8: Intermittent rivers - minor.
            9: Major canals.
            10: Minor canals.
            11: Irrigation canals.

        :param river: if True return path to river data rather than shoreline data.
        :return: path as string suitable for feeding to fiona.load or cartopy.io.shapereader.Reader
        """
        scale = self.check_trans_scale(scale)
        # using pathlib module so don't have to worry about seperators.
        pth = pathlib.Path(super().path(dict()))# use the super-path. No translation going on so pass in an empty dir.
        if river:
            pth = pth  / f'WDBII_shp/{scale}/WDBII_river_{scale}_L{level:02d}.shp'
        else:
            pth = pth / f'GSHHS_shp/{scale}/GSHHS_{scale}_L{level:1d}.shp'

        pth = r'zip://' + pth.as_posix()  # so fiona recognizes it is in a zip file.
        # TODO: check pth exists and raise error??
        return pth
    @functools.lru_cache()
    def feature(self, scale='h', level=1, river=False, min_area=None, default_name=None, **kwargs):
        """
        Return a shapely feature from GSHSS/WDBII data

        :param scale: scale to be used. Default is 'h'
        :param level: level to be used.
        :param river: True then return "river". If False return "shoreline" data.
        See path method for details of scale & level
        :param  min_area: If not None then min area in decimal degrees**2 for feature to be kept.
        :param default_name: name in self.default_plot_kws to use to set up default values for Shapely Feature.
        Remaining kwargs are passed into  cartopy.feature.ShapelyFeature
        :return: a ShapelyFeature.
        """
        filename = self.path(scale=scale, level=level, river=river)
        shape = cartopy.io.shapereader.Reader(filename)
        if default_name is not None:  # have a default name
            # lets update kwargs if they key exists in the default_plot_kws
            for k, v in self.default_plot_kws[default_name].items():
                if kwargs.get(k) is None:
                    kwargs[k] = v

        if min_area is not None:  # filter by min_area
            geometries = [g for g in shape.geometries() if g.area > min_area]
        else:
            geometries = shape.geometries()
        # create the feature
        feature = cartopy.feature.ShapelyFeature(geometries, crs=cartopy.crs.PlateCarree(), **kwargs)
        return feature

    def min_area_scale(self, scale):
        """
        Returns min_area in degree^2 for a feature depending on the resoln of the dataset.
        Scale  min_area
           i     0.01
           l     0.05
           c     0.1
        :param scale -- scale wanted
        :return min_area
        """
        min_area = None
        if scale == 'i':
            min_area = 0.01
        elif scale == 'l':
            min_area = 0.05
        elif scale == 'c':
            min_area = 0.1

        return min_area

    def coastlines(self, scale=None, min_area=None, **kwargs):
        """
        Generate GSHHS coastline.
            Union of land/sea coastline (level=1) + Antarctica ice/ocean line (level=5) (circa 2010?)
        :param scale: scale wanted
        :param min_area: minimum area wanted to retain feature. If None then will be computed using min_area_scale
                    If you want all features regardless of area set min_area to 0.0
        : all other kwargs are passed to ShapelyFeature and used there.
        : Default colors (and other attributes) are set in default_plot_kws['coastline']
        :return:cartopy feature of GSHSS coastlines
        """
        scale = self.check_trans_scale(scale)
        for k, v in self.default_plot_kws['coastline'].items():
            if kwargs.get(k) is None:
                kwargs[k] = v

        if min_area is None:  # compute the min_area.
            min_area = self.min_area_scale(scale)

        coastline = self.feature(scale=scale, level=1, min_area=min_area)
        antarctica = self.feature(scale=scale, level=5, min_area=min_area)
        geoms = itertools.chain(coastline.geometries(), antarctica.geometries())

        coast = cartopy.feature.ShapelyFeature(geoms, crs=cartopy.crs.PlateCarree(), **kwargs)
        return coast

    def lakes(self, scale=None, min_area=None, **kwargs):
        """
        Generate lake feature.

        :param scale: scale wanted
        :param min_area -- if set then this is the min_area is degrees**2  for lake to be shown.
            If not set then value will be computed using min_area_scale.
            If you want all features regardless of area set min_area to 0.0
        :param kwargs passed to ShapelyFeature and used there.
                : Default colors (and other attributes) are set in default_plot_kws['lake']
        :return:cartopy feature of GSHSS lakes
        """
        scale = self.check_trans_scale(scale)
        if min_area is None:  # compute the min_area.
            min_area = self.min_area_scale(scale)

        lake = self.feature(scale=scale, level=2, min_area=min_area, default_name='lake', **kwargs)
        return lake

    def water_feature(self, scale, levels, default_name=None, **kwargs):
        """
        Return single cartopy shapely feature of desired water features
        :param scale: scale wanted
        :param levels: iterable of levels wanted
        :param default_name: name used to augment kwargs
        :param kwargs: kwargs passed through to ShapelyFeature.
        :return: ShapelyFeature
        """
        scale = self.check_trans_scale(scale)
        if default_name is not None:
            for k, v in self.default_plot_kws[default_name].items():
                if kwargs.get(k) is None:
                    kwargs[k] = v
        geoms = []
        for level in levels:  # loop over levels getting the geometries
            geoms.append(self.feature(scale=scale, level=level, river=True).geometries())

        geomsi = itertools.chain(*geoms)  # convert into one chained iterable
        waterf = cartopy.feature.ShapelyFeature(geomsi, crs=cartopy.crs.PlateCarree(), **kwargs)
        return waterf

    def rivers(self, scale=None, levels=None, **kwargs):
        """
        Return scale relevant rivers
        :param scale: Scale wanted:
            if 'f' or 'h' then all rivers will be returned. (levels 1 to 5)
            if 'i' then rivers levels 1 to 4 will be returned
            if 'l' then rivers levels 1 to 3 returned
            if 'c' then rivers levels 1 to 2 returned
        :param levels: if not None an iterable of the levels wanted which over rules the scale decision
        : Default colors (and other attributes) are set in default_plot_kws['river']
        :param kwargs passed to ShapelyFeature and used there.
        :return:cartopy feature of WDBII rivers
        """
        scale = self.check_trans_scale(scale)
        if scale not in self.allowed_scales:
            raise ValueError(f"scale {scale} not in {self.allowed_scales}")
        if levels is None:
            if scale in ['f', 'h']:
                levels = [1, 2, 3, 4, 5]
            elif scale == 'i':
                levels = [1, 2, 3, 4]
            elif scale == 'l':
                levels = [1, 2, 3]
            elif scale == 'c':
                levels = [1, 2]
            else:
                pass

        rivers = self.water_feature(scale, levels, default_name='river', **kwargs)
        return rivers

    def inter_rivers(self, scale=None, levels=None, **kwargs):
        """
        Return scale relevant intermittent rivers

        :param scale: Scale wanted:
            if 'f' or 'h' then all intermittent rivers will be returned. (levels 6 to 8)
            if 'i' then intermittent rivers levels 6 to 7 will be returned
            if 'l' then intermittent rivers levels 6 returned
            if 'c' then no intermittent rivers returned
        :param levels: If not None an iterable of the levels wanted which over rules the scale decision
        : Default colors (and other attributes) are set in default_plot_kws['intermittent_river']
        :param kwargs passed to ShapelyFeature and used there.
        :return:cartopy feature of WDBII rivers
        """
        scale = self.check_trans_scale(scale)

        if levels is None:
            if scale in ['f', 'h']:
                levels = [6, 7, 8]
            elif scale == 'i':
                levels = [6, 7]
            elif scale == 'l':
                levels = [6]
            elif scale == 'c':
                # return empty feature.
                return cartopy.feature.ShapelyFeature((), crs=cartopy.crs.PlateCarree())
            else:
                pass

        inter_rivers = self.water_feature(scale, levels, default_name='intermittent_river', **kwargs)
        return inter_rivers

    def canals(self, scale=None, levels=None, **kwargs):
        """
        Return scale relevant  canals
        :param scale: Scale wanted:
            if 'f' or 'h' then all canals will be returned. (levels 9 to 11)
            if 'i' then canals levels 9 to 10 will be returned
            if 'l' then canals levels 9 returned
            if 'c' then no canals returned (an empty ShapelyFeature)
        :param levels (default None). If not None an iterable of the levels wanted which over rules the scale decision
        : Default colors (and other attributes) are set in default_plot_kws['canal']
        :param kwargs passed to ShapelyFeature and used there.
        :return:cartopy feature of WDBII rivers
        """

        scale = self.check_trans_scale(scale)
        if levels is None:  # need levels.
            if scale in ['f', 'h']:
                levels = [9, 10, 11]
            elif scale == 'i':
                levels = [9, 10]
            elif scale == 'l':
                levels = [9]
            elif scale == 'c':  # nothing to be returned.
                return cartopy.feature.ShapelyFeature((), crs=cartopy.crs.PlateCarree(), **kwargs)
            else:
                pass

        inter_rivers = self.water_feature(scale, levels, default_name='canal', **kwargs)
        return inter_rivers
