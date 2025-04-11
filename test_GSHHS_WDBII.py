import itertools
import unittest
import GSHHS_WDBII
import pathlib
import cartopy.io
import itertools


class testGSHHS_WDBII(unittest.TestCase):
    """
    Test cases for GSHSS_WDBII
    """

    def setUp(self) -> None:
        self.gshhs = GSHHS_WDBII.GSHHS_WDBII()

    def assertFeatureEqual(self, feature1, feature2, msg=None):
        """
        Assert that two cartopy shapely features are the same,
        Does by checking geometries and crs
        :param feature1: feature#1
        :param feature2: feature#2
        :param msg: mesg to be printed out if not equal.
        :return: nada
        """
        # suck up geometries making them lists.
        g1 = list(feature1.geometries())
        g2 = list(feature2.geometries())
        self.assertEqual(g1, g2 , msg)
        self.assertEqual(feature1.crs, feature2.crs, msg)

    def test_check_trans_scale(self):
        """
        test that check_trans_scale works
        :return: nada
        """
        # 1) test that std scales works.
        for scale in self.gshhs.allowed_scales:
            self.assertEqual(scale, self.gshhs.check_trans_scale(scale))
        # 2) test translation works
        for scale in ['full', 'high', 'intermediate', 'low', 'crude']:
            self.assertEqual(scale[0], self.gshhs.check_trans_scale(scale))
        # 3) test bad names causes ValueError
        with self.assertRaises(ValueError):
            self.gshhs.check_trans_scale('shouldFail')

    def test_path(self):
        """
        test path method
        :return: nada
        """

        GSHHS_URL = cartopy.io.shapereader.GSHHSShpDownloader._GSHHS_URL_TEMPLATE  # Note this is not a public method so could change.
        root_name = pathlib.Path(cartopy.config['data_dir']) / 'GSHHS_WDBII' / pathlib.Path(GSHHS_URL).name

        # 1) Crude GSHHS coastline level 1
        pth = self.gshhs.path(scale='crude', level=1)
        expect_name = 'zip://' + (root_name / 'GSHHS_shp/c/GSHHS_c_L1.shp').as_posix()
        self.assertEqual(expect_name, pth)

        # 2) Crude river level 1
        pth = self.gshhs.path(scale='crude', level=1, river=True)
        expect_name = 'zip://' + (root_name / 'WDBII_shp/c/WDBII_river_c_L01.shp').as_posix()
        self.assertEqual(expect_name, pth)

    def test_feature(self):
        """
        Test feature works
        :return:
        """
        for scale in ['c']:
            for level in [1, 2]:
                feature = self.gshhs.feature(scale=scale, level=level, min_area=0.0)
                path = self.gshhs.path(scale=scale, level=level)
                shape = cartopy.io.shapereader.Reader(path)
                expect = cartopy.feature.ShapelyFeature(shape.geometries(),
                                                        crs=cartopy.crs.PlateCarree())
                self.assertFeatureEqual(feature, expect)

    def test_min_area_scale(self):
        """
        test min_area_scale
        :return:
        """

        expect = dict(f=None, h=None, i=0.01, l=0.05, c=0.1)
        for scale, ma in expect.items():
            got = self.gshhs.min_area_scale(scale)
            self.assertEqual(got, ma)

    def test_coastlines(self):
        """
        Test coastline
        :return:
        """
        # test case -- read in L1 & l5 C level, with no filtering, merge them and compare them.
        coastline = self.gshhs.coastlines(scale='c', min_area=0.0)
        geoms = []
        for l in [1, 5]:
            feature = self.gshhs.feature(scale='c', level=l)
            geoms.append(feature.geometries())
        geomsi = itertools.chain(*geoms)
        expect = cartopy.feature.ShapelyFeature(geomsi, crs=cartopy.crs.PlateCarree())
        self.assertFeatureEqual(coastline, expect)

    def test_lakes(self):
        """
        Test lakes are as expected
        :return: nada
        """
        lakes = self.gshhs.lakes(scale='c', min_area=0.0)
        expect = self.gshhs.feature(scale='c', level=2)
        self.assertFeatureEqual(lakes, expect)

    def test_water_feature(self):
        """
          test water_feature works
          :return: nada
        """

        waterf = self.gshhs.water_feature(scale='c', levels=[1, 2])
        geoms = []
        for l in [1, 2]:
            feature = self.gshhs.feature(scale='c', level=l, river=True)
            geoms.append(feature.geometries())
        geomsi = itertools.chain(*geoms)
        expect = cartopy.feature.ShapelyFeature(geomsi, crs=cartopy.crs.PlateCarree())
        self.assertFeatureEqual(waterf, expect)

    def test_rivers(self):
        """
        Test that rivers works.
        :return: Nothing
        """

        # crude resoln => levels 1 & 2
        expect = self.gshhs.water_feature(scale='crude', levels=[1, 2])
        river = self.gshhs.rivers(scale='crude')
        self.assertFeatureEqual(expect, river)

    def test_inter_rivers(self):
        """
        Test intermittent river method
        :return: Nothing
        """
        # Intermittent rivers don't show at scale=c so use low resoln instead
        inter_river = self.gshhs.inter_rivers(scale='low')
        expect = self.gshhs.water_feature(scale='low',levels=[6])
        self.assertFeatureEqual(inter_river,expect)

    def test_canals(self):
        """
        test canals method
        :return:
        """

        # canals don't show at scale=c so use low resoln instead
        canal = self.gshhs.canals(scale='low')
        expect = self.gshhs.water_feature(scale='low',levels=[9])
        self.assertFeatureEqual(canal,expect)
