from unittest import TestCase
from models.toolbox.geometry import Geometry
import numpy as np


class TestGeometry(TestCase):

    def test_stratified_angle_approximate(self):

        cases = [
            {"alphaL": 0 , "betha": 0},
            {"alphaL": 0.25 , "betha": 2.30988146},
            {"alphaL": 0.5, "betha": np.pi},
            {"alphaL": 0.75 , "betha": 3.97330385},
            {"alphaL": 1, "betha": 2*np.pi},
        ]

        for case in cases:
            betha = Geometry.stratified_angle_approximate(case["alphaL"])
            self.assertAlmostEqual(betha, case["betha"], places=1)

    def test_stratified_angle_approximate_array(self):

        alphaL = np.array([0,0.25,0.5,0.75,1])
        betha = np.array([0,2.30988146,np.pi,3.97330385, 2*np.pi])

        betha_calc = Geometry.stratified_angle_approximate(alphaL)
        np.testing.assert_almost_equal(betha_calc, betha, decimal=1)

    def test_stratified_angle(self):

        cases = [
            {"alphaL": 0 , "betha": 0},
            {"alphaL": 0.25 , "betha": 2.30988146},
            {"alphaL": 0.5, "betha": np.pi},
            {"alphaL": 0.75 , "betha": 3.97330385},
            {"alphaL": 1, "betha": 2*np.pi},
        ]

        for case in cases:
            betha = Geometry.stratified_angle(case["alphaL"])
            self.assertAlmostEqual(betha, case["betha"], places=7)

    def test_stratified_angle_array(self):

        alphaL = np.array([0,0.25,0.5,0.75,1])
        betha = np.array([0,2.30988146,np.pi,3.97330385, 2*np.pi])

        betha_calc = Geometry.stratified_angle(alphaL)
        np.testing.assert_almost_equal(betha_calc, betha, decimal=1)
