from PyMUSE.ds9 import Ds9RegionString
import unittest

import numpy as np

class TestEllipse(unittest.TestCase):
    def test_ellipse_params_to_region_string_returns_string(self):
        x_c, y_c, a, b = np.random.randint(0, 100, 4)
        theta = np.random.random()
        region_string = Ds9RegionString.ellipse_params_to_region_string(x_c, y_c, a, b, theta)
        self.assertIsInstance(region_string, str)





if __name__ == '__main__':
    unittest.main()
