import unittest
from PyMUSE.musecube_new import MuseCube
from PyMUSE.spectrum import Spectrum
import numpy as np


class TestGetItem(unittest.TestCase):
    def setUp(self):
        self.filename = "../minicube.fits"
        self.cube = MuseCube(self.filename)
        self.spectrum = self.cube[:, 0, 0]

    def test_get_item_slice_return_a_value(self):
        value = self.spectrum[0]
        self.assertIsInstance(value, np.float)

    def test_slice_return_a_spectrum(self):
        # corregir en la clase
        subspectrum = self.spectrum[0:3]
        self.assertIsInstance(subspectrum, Spectrum)


class TestIter(unittest.TestCase):
    def setUp(self):
        self.filename = "../minicube.fits"
        self.cube = MuseCube(self.filename)
        self.spectrum = self.cube[:, 0, 0]


class TestAritmetic(unittest.TestCase):
    def setUp(self):
        self.filename = "../minicube.fits"
        self.cube = MuseCube(self.filename)

    def test_add_two_spectra_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0, 0]
        spectrum2 = self.cube[:, 0, 0]
        sum = spectrum1 + spectrum2
        self.assertIsInstance(sum, Spectrum)

    def test_add_spectrum_and_number_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0 ,0]
        sum = spectrum1 + 1
        self.assertIsInstance(sum, Spectrum)

    def test_sub_two_spectra_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0, 0]
        spectrum2 = self.cube[:, 0, 0]
        sub = spectrum1 - spectrum2
        self.assertIsInstance(sub, Spectrum)

    def test_sub_spectrum_and_number_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0, 0]
        sub = spectrum1 - 1
        self.assertIsInstance(sub, Spectrum)

    def test_mul_two_spectra_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0, 0]
        spectrum2 = self.cube[:, 0, 0]
        mul = spectrum1 * spectrum2
        self.assertIsInstance(mul, Spectrum)

    def test_mul_spectrum_and_number_returns_Spectrum(self):
        spectrum1 = self.cube[:, 0, 0]
        mul = spectrum1 * 2
        self.assertIsInstance(mul, Spectrum)


if __name__ == '__main__':
    unittest.main()