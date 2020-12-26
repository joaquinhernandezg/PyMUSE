from PyMUSE.musecube_new import MuseCube
from PyMUSE.spectrum import Spectrum
from abc import ABC, abstractmethod

import numpy as np

class SpaxelCombiner(ABC):
    @abstractmethod
    @staticmethod
    def extract_spectrum(self, cube):
        pass

class MedianCombiner(SpaxelCombiner):
    @staticmethod
    def extract_spectrum(self, cube):
        flux = np.median(cube.flux, axis=(1, 2))
        stat = np.median(cube.stat, axis=(1, 2))#should be changed
        return cube._to_obj(Spectrum, flux=flux, stat=stat)

class WeightedSpaxelCombiner(SpaxelCombiner):

    @abstractmethod
    @staticmethod
    def get_weights(*args, **kwargs):
        pass

    @classmethod
    def extract_spectrum(cls, cube):
        weigths = cls.get_weights(cube)
        flux = np.sum(cube.flux * weigths, axis=(1, 2))
        # propagate variances
        return cube._to_obj(Spectrum, flux=flux)

class SumCombiner(WeightedSpaxelCombiner):
    @staticmethod
    def get_weights(cube):
        return np.ones_like(cube.flux)

class MeanCombiner(WeightedSpaxelCombiner):
    @staticmethod
    def get_weights(cube):
        shape = cube.shape
        n_spaxels = shape[1]*shape[2]
        weights = np.ones_like(cube.flux)/n_spaxels
        return weights









