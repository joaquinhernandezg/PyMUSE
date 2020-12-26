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
        # propagate variances pendient
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


class InverseVarianceCombiner(WeightedSpaxelCombiner):
    @staticmethod
    def get_weights(cube):
        if cube.stat is None:
            raise NotImplementedError("Variances cant be estimated empirically for the moment")
        var = cube.stat
        # collapse the variances in the spectral direction
        white_var = np.nansum(var, axis=0)
        inv_var = np.nan_to_num(1/white_var, nan=0, posinf=0, neginf=0)

        # normalizes the weights to make them sum 1
        weights = inv_var/np.sum(inv_var)
        #verify what happens when sum(inv_var)==0
        return weights

class













