import astropy.units as u

class Image:
    def __init__(self, filename=None, pixelsize=0.2*u.arcsec,  flux_units=1E-20 * u.erg / u.s / u.cm ** 2 / u.angstrom, data=None, var=None):
        self.flux_units = flux_units
        self.pixelsize=pixelsize

        self.filename = filename

        self.data = data
        self.var = var

        if filename:
            self.load_data()



    def load_data(self):
        pass






