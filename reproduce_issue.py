from cosipy import UnBinnedData
from cosipy.spacecraftfile import SpacecraftFile
from cosipy.polarization.conventions import MEGAlibRelativeX, MEGAlibRelativeY, MEGAlibRelativeZ, IAUPolarizationConvention
from cosipy.polarization.polarization_asad import PolarizationASAD
from cosipy.threeml.custom_functions import Band_Eflux
from astropy.time import Time
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u

grb_background = UnBinnedData('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb.yaml')
#grb_background.select_data_time(unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background.fits.gz', output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background_source_interval')
#grb_background.select_data_energy(200., 10000., output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background_source_interval_energy_cut', unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background_source_interval.fits.gz')
data = grb_background.get_dict_from_fits('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background_source_interval_energy_cut.fits.gz')

background_before = UnBinnedData('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_before.yaml')
#background_before.select_data_time(unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background.fits.gz', output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_before') 
#background_before.select_data_energy(200., 10000., output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_before_energy_cut', unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_before.fits.gz')
background_1 = background_before.get_dict_from_fits('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_before_energy_cut.fits.gz') 

background_after = UnBinnedData('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_after.yaml')
#background_after.select_data_time(unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/grb_background.fits.gz', output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_after')
#background_after.select_data_energy(200., 10000., output_name='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_after_energy_cut', unbinned_data='/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_after.fits.gz')
background_2 = background_after.get_dict_from_fits('/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/tutorial/background_after_energy_cut.fits.gz')

background = [background_1, background_2]

response_file = '/Users/eneights/work/COSI/data-challenges/data_challenge_3/polarization_asad/psr/ResponseContinuum.o3.pol.e200_10000.b4.p12.s10396905069491.m441.filtered.nonsparse.binnedpolarization.11D_nside8.area.h5'

sc_orientation = SpacecraftFile.parse_from_file('/Users/eneights/work/COSI/data-challenges/data_challenge_3/files/DC3_SAA.ori')
sc_orientation = sc_orientation.source_interval(Time(1835493492.2, format = 'unix'), Time(1835493492.8, format = 'unix'))

source_direction = SkyCoord(l=23.53, b=-53.44, frame='galactic', unit=u.deg)

a = 100. * u.keV
b = 10000. * u.keV
alpha = -0.7368949
beta = -2.095031
ebreak = 622.389 * u.keV
K = 300. / u.cm / u.cm / u.s

spectrum = Band_Eflux(a = a.value,
                      b = b.value,
                      alpha = alpha,
                      beta = beta,
                      E0 = ebreak.value,
                      K = K.value)

spectrum.a.unit = a.unit
spectrum.b.unit = b.unit
spectrum.E0.unit = ebreak.unit
spectrum.K.unit = K.unit

asad_bin_edges = Angle(np.linspace(-np.pi, np.pi, 17), unit=u.rad)

grb_polarization = PolarizationASAD(source_direction, spectrum, asad_bin_edges, data, background, sc_orientation, response_file, response_convention='RelativeZ', show_plots=False)