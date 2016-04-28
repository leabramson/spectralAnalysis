from specfit import *
from spectrum import *
from sfhdust import *
from mcmc import *
import util, sys, os

# Define the model and priors
# Parameters in first group pertain to the spectrum as a whole (z, sigma)
# while those in the second pertain to a single stellar population (of
# which there can, in principle, be more than one).
priors = [{'z':     GaussianPrior(0.3, 0.01),
           'sigma': UniformPrior(150., 250.)},
          {'age':      LogUniformPrior(1e9, 5e9),     # yr
           'ltau':     UniformPrior(8.0, 10.0),       # log yr
           'metal':    LogUniformPrior(0.01, 0.03),   # Z
           'A_V':      UniformPrior(0.0, 1.0)}]       # mag
csp = CSP(sfh=ExponentialSFH(), dust=CalzettiDust())
model = SpectrumModel(lsf=GaussianVelocityLSF(300.), csp=[csp])

# Options
options = util.options()
options.nlive = 100           # Number of MultiNest live points.
options.poly_order = 0        # Effectively, filter continuum by polynomial of this order.
                              # Here, assume the flux calibration is reliable.
options.output_prefix = "output/fake"   # Prefix for all output files.
# Set up BC03 grid. Must use hr models. 
options.grid_path = os.environ['GRID_PATH']+"/Padova1994/salpeter/bc2003_hr_stelib_m[metal]_salp_ssp.ised"
options.grid_param = ["metal"]
options.grid_paramlabels = [["52", "62", "72"]]
options.grid_paramvalues = [[0.008, 0.02, 0.05]]
options.grid_paramlogint = [True]
options.air2vac = False       # Stay in air wavelengths.

# Define constraints
photdata = None
specdata = ['data.fits']
specerror = ['sigma.fits']
specmaskfiles = None  
# Ignoring the instrumental resolution in this example.
inst_lsf = [None]

# Make sure the output directory "output/" exists before running.

specfit(specdata, specerror, photdata, priors, model,
        options, inst_lsf, specmaskfiles=specmaskfiles)
