from specfit import *
from spectrum import *
from sfhdust import *
from mcmc import *
import util, sys, os

# Fix solar abundances, vary age and Ha EW and redshift.
# syntax is python xxxxxx.py 'objPREFIX' 'specREGION' redshift z_prior_width outdir

bin = sys.argv[2].upper() # 'INNER', 'OUTER', 'INTER'
outname = bin.lower()

bluSpec = sys.argv[1]+'_B_'+outname+'.pyspec'
redSpec = sys.argv[1]+'_R_'+outname+'.pyspec'

bluLSF = sys.argv[1]+'_B_'+outname+'.lsf'
redLSF = sys.argv[1]+'_R_'+outname+'.lsf'

#
# Set up a simple stellar population model (no SFH or dust)
#
priors = [{'z':     GaussianPrior(float(sys.argv[3]), float(sys.argv[4])),
           'HpsEW': UniformPrior(0., 100.),
           'OIIIpsEW': UniformPrior(0., 1000.)},
          {'age':   LogUniformPrior(1e7, 10e9)}]
csp = CSP(sfh=None,dust=None)
csp.param['metal'] = 0.02
# A spectral model consists of one or more composite stellar populations (i.e. with
# star formation histories and dust attenuation), with a common redshift and 
# velocity dispersion. The velocity dispersion is neglected in this case (lsf=None)
# due to the low resolution of the grism. Parameters for the whole model are 
# stored in the first dictionary in the 'priors' array (here only 'z'), and parameters
# for the CSP (only one in this case) are stored in the second (here 'age' and 'metal'=Z/H).
model = SpectrumModel(lsf=None, csp=[csp])
model.param['emsigma'] = 200.          # Should be irrelevant since morphologically dominated.

options = util.options()
options.nlive = 100                    # Num. live points for MultiNest
options.poly_order = 1                 # Allow the model to be divided by a polynomial.
                                       # of this order to flatten the continuum.
options.individual_polynomials = False
# Specify parts of the BC03 grid that we want to load (isochrones, IMF, metallicities)
options.grid_path = os.environ['GRID_PATH']+"/Padova1994/salpeter/bc2003_hr_stelib_m[metal]_salp_ssp.ised"
options.grid_param = ["metal"]
options.grid_paramlabels = [["52", "62", "72"]]
options.grid_paramvalues = [[0.008, 0.02, 0.05]]
options.grid_paramlogint = [True]
options.trim_wavelength = [3500., 10000.]
options.air2vac = False                 # Set to True if model grid should be converted to vac.
options.rescale_spectrum_errors = 1.   # Convenient way to rescale errors if necessary.
options.output_prefix = sys.argv[-1]   # Where the output goes.
options.age_lt_universe = True

# Input files
specdata = [bluSpec, redSpec]
specerror, specmaskfiles, photdata = None, None, None
options.tabulated_data = True

# Line spread function
# Here using a stacked profile of the target galaxies. There should be an odd number
# of lines, such that the central line is the center of the LSF, with lines at one
# pixels intervals.
lsfdata = np.loadtxt(bluLSF)
lsf1 = ArbitraryLSF(lsfdata[:,1])
lsfdata = np.loadtxt(redLSF)
lsf2 = ArbitraryLSF(lsfdata[:,1])
inst_lsf = [lsf1, lsf2]

if not os.path.exists(sys.argv[-1]): os.mkdir(sys.argv[-1])
specfit(specdata, specerror, photdata, priors,
        model, options, inst_lsf, specmaskfiles=specmaskfiles)
