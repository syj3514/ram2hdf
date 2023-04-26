# Simulation output
'''
> base:         the directory that contains `snaprepo`.
> snaprepo:     in `base`, `snaprepo` has raw RAMSES output data.
> halorepo:     in `base`, `halorepo` has the HaloMaker/GalaxyMaker data.
> ramsesver:    Ramses version (Ra, Ra3, Ra4). HAGN is Ra3. From Ra4, particles have `family` to identify category
> full_path:    For arbitrary path of any RAMSES
'''
base:str        = '/storage4/Horizon_AGN'
snaprepo:str    = 'snapshots'
halorepo:str    = 'galaxy'
ramsesver:str   = 'Ra3'
full_path:str   = None

# HAGN to hdf5 params
'''
> boxcut:       From galactic center, cutout range
> minmass:      Ignore galaxies below `minmass`
> maxmass:      Ignore galaxies above `maxmass`
> savedir:      the directory where HDF5 data will be saved.
> verbose:      log level
> nthread:      If fortran mode, omp_num_thread. If python mode, disabled
> calc_total:   calculate all particle numbers. (time-consuming)
> maxGB:        Memory limit. If exceed, it will flush loaded data
> functype:     For some functions, choose Fortran or Python mode
'''
boxcut:float    = 50 # kpc 
minmass:float   = 10**10 # Msol
maxmass:float   = 10**11.5 # Msol
savedir:str     = '/storage1/jeon/MAGPI/hdf5'
verbose:int     = 1
nthread:int     = 20
calc_total:bool = True
maxGB:float     = 10 #GB
functype:str    = 'fortran'#'python' or 'fortran'

# Unit conversion
Mpc2cm              = 3.08568e+24
kpc2cm              = 3.08568e+21
speed_of_light      = 299792.458
cm2kpc              = 3.24078e-22
cms2kms             = 1e-5
g2Msol              = 5.02785e-34
gcm32Msolkpc3       = 1.477e+31
g_constant_cgs      = 6.67430e-11
g_in_kpcMsolkms2    = 4.3009e-6
yr2sec              = 365*24*60*60