# ram2hdf
by Seyoung Jeon (syj3514@yonsei.ac.kr)  
The `ram2hdf` reads RAMSES raw output and converts to HDF5 which is compatible with `SimSpin`.  

---

## Preparation
Install `rur` from https://github.com/sanhancluster/rur  
Also, you should have proper RAMSES output and `GalaxyMaker` data  
(Currently support: Horizon-AGN, YZiCS, NewHorizon, NewHorizon2, FORNAX, NewCluster)  

---

## Install
Download the code from github repository  
```bash
~$ git clone https://github.com/syj3514/ram2hdf.git
```

Add path:
```bash
~$ conda develop ~/ram2hdf
```

(Optional) If you want to use fortran subroutine as python function, run `f2py.sh`:
```bash
~/ram2hdf$ ./f2py.sh
```
:warning: **Please confirm your fortran compiler before f2py**

---

## Parameter setting
Set parameters in `rparams.py`  
####Example:
```python
# Simulation output
'''
> base:         the directory that contains `snaprepo`.
> snaprepo:     in `base`, `snaprepo` has raw RAMSES output data.
> halorepo:     in `base`, `halorepo` has the HaloMaker/GalaxyMaker data.
> ramsesver:    Ramses version (Ra, Ra3, Ra4). HAGN is Ra3. From Ra4, particles have `family` to identify category
> full_path:    For arbitrary path of any RAMSES output
'''
mode:str        = 'hagn' # See rur/config.py (hagn, yzics, nh, nh2, nc...)
simname:str     = 'HorizonAGN'
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
savedir:str     = '/storage1/jeon/MAGPI/testram'
verbose:int     = 0
nthread:int     = 20
calc_total:bool = True
maxGB:float     = 50 #GB
functype:str    = 'python'#'python' or 'fortran'

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
```
### Simulation output
Desired RAMSES output data file structure is like:
```plain
base
|---snapshots
|   |---output_00001
|   |   |   info_00001.txt
|   |   |   amr_00001.out00001
|   |   |   amr_00001.out00002
|   |   |   ...
|   |   |   part_00001.out04096
|   |   
|   |---output_00002
|   |   ...
|   |---output_00782
|
|---galaxy
|   |   tree_bricks00001
|   |   tree_bricks00002
|   |   ...
|   |   tree_bricks00782
|
|---halo
|   |   tree_bricks00001
|   |   tree_bricks00002
|   |   ...
|   |   tree_bricks00782

```
Following the structure, set `base`, `snaprepo`, `halorepo`, `full_path` properly.  
`ramsesver` indiciates the ramses code version.  
`Horizon-AGN` and `NewHorizon` is `Ra3`, and `NewHorizon2` and `NewCluster` is `Ra4`.  
`Ra3` classifies particles in `part` file using the sign of ID and epoch of particles.  
`Ra4` does it using `family`.  

### RAMSES to HDF5 params
See https://kateharborne.github.io/SimSpin/

### Running params
The result hdf5 files are saved in `savedir`.  
`verbose` for debugging level.  
`nthread` is used only for fortran mode as `omp_num_threads`.  
If memory usage exceeds `maxGB`, cached data will be deallocated.  

---

## Run 
```bash
~/ram2hdf$ python3 RAM2HDF.py <iout>
```
You can see the log, `ram2hdf_{mode}_{iout:03d}.log`, in `savedir`.  
Output file will saved as `{savedir}/{iout:03d}/{mode}_{iout:03d}_MAGPI_{Galaxy_ID:06d}.hdf5`

---

## HDF5 format
See https://kateharborne.github.io/SimSpin/examples/generating_hdf5.html  
or See `save_hdf` function in `rfunc.py`
