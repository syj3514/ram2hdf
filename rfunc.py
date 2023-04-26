# Seyoung Jeon
# syj3514@yonsei.ac.kr

import numpy as np
import os
import sys
import time
import h5py
import pickle
from IPython import get_ipython
from IPython.display import clear_output
import rparams
from rur import uri
import logging

class DotDict(dict):
    """dot.notation access to dictionary attributes"""
    def __getattr__(self, attr):
        return self.get(attr)
    __setattr__= dict.__setitem__
    __delattr__= dict.__delitem__

    def __getstate__(self):
        return self

    def __setstate__(self, state):
        self.update(state)
        self.__dict__ = self
p = {}
for key in rparams.__dict__.keys():
#    if not "_" in key:
    p[key] = rparams.__dict__[key]
p = DotDict(p)

if p.functype=='python':
    import read_ramses_python as readramses
    print("Read ramses using python")
elif p.functype=='fortran':
    from read_ramses_fortran import readramses
    print("Read ramses using fortran")
else:
    raise ImportError(f"`functype` in `rparams.py` should be `fortran` or `python`! (`{p.functype}` is given)")
chems = {
    'hagn':['H','O','Fe', 'C', 'N', 'Mg', 'Si'], 
    'yzics':[], 'nh':[], "fornax":[], "y2":[], "dm_only":[],
    "y3":['H', 'O', 'Fe', 'Mg', 'C', 'N', 'Si', 'S'], 
    "y4":['H', 'O', 'Fe', 'Mg', 'C', 'N', 'Si', 'S', 'D'], 
    "nc":['H', 'O', 'Fe', 'Mg', 'C', 'N', 'Si', 'S', 'D'],
    'nh2':['H', 'O', 'Fe', 'Mg', 'C', 'N', 'Si', 'S', 'D']}
chem = chems[p.mode]

###################################################################
###        Useful functions                                     ###
###################################################################
kernel=str(get_ipython())
if kernel is None:
    ipy = False
else:
    if ('zmqshell' in kernel) or ('jupyter' in kernel):
        ipy = True
    else:
        ipy = False

from logging.handlers import RotatingFileHandler
def custom_debugger(fname:str):
    logger_file_handler = RotatingFileHandler(fname, mode='a')

    logger_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(fmt=u'%(asctime)s [%(levelname)8s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    
    logger_file_handler.setFormatter(formatter)

    logging.captureWarnings(True)

    root_logger = logging.getLogger(fname)
    warnings_logger = logging.getLogger("py.warnings")
    root_logger.handlers = []
    warnings_logger.handlers = []
    root_logger.addHandler(logger_file_handler)
    warnings_logger.addHandler(logger_file_handler)
    root_logger.setLevel(logging.DEBUG)
    root_logger.info("Debug Start")
    root_logger.propagate=False
    return root_logger

def dprint(msg:str, debugger:logging.Logger, level='debug'):
    if debugger is not None:
        if level=='debug':
            debugger.debug(msg)
        elif level == 'info':
            debugger.info(msg)
        elif level == 'warning':
            debugger.warning(msg)
        elif level == 'critical':
            debugger.critical(msg)
        elif level == 'error':
            debugger.error(msg)
        else:
            debugger.info(msg)
    else:
        print(msg)

def cleanprint(num=1):
    "Use this function to delete the last line in the STDOUT"
    global ipy
    if ipy:
        clear_output()
    else:
        if not isinstance(num, int):
            num = int(num)
        for _ in range(num):
            # cursor up one line
            sys.stdout.write('\x1b[1A')
            # delete last line
            sys.stdout.write('\x1b[2K')

def printgal(gal, mode=p.mode):
    '''Summary of galaxy information'''
    iout = gal['timestep']
    made = 'GalaxyMaker'
    print(f"[{made}: {mode}] ID={gal['id']}, iout(istep)={iout}, logM={np.log10(gal['m']):.2f}")
    return f"[{made}: {mode}] ID={gal['id']}, iout(istep)={iout}, logM={np.log10(gal['m']):.2f}"

def pklload(fname):
    '''
    >>> array = pklload('path/fname.pickle')
    '''
    with open(fname, 'rb') as handle:
        try:
            arr = pickle.load(handle)
        except EOFError:
            arr = pickle.load(handle.read())
            # arr = {}
            # unpickler = pickle.Unpickler(handle)
            # # if file is not empty scores will be equal
            # # to the value unpickled
            # arr = unpickler.load()
    return arr

class timestamp():
    def __init__(self, msg, logger=None):
        self.ref = time.time()
        self.logger=logger
        dprint(msg, self.logger)

    def done(self, unit='sec', prefix="::"):
        units = {'ms':1000, 'sec':1, 'min':1/60}
        elapse = time.time() - self.ref
        dprint(f"{prefix} Done ({elapse*units[unit]:.2f} {unit})",self.logger)


###################################################################
###        Post processing                                      ###
###################################################################
# For values which are not saved in default

# minmass = 4061240 # Msol
def massinit(m, minmass=4061240):
    '''Calculate initial mass of star before mass loss'''
    return (m//minmass + 1) * minmass

def sfr_in_cell(snap, timewindow=0.1):
    global p
    c = snap.cell
    s = snap.part['star']
    readramses.sfrincell(
        c['x'], c['y'], c['z'], c['dx'],
        s['age', 'Gyr'], s['x'], s['y'], s['z'], s['m', 'Msol'],
        timewindow, p.nthread
        )
    return readramses.sfrs

nstar=None; ndm=None; ncell=None
def calc_totnum(snap):
    global p
    global nstar, ndm, ncell
    if ncell is None:
        if p.calc_total:
            cpulist = np.arange(snap.ncpu)+1
            readramses.read_part(f"{p.base}/{p.snaprepo}", snap.iout, cpulist, p.verbose, p.nthread)
            nstar = int(readramses.nstar)
            ndm = int(readramses.ndm)
            readramses.read_cell(f"{p.base}/{p.snaprepo}", snap.iout, cpulist, p.verbose)
            ncell = int(readramses.ncell_tot)
        else:
            nstar, ndm, ncell = 0, 0, 0
    return nstar, ndm, ncell

def set_attribute(group:h5py.Dataset, key:str, conversion, unit=None, readme=None):
    group[key].attrs["CGSConversionFactor"] = conversion
    group[key].attrs["aexp-scale-exponent"] = 0
    group[key].attrs["h-scale-exponent"] = 0
    if unit is None: unit="Code unit"
    group[key].attrs["unit"] = unit
    if readme is None: readme = f"{key} in {unit}"
    group[key].attrs["README"] = readme


###################################################################
###        Main function                                        ###
###################################################################
def save_hdf(snap:uri.RamsesSnapshot, target:np.ndarray, partids:np.ndarray, calc_total:bool=False, maxGB:float=50, logger:logging.Logger=None):
    global p, ipy
    ref = time.time()
    fname = f"{p.savedir}/{p.mode}/{snap.iout:03d}/{p.mode}_{snap.iout:03d}_MAGPI_{target['id']:06d}.hdf5"
    if p.verbose>0: dprint(printgal(target, mode=p.mode), logger)

    if os.path.isfile(fname):
        msg = f"`{fname}` is already saved!"
        # if p.verbose==1: cleanprint(1)
        if p.verbose>1:
            with h5py.File(fname, 'r') as f:
                dprint(f.keys(),logger)
            if p.verbose>0: dprint(msg, logger)
        return time.time()-ref, msg

    else:
        if p.verbose>0: timer0=timestamp(f" > Making `{fname}`...", logger=logger)

        # Set box & load part & load cell
        unitkpc = p.kpc2cm / snap.params['unit_l']
        unitMsol = 1 / p.g2Msol / snap.params['unit_m']
        snap.set_box_halo(target, p.boxcut*unitkpc, use_halo_radius=False)
        if p.verbose>1: timer = timestamp(" > Load part...",logger=logger)
        snap.get_part(nthread=p.nthread)
        star = snap.part['star']
        if p.verbose>1: timer.done(prefix="  ")
        if target['m'] > 1e11:
            if np.sum(star['m']/unitMsol) < target['m']/2:
                msg = f"Galaxy {target['id']} is too big ({np.log10(target['m']):.2f})!"
                if p.verbose>0: dprint(msg,logger)
                return time.time()-ref, msg
        if p.verbose>1: timer = timestamp(" > Load cell...",logger=logger)
        snap.get_cell()
        if p.verbose>1: timer.done(prefix="  ")


        if p.verbose>1: dprint(" > Convert to HDF5 format...",logger)
        
        with h5py.File(fname, 'w') as f:
            compression = 'gzip'

            ######################
            #       HEADER       #
            ######################
            header = f.create_group("Header")
            if p.verbose>0: dprint("     > Writing Header...",logger)
            
            # From Snapshot data
            header.attrs["BoxSize"] = snap.params['boxsize']
            header.attrs["BoxSize_physical"] = snap.params['boxsize_physical']
            header.attrs["BoxSize_comoving"] = snap.params['boxsize_comoving']
            header.attrs["Redshift"] = snap.params['z']
            header.attrs["Lookback_Time"] = snap.params['lookback_time']
            header.attrs["HubbleParam"] = snap.params['h']
            header.attrs["MassTable"] = np.array([0., 0., 0., 0., 0., 0.])
            nstar, ndm, ncell = calc_totnum(snap)
            header.attrs["NumPart_Total"] = np.array([ncell, ndm, 0, 0, nstar, 0])
            header.attrs["NumPart_This"] = np.array([len(snap.cell.table), 0, 0, 0, len(star.table), 0])
            header.attrs["RunLabel"] = p.simname #dtype=h5py.string_dtype(encoding='utf-8')
            header.attrs["ExpansionFactor"] = snap.params['aexp']
            header.attrs["Omega0"] = snap.params['omega_m']
            header.attrs["OmegaBaryon"] = snap.params['omega_b']
            header.attrs["OmegaLambda"] = snap.params['omega_l']
            header.attrs["unit_l"] = snap.params['unit_l']
            header.attrs["unit_m"] = snap.params['unit_m']
            header.attrs["unit_d"] = snap.params['unit_d']
            header.attrs["unit_t"] = snap.params['unit_t']
            # From GalaxyMaker
            header.attrs["id"] = target['id']
            header.attrs["m_msol"] = target['m']
            header.attrs["x_kpc"] = target["x"]/unitkpc
            header.attrs["y_kpc"] = target["y"]/unitkpc
            header.attrs["z_kpc"] = target["z"]/unitkpc
            header.attrs["vx_kms"] = target["vx"]
            header.attrs["vy_kms"] = target["vy"]
            header.attrs["vz_kms"] = target["vz"]
            header.attrs["hostid"] = target['host']
            header.attrs["nbsub"] = target['nbsub']
            header.attrs["_README"] = "\
                # BoxSize: The linear extent of the simulation cube in Mpc/h.\n\
                # Redshift: The redshift of the simulation at this snapshot.\n\
                # HubbleParam: The Hubble Parameter (H0 / 100 km/s/Mpc).\n\
                # MassTable: The mass of particles for a given particle type. If 0, see the Mass table within the respective PartType group.\n\
                # NumPart_Total: The number of particles of each particle type within the simulation.\n\
                # NumPart_This: The number of particles of each particle type within this file.\n\
                # RunLabel: Name of the simulation.\n# ExpansionFactor: The expansion scale factor.\n\
                # Omega0: The matter density parameter.\n# OmegaBaryon: The baryon density parameter.\n\
                # OmegaLambda: The dark energy density parameter.\n\
                # unit_l: length conversion of codeunit into cm (cgs) unit\n\
                # unit_m: mass conversion of codeunit into gram (cgs) unit\n\
                # unit_d: density conversion of codeunit into g/cc (cgs) unit\n\
                # unit_t: time conversion of codeunit into sec (cgs) unit\n\
                # hostid: If `hostid`!=`id`, subhalo(`id`) of hosthalo(`hostid`)\n\
                # nbsub: The number of subhalos"


            ######################
            #     PARTTYPE0      #
            ######################
            # phy_param = simul_param * a**aexp-scale-exponent * h**h-scale-exponent * CGSConversionFactor
            part0 = f.create_group("PartType0")
            if p.verbose>1: dprint("     > Writing PartType0...",logger)
            
            part0.create_dataset(f"Coordinates", data=snap.cell['pos'], compression=compression) # ★★★
            set_attribute(part0, "Coordinates", snap.params['unit_l'], readme="The (x,y,z) coordinates of each cell center")
            part0.create_dataset(f"Density", data=snap.cell['rho'], compression=compression) # ★★★
            set_attribute(part0, "Density", snap.params['unit_d'], readme="The gas cell density")
            part0.create_dataset(f"dx", data=snap.cell['dx'], compression=compression)
            set_attribute(part0, "dx", snap.params['unit_l'], readme="The cell size of a side")
            part0.create_dataset(f"Mass", data=snap.cell['m'], compression=compression) # ★★★
            set_attribute(part0, "Mass", snap.params['unit_m'], readme="Gas cell mass")
            part0.create_dataset(f"ParticleIDs", data=np.zeros(len(snap.cell['x'])), compression=compression) # ★★★
            set_attribute(part0, "ParticleIDs", 0, unit='None', readme="AMR cells don't have unique IDs because they are changed at every timestep")

            chem = part0.create_group("ElementAbundance")
            dat = snap.cell['C'] if 'C' in chem else np.zeros(len(snap.cell['x']))-1
            chem.create_dataset(f"Carbon", data=dat, compression=compression)
            set_attribute(chem, "Carbon", 0, unit="Mass fraction", readme="Mass fraction of Carbon for a given cell mass")
            dat = snap.cell['O'] if 'O' in chem else np.zeros(len(snap.cell['x']))-1
            chem.create_dataset(f"Oxygen", data=dat, compression=compression)
            set_attribute(chem, "Oxygen", 0, unit="Mass fraction", readme="Mass fraction of Oxygen for a given cell mass")
            dat = snap.cell['H'] if 'H' in chem else np.zeros(len(snap.cell['x']))-1
            chem.create_dataset(f"Hydrogen", data=dat, compression=compression)
            set_attribute(chem, "Hydrogen", 0, unit="Mass fraction", readme="Mass fraction of Hydrogen for a given cell mass")

            part0.create_dataset(f"Metallicity", data=snap.cell['metal'], compression=compression) # ★★★
            set_attribute(part0, "Metallicity", 0, unit="Mass fraction", readme="Mass fraction of all elements heavier than Helium")

            if p.verbose>1: timer = timestamp("     > Calculate SFR...",logger=logger)
            SFRs = sfr_in_cell(snap, timewindow=0.1)
            if p.verbose>1: timer.done(prefix="      ")
            part0.create_dataset(f"StarFormationRate", data=SFRs, compression=compression) # ★★★
            set_attribute(part0, "StarFormationRate", 1 / p.g2Msol / p.yr2sec, unit="Msol/yr", readme="(Last 100Myr) star formation rate")
            part0.create_dataset(f"Temperature", data=snap.cell['T', 'K'], compression=compression)
            set_attribute(part0, "Temperature", 1, unit="Kelvin", readme="Temperature: Pressure*density. k_B = 1.38064852E-16. m_H = 1.6737236E-24")
            part0.create_dataset(f"Velocity", data=snap.cell['vel'], compression=compression) # ★★★
            set_attribute(part0, "Velocity", snap.params['unit_l'] / snap.params['unit_t'], unit="code unit", readme="The flow velocity (vx,vy,vz) coordinates of each cell")
            part0.create_dataset(f"Pressure", data=snap.cell['P', 'Ba'], compression=compression)
            set_attribute(part0, "Pressure", 1, unit="Ba", readme="The gas pressure of each cell")
                
                
            ######################
            #     PARTTYPE4      #
            ######################
            # phy_param = simul_param * a**aexp-scale-exponent * h**h-scale-exponent * CGSConversionFactor
            part4 = f.create_group("PartType4")
            if p.verbose>1: dprint("     > Writing PartType4...",logger)
            
            part4.create_dataset(f"Coordinates", data=star['pos'], compression=compression)
            set_attribute(part4, "Coordinates", snap.params['unit_l'], unit="code unit", readme="The (x,y,z) coordinates of each particle.")
            initialmass = massinit(star['m','Msol']) /p.g2Msol /snap.params['unit_m'] # in code unit
            part4.create_dataset(f"InitialMass", data=initialmass, compression=compression)
            set_attribute(part4, "InitialMass", snap.params['unit_m'], unit="code unit", readme="The particle mass at formation time (but not exact)")
            part4.create_dataset(f"Mass", data=star['m'], compression=compression)
            set_attribute(part4, "Mass", snap.params['unit_m'], unit="code unit", readme="Particle mass")
            part4.create_dataset(f"ParticleIDs", data=np.abs(star['id']), compression=compression)
            set_attribute(part4, "ParticleIDs", 1, unit="None", readme="Unique particle identification number")
            part4.create_dataset(f"Metallicity", data=star['metal'], compression=compression)
            set_attribute(part4, "Metallicity", 1, unit="Mass fraction", readme="Mass fraction of all elements heavier than Helium")
            vals = np.interp(star['epoch'], snap.cosmo_table['u'], snap.cosmo_table['aexp'])
            part4.create_dataset(f"StellarFormationTime", data=vals, compression=compression)
            set_attribute(part4, "StellarFormationTime", 1, unit="None", readme="Expansion factor, a, when the star was born")
            part4.create_dataset(f"Velocity", data=star['vel'], compression=compression)
            set_attribute(part4, "Velocity", snap.params['unit_l'] / snap.params['unit_t'], unit="code unit", readme="The velocity (vx,vy,vz) coordinates of particle")
            vals = np.isin(np.abs(star['id']), partids, assume_unique=True)
            part4.create_dataset(f"Member", data=vals, compression=compression)
            set_attribute(part4, "Member", 1, unit="None", readme="Membership boolean by GalaxyMaker")

        size = (snap.part_data.nbytes+snap.cell_data.nbytes) / 2**30 #GB
        msg = f" > `Save (cache {size:.2f} GB)"
        if p.verbose>0: timer0.done(prefix=msg)
            # if p.verbose==1:
            #     if not ipy:
            #         cleanprint(10)
            #     else:
            #         clear_output()
            
        if size > maxGB:
            snap.clear()
        return time.time()-ref, msg
