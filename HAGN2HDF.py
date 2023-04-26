from hfunc import *
from rur import uri, uhmi
uri.timer.verbose = 0
from tqdm import tqdm
import time
import os
import traceback

# Read Input
args = sys.argv
if len(args)<2:
    print(f"\nYou should specify iout!:")
    print(f" $ python3 run_hagn2hdf.py <iout>\n")
    raise AssertionError("STOP")

# Make save directory and log file
iout = int(args[1])
if not os.path.isdir(f"{p.savedir}/{iout:03d}"):
    os.mkdir(f"{p.savedir}/{iout:03d}")
logname = f"{p.savedir}/hagn2hdf_{iout:03d}.log"
logger = custom_debugger(logname)

# Load galaxy
snap = uri.RamsesSnapshot(p.base, iout, mode='hagn', box=None, path_in_repo=p.snaprepo, longint=False)
gals, gpids = uhmi.HaloMaker.load(snap, path_in_repo=p.halorepo, galaxy=True, double_precision=False, full_path=p.full_path, load_parts=True)
mask = (gals['m'] >= p.minmass) & (gals['m'] <= p.maxmass)
gals, gpids = uhmi.HaloMaker.cut_table(gals, gpids, mask)
print(f"Ngal = {len(gals)}, ({np.min(gals['id'])} ~ {np.max(gals['id'])})\nNparts = {len(gpids)}\n")
dprint(f"Ngal = {len(gals)}, ({np.min(gals['id'])} ~ {np.max(gals['id'])})\nNparts = {len(gpids)}\n", logger)
table = np.vstack((gals['id'], gals['m'], gals['x']*snap.unit_l/p.kpc2cm, gals['y']*snap.unit_l/p.kpc2cm, gals['z']*snap.unit_l/p.kpc2cm)).T
np.savetxt(f"{p.savedir}/catalogue_{iout:03d}.csv", table, fmt='%d,%e,%e,%e,%e', header='id,m_Msol,x_kpc,y_kpc,z_kpc', comments='#')

nparts = gals['nparts']
cparts = np.insert( np.cumsum(nparts), 0, 0 )

# Iteration
for i, gal in enumerate(gals):
    try:
        partids = gpids[ cparts[ i : i+1 ] ]
        elapse, cache = save_hdf(snap, gal, partids, calc_total=True, maxGB=p.maxGB, logger=logger)
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        logger.error(traceback.format_exc())
        logger.error(e)
        print("Iteration is terminated")
        os._exit(1)
