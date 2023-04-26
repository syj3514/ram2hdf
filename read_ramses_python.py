from rur.fortranfile import FortranFile
import numpy as np
import glob
import gc, os

# def read_cell(iout, base='.', path_in_repo='snapshots', verbose=True):
ncell_tot = 0
def read_cell(repo, iout, cpu_list, verbose):
    global ncell_tot
    ncell_tot = 0
    fname = f"{repo}/output_{iout:05d}/amr_{iout:05d}.out00001"
    # if verbose: print(f"\tfrom `{fname}` ...")
    with FortranFile(fname, mode='r') as f:
        ncpu,                   = f.read_ints()
        ndim,                   = f.read_ints(); assert ndim==3
        f.skip_records(1)
        nlevelmax,              = f.read_ints()
        f.skip_records(1)
        nboundary,              = f.read_ints()
        ngridfile = np.empty((ncpu+nboundary, nlevelmax), dtype='i4')
    twotondim=2**ndim
    skip_amr = 3 * (2**ndim + ndim) + 1

    # if verbose: print(f"2. Check total number of grids")
    # ncell_tot = 0
    # iterobj = tqdm(range(1, ncpu+1), "\tLoad amr files") if verbose else range(1, ncpu+1)
    # for icpu in iterobj:
    for icpu in cpu_list:
        fname = f"{repo}/output_{iout:05d}/amr_{iout:05d}.out{icpu:05d}"
        with FortranFile(fname, mode='r') as f:
            f.skip_records(21)
            numbl                   = f.read_ints()
            for ilevel in range(nlevelmax):
                ngridfile[:,ilevel]=numbl[ncpu*ilevel : ncpu*(ilevel+1)]
            f.skip_records(7)
            if nboundary>0: f.skip_records(3)
            levels, cpus = np.where(ngridfile.T>0)
            for ilevel, jcpu in zip(levels, cpus+1):
                f.skip_records(3)
                if jcpu==icpu:
                    f.skip_records(3*ndim+1)
                    for _ in range(twotondim):
                        son = f.read_ints()
                        if 0 in son:
                            ncell_tot += len(son.flatten())-np.count_nonzero(son)
                    f.skip_records(2*twotondim)
                else:
                    f.skip_records(skip_amr)

    # return ncell_tot

nstar = 0
ndm = 0
def read_part(repo, iout, cpu_list, verbose, nthread):
    global nstar, ndm
    nstar = 0; ndm = 0
    # files = glob.glob(f"{repo}/output_{iout:05d}/part*out*")
    files = [f"{repo}/output_{iout:05d}/part_{iout:05d}.out{icpu:05d}" for icpu in cpu_list]
    header = f"{repo}/output_{iout:05d}/header_{iout:05d}.txt"
    sinkinfo = f"{repo}/output_{iout:05d}/sink_{iout:05d}.info"

    isfamily = False
    # (NH, NH2, Fornax, NC)
    if os.path.isfile(header): 
        with open(header, "rt") as f:
            temp = f.readline()
            # (Fornax, NH2, NC)
            if "Family" in temp: 
                isfamily = True
                ntracer_tot = int( f.readline()[14:] ) # other_tracer
                for _ in range(5): # tracers of debris, cloud, star, other, gas
                    ntracer_tot += int( f.readline()[14:] )
                ndm = int( f.readline()[14:] )
                nstar = int( f.readline()[14:] )
                ncloud_tot = int( f.readline()[14:] )
                npart_tot = ntracer_tot + ndm + nstar + ncloud_tot
                for _ in range(3): # debris, other, undefined
                    npart_tot += int( f.readline()[14:] )
        
                if os.path.isfile(sinkinfo): # (NH2, NC)
                    with open(sinkinfo, 'rt') as f:
                        nsink_tot = int(f.readline().split()[-1])
                else: # (Fornax)
                    with FortranFile(f"{files[0]}", mode='r') as f:
                        f.skip_records(7)
                        nsink_tot = f.read_ints(np.int32)[0]
            # (NH)
            else: 
                npart_tot = int(f.readline()); f.readline()
                ndm = int(f.readline()); f.readline()
                nstar = int(f.readline()); f.readline()
                nsink_tot = int(f.readline()); f.readline()
                ncloud_tot = nsink_tot * 2109
                ntracer_tot = 0
    # (hagn, yzics)
    else: 
        with FortranFile(f"{files[0]}", mode='r') as f:
            f.skip_records(4) # ncpu, ndim, npart, localseed(+tracer_seed)
            nstar = f.read_ints(np.int32)[0]
            f.skip_records(2) # mstar_tot, mstar_lost
            nsink_tot = f.read_ints(np.int32)[0]
            ncloud_tot = 2109 * nsink_tot
        ntracer_tot = 0
        npart_tot = 0
        for file in files:
            with FortranFile(f"{file}", mode='r') as f:
                f.skip_records(2)
                npart_tot += f.read_ints(np.int32)[0]
        ndm = npart_tot - nstar - ncloud_tot

    # with FortranFile(f"{files[0]}", mode='r') as f:
    #     f.skip_records(4) # ncpu, ndim, npart, localseed(+tracer_seed)
    #     nstar = f.read_ints(np.int32)[0]
    #     f.skip_records(2) # mstar_tot, mstar_lost
    #     nsink = f.read_ints(np.int32)[0]
    #     ncloud = 2109 * nsink
    # npart_tot = 0
    # for file in files:
    #     with FortranFile(f"{file}", mode='r') as f:
    #         f.skip_records(2)
    #         npart_tot += f.read_ints(np.int32)[0]
    # ndm = npart_tot - nstar - ncloud


sfrs = np.array([])
def sfrincell(cx, cy, cz, cdx, sage, sx, sy, sz, sm, timewindow, nthread):
    '''Calculate SFR of each cell'''
    global sfrs
    cell_xmin = cx - cdx/2
    cell_xmax = cx + cdx/2
    cell_ymin = cy - cdx/2
    cell_ymax = cy + cdx/2
    cell_zmin = cz - cdx/2
    cell_zmax = cz + cdx/2
    
    # Young star : age < `timewindow` Gyr
    yarg = sage < timewindow
    if True in yarg:
        ystar_x = sx[yarg]
        ystar_y = sy[yarg]
        ystar_z = sz[yarg]
        ystar_m = sm[yarg]/1e8
        
        condx = (cell_xmin[:, None] <= ystar_x) & (cell_xmax[:, None] > ystar_x)
        condy = (cell_ymin[:, None] <= ystar_y) & (cell_ymax[:, None] > ystar_y)
        condz = (cell_zmin[:, None] <= ystar_z) & (cell_zmax[:, None] > ystar_z)
        cond = condx & condy & condz # Boolean array; find a cell where each young star is
        where = np.where(cond == True)
        cell_indices = where[0]
        ystar_indices = where[1]
        
        SFRs = np.zeros(len(cell_xmin))
        np.add.at(SFRs, cell_indices, ystar_m[ystar_indices])
    else:
        print("No Young stars")
        SFRs = np.zeros(len(cell_xmin))
    sfrs = SFRs
    # return SFRs


def close():
    gc.collect()