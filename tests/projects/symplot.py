import json
import os
import numpy as np
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

def angle(v1, v2):
    """Return unsigned angle between two 3d vectors

    Parameters
    ----------
    v1 : 1x3 array
    v1 : 1x3 array

    Returns
    -------
    float

    """
    angle = np.arccos(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle * 180 / np.pi


def transform_energy_to_surface_energy(E,avec,bvec):
    area=np.linalg.norm(np.cross(avec,bvec))
    ev_to_joules=1.602176565e-19
    angstrom_to_meter=1e-10
    return ev_to_joules*E/(area*angstrom_to_meter*angstrom_to_meter)


def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

    return ax


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    return set_axes_radius(ax, origin, radius)


def load_slab_plane_vectors(slabfile):
    with open(slabfile) as f:
        content = f.readlines()

    aline = content[2]
    bline = content[3]
    cline = content[3]

    a = np.array([float(x) for x in aline.split()])
    b = np.array([float(x) for x in bline.split()])
    c = np.array([float(x) for x in cline.split()])

    return a, b


def load_record(record_file):
    with open(record_file) as f:
        records = json.load(f)
    return records


def write_json(data,target_file):
    with open(target_file,'w') as f:
        json.dump(data,f)
    return


def scatter_equivalent_shifts(ax,data_slice,avec,bvec):
    awall=data_slice.loc[data_slice["a_index"]==0].copy()
    bwall=data_slice.loc[data_slice["b_index"]==0].copy()

    awall["a_index"]=len(bwall)
    bwall["b_index"]=len(awall)

    corner=data_slice.loc[(data_slice["b_index"]==0) & (data_slice["a_index"]==0)].copy()
    corner["a_index"]=len(bwall)
    corner["b_index"]=len(awall)

    shift_records=pd.concat([data_slice,awall,bwall,corner])

    aix = shift_records["a_index"]
    bix = shift_records["b_index"]

    acount = max(set(aix))
    bcount = max(set(bix))

    eqs = shift_records["equivalent_shifts"].astype(str)

    unq = list(set(eqs))
    unq.sort()
    colors = np.arange(len(unq)) / len(unq)
    # np.random.shuffle(colors)
    cmap = {u: c for u, c in zip(list(unq), colors)}

    pts = np.array(
        [avec * a / acount + bvec * b / bcount for a, b in zip(aix, bix)])

    c = [cmap[eq] for eq in eqs]

    ax.scatter(pts[:, 0], pts[:, 1], c=c,s=85, cmap=matplotlib.cm.get_cmap("twilight"))

    ax.set_aspect('equal')

    ax.set_xlabel("x $\mathrm{[\AA]}$")
    ax.set_ylabel("y $\mathrm{[\AA]}$")

    return ax


def load_oszicar_energy(oszicar_path):
    with open(oszicar_path) as f:
        content = f.readlines()

    return float(content[-1].split()[2])


def unwind_data(records):
    paths = [p for p in records["shift-cleave"]]
    avals = [
        records["shift-cleave"][p]["a_index"] for p in records["shift-cleave"]
    ]
    bvals = [
        records["shift-cleave"][p]["b_index"] for p in records["shift-cleave"]
    ]
    angles = [
        records["shift-cleave"][p]["angle"] for p in records["shift-cleave"]
    ]
    ids = [records["shift-cleave"][p]["id"] for p in records["shift-cleave"]]
    cleaves = [
        records["shift-cleave"][p]["cleavage"] for p in records["shift-cleave"]
    ]
    equivs = [
        records["shift-cleave"][p]["equivalent_structures"]
        for p in records["shift-cleave"]
    ]

    sequivs=[[[int(x) for x in eq.split(':')[0:2]] for eq in group] for group in equivs]

    unwinded = pd.DataFrame.from_dict({
        "path": paths,
        "a_index": avals,
        "b_index": bvals,
        "angle": angles,
        "id": ids,
        "cleavage": cleaves,
        # "equivalent_structures": equivs
        "equivalent_shifts": sequivs
    })

    unwinded = unwinded.sort_values(["a_index", "b_index"])
    return unwinded


def grid_divisions(unwinded):
    return max(unwinded["a_index"]) + 1, max(unwinded["b_index"]) + 1


def fractional_coordinates(unwinded):
    amax, bmax = grid_divisions(unwinded)
    A = np.array(unwinded["a_index"]) / amax
    B = np.array(unwinded["b_index"]) / bmax

    return A,B


def bin_data_slice(unwinded, avec, bvec, column):
    L = np.array([avec[0:2], bvec[0:2]]).T


    amax, bmax = grid_divisions(unwinded)
    A, B = fractional_coordinates(unwinded)
    X, Y = L.dot(np.array([A, B]))

    X = X.reshape(amax, bmax)
    Y = Y.reshape(amax, bmax)
    Z = np.array(unwinded[column]).reshape(X.shape)

    return X, Y, Z


def data_slice_orbits(unwinded):
    id_orbits = []
    for s in unwinded["equivalent_shifts"]:
        if s not in id_orbits:
            id_orbits.append(s)
    return id_orbits


def load_energies(unwinded,avec,bvec,project):
    unwinded["energy"] = np.nan

    # oszis = glob.glob("./Al-FCC.chain/shift__*/cleave__*/OSZICAR")
    # oszis = glob.glob("./Mg-prismatic.chain/shift__*/cleave__*/OSZICAR")
    # oszis = glob.glob("./Mg-pyramidal1.chain/shift__*/cleave__*/OSZICAR")
    oszis = glob.glob(os.path.join(project+".chain","shift__*","cleave__*","OSZICAR"))

    energies = {f[0:-8]: load_oszicar_energy(f) for f in oszis}
    for p in energies:
        unwinded.loc[unwinded["path"] == p, "raw_energy"] = energies[p]
    unwinded["energy"]=transform_energy_to_surface_energy(unwinded["raw_energy"],avec,bvec)

    return unwinded


def plot_gamma_surface(unwinded_slice,slab_file):
    avec, bvec = load_slab_plane_vectors(slab_file)

    X, Y, Z = bin_data_slice(unwinded_slice, avec, bvec)
    X1, Y1, Z1 = X, Y, Z

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z)
    ax = set_axes_equal(ax)
    plt.show()


def unfold_orbits(unwinded_slice,column):
    orbits = data_slice_orbits(unwinded_slice)

    for o in orbits:
        equivs=unwinded_slice.loc[unwinded_slice["equivalent_shifts"].astype(str)
                           == str(o),column]
        not_nan=equivs[equivs.notnull()]
        assert(len(not_nan)>0)

        e = float(not_nan)
        unwinded_slice.loc[unwinded_slice["equivalent_shifts"].astype(str)
                           == str(o),column] = e
    return unwinded_slice



def equilibrium_energy(unwinded):
    return float(
        unwinded.loc[(unwinded["a_index"] == 0) & (unwinded["b_index"] == 0) &
                     (unwinded["cleavage"] == 0.0) &
                     (unwinded["angle"] == 0.0)]["energy"])


def surface_energy(unwinded):
    eqe = equilibrium_energy(unwinded)
    maxcleave = max(unwinded["cleavage"])

    surface_energies = (unwinded.loc[unwinded["cleavage"] ==
                                    maxcleave]["energy"] - eqe)*0.5
    if max(surface_energies)-min(surface_energies)>0.001:
        raise ValueError("Surface energies don't match at different shifts. Did you separate enough?")

    return np.mean(surface_energies)


def shift_energies_to_surface_energy(unwinded):
    eqe = equilibrium_energy(unwinded)
    gamma=surface_energy(unwinded)
    unwinded["energy"]-=(2*gamma+eqe)
    return unwinded


def uber(d, d0, gamma, lamb, shift):
    delta=d-d0
    return -2 * gamma * (1 + delta / lamb + 0) * np.exp(-delta / lamb) + shift


def uber_fit(unwinded_slice):
    # -25.359233540234765 -23.931949647554593
    # unwinded_slice=unwinded_slice.loc[unwinded_slice["cleavage"]>-0.01]
    # unwinded_slice=unwinded_slice.loc[unwinded_slice["energy"]<-24.0]
    unwinded_slice = unwinded_slice.sort_values(["cleavage"])
    sigma=np.ones(len(unwinded_slice))
    # sigma[0::10]=0.001

    popt, pcov = curve_fit(uber, unwinded_slice["cleavage"], unwinded_slice["energy"],sigma=sigma,absolute_sigma=True)

    return popt,pcov


def plot_shifted_uber_fit(ax,unwinded_slice,popt):
    deltas = np.linspace(
        min(unwinded_slice["cleavage"]), max(unwinded_slice["cleavage"]), 1000)
    fitenergy = uber(deltas, *popt)
    ax.set_xlabel("$d$ $\mathrm{[\AA]}$")
    ax.set_ylabel("Energy $\mathrm{[J/m^2]}$")
    # ax.scatter(unwinded_slice["cleavage"]-popt[0], unwinded_slice["energy"])
    # ax.plot(deltas-popt[0], fitenergy, "g--")
    ax.plot(deltas, fitenergy-popt[3], "g--", color='k',zorder=1)
    ax.scatter(unwinded_slice["cleavage"], unwinded_slice["energy"]-popt[3],zorder=2)

    # ax.set_ylim([-1.8,0.1]) #Al
    # ax.set_ylim([-1.4,0.1]) #Mg prismatic
    ax.set_ylim([-1.5,0.1]) #Mg pyramidal2
    # -25.359233540234765 -23.931949647554593

    vertical_dash=False
    if vertical_dash:
        ax.axvline(x=popt[0],ls="-.",c='k')
        ax.annotate("$d_0$",(popt[0]+0.2,0.0-0.05))

    return ax


def plot_surface_heatmap(ax,X,Y,Z,avec,bvec,xlabel,ylabel,cbar_label):
    pc=ax.pcolormesh(X, Y, Z)
    ax.set_aspect('equal')

    im_ratio=(max(Y.ravel())-min(Y.ravel()))/(max(X.ravel())-min(X.ravel()))
    cbar=plt.colorbar(pc,ax=ax,fraction=0.047*im_ratio,pad=0.04)
    cbar.set_label(cbar_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax


def plot_gamma_surface_heatmap(ax,unwinded_slice,avec,bvec):
    X, Y, Z = bin_data_slice(unwinded_slice, avec, bvec, "energy")
    return plot_surface_heatmap(ax,X,Y,Z-max(Z.ravel()),avec,bvec,"x $\mathrm{[\AA]}$","y $\mathrm{[\AA]}$","Energy [eV]")


def main():

    project="Al-FCC"
    # project="Mg-pyramidal1"
    # project="Mg-prismatic"
    # project="Mg-pyramidal2"

    records = load_record("./Al-FCC.chain/record.json")
    # records = load_record("./Mg-prismatic.chain/record.json")
    # records = load_record("./Mg-pyramidal1.chain/record.json")
    # records = load_record(os.path.join(project+".chain","record.json"))
    # records = load_record("./Mg-pyramidal2.chain/record.json")

    # print(json.dumps(records["shift-cleave"],indent=4,sort_keys=True))

    avec, bvec = load_slab_plane_vectors("./Al-FCC.chain/slab.vasp")
    # avec, bvec = load_slab_plane_vectors("./Mg-prismatic.slices/aligned_sliced_prim.vasp")
    # avec, bvec = load_slab_plane_vectors(os.path.join(project+".slices","aligned_sliced_prim.vasp"))
    # avec, bvec = load_slab_plane_vectors("./Mg-pyramidal2.slices/aligned_sliced_prim.vasp")

    unwinded = unwind_data(records)
    nocleave = unwinded.loc[unwinded["cleavage"] == 0.0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    scatter_equivalent_shifts(ax,nocleave,avec,bvec)
    plt.savefig("./figs/shift_symmetry.pdf",bbox_inches='tight',pad_inches=0)

    unwinded = load_energies(unwinded,avec,bvec,project)
    unfolded=[unfold_orbits(unwinded.loc[unwinded["cleavage"]==cleave].copy(),"energy") for cleave in set(unwinded["cleavage"])]
    unwinded=pd.concat(unfolded)

    # unwinded = shift_energies_to_surface_energy(unwinded)
    # gamma=surface_energy(unwinded)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax=plot_gamma_surface_heatmap(ax,nocleave,avec,bvec)
    # plt.show()


    noshift = unwinded.loc[(unwinded["a_index"] == 0) &
                           (unwinded["b_index"] == 0)]

    id_orbits=data_slice_orbits(nocleave)
    prototypes=[o[0] for o in id_orbits]

    uberfits=nocleave.copy()


    # for a,b in prototypes:
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     shiftdata=unwinded.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b)]
    #     ax.scatter(shiftdata["cleavage"],shiftdata["energy"])
    #     ax.set_ylim(min(noshift["energy"]),max(noshift["energy"]))
    #     plt.savefig("./figs/{}.raw{}.{}.pdf".format(project,a,b),bbox_inches='tight',pad_inches=0)

    #stack 0,0; 4,4; 8,8 UBER plots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for a,b in prototypes:
        if a != b:
            continue
        if a not in [0,4,8]:
            continue

        shiftdata=unwinded.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b)]
        popt, pcov = uber_fit(shiftdata)

        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_d0"]=popt[0]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_gamma"]=popt[1]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_lambda"]=popt[2]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_shift"]=popt[3]

        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_tmp"]=uber(popt[0],*popt)

        ax=plot_shifted_uber_fit(ax,shiftdata,popt)
        # ax.set_ylim(min(noshift["energy"]),max(noshift["energy"]))
    ax.set_ylim((-1.75,0.1))
    plt.savefig("./figs/{}.uberstack.pdf".format(project),bbox_inches='tight',pad_inches=0)
    plt.show()
    exit()

    for a,b in prototypes:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        shiftdata=unwinded.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b)]
        popt, pcov = uber_fit(shiftdata)

        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_d0"]=popt[0]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_gamma"]=popt[1]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_lambda"]=popt[2]
        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_shift"]=popt[3]

        uberfits.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b),"_tmp"]=uber(popt[0],*popt)

        print(uber(popt[0],*popt)-popt[3],-2*popt[1])
        print(a,b)

        ax=plot_shifted_uber_fit(ax,shiftdata,popt)
        # ax.set_ylim(min(noshift["energy"]),max(noshift["energy"]))
        plt.savefig("./figs/{}.uber{}.{}.pdf".format(project,a,b),bbox_inches='tight',pad_inches=0)


    # plt.show()
    uberfits=unfold_orbits(uberfits,"_d0")
    uberfits=unfold_orbits(uberfits,"_gamma")
    uberfits=unfold_orbits(uberfits,"_lambda")
    uberfits=unfold_orbits(uberfits,"_shift")
    uberfits=unfold_orbits(uberfits,"_tmp")

    A, B = fractional_coordinates(uberfits)

    E=uberfits["_gamma"].values
    ipolable={"a_frac":A.tolist(),"b_frac":B.tolist(),"values":E.tolist()}
    write_json(ipolable,"./ipol/gamma.json")

    E=uberfits["_d0"].values
    ipolable={"a_frac":A.tolist(),"b_frac":B.tolist(),"values":E.tolist()}
    write_json(ipolable,"./ipol/d0.json")


    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y, Z = bin_data_slice(uberfits, avec, bvec, "_gamma")
    Z=transform_energy_to_surface_energy(Z,avec,bvec)
    ax=plot_surface_heatmap(ax,X,Y,Z,avec,bvec,"x","y","z")
    plt.show()

if __name__ == "__main__":
    main()
