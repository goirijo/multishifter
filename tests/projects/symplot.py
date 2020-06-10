import json
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

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

    a = np.array([float(x) for x in aline.split()])
    b = np.array([float(x) for x in bline.split()])

    return a, b


def load_record(record_file):
    with open(record_file) as f:
        records = json.load(f)
    return records


def scatter_equivalent_shifts(record_file,slab_file):
    avec, bvec = load_slab_plane_vectors(slab_file)

    records = load_record(record_file)

    # shift_records = records["shift"]
    shift_records=records

    aix = [shift_records[s]["a_index"] for s in shift_records]
    bix = [shift_records[s]["b_index"] for s in shift_records]

    acount = len(set(aix))
    bcount = len(set(bix))

    eqs = [
        str(shift_records[s]["equivalent_shifts"]) for s in shift_records
    ]

    unq = set(eqs)
    colors = np.arange(len(unq)) / len(unq)
    cmap = {u: c for u, c in zip(list(unq), colors)}

    # print(len(unq),len(eqs))
    # exit()

    pts = np.array(
        [avec * a / acount + bvec * b / bcount for a, b in zip(aix, bix)])
    c = [cmap[eq] for eq in eqs]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.scatter(pts[:, 0], pts[:, 1], c=c)

    ax.set_aspect('equal')
    plt.show()


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

def bin_data_slice(unwinded, avec, bvec):
    L = np.array([avec[0:2], bvec[0:2]]).T


    amax, bmax = grid_divisions(unwinded)
    A, B = fractional_coordinates(unwinded)
    X, Y = L.dot(np.array([A, B]))

    X = X.reshape(amax, bmax) / amax    #already divided?
    Y = Y.reshape(amax, bmax) / bmax
    Z = np.array(unwinded["energy"]).reshape(X.shape)

    return X, Y, Z


def data_slice_orbits(unwinded):
    id_orbits = []
    for s in unwinded["equivalent_shifts"]:
        if s not in id_orbits:
            id_orbits.append(s)
    return id_orbits


def load_energies(unwinded):
    unwinded["energy"] = np.nan

    oszis = glob.glob("./Al-FCC.chain/shift__*/cleave__*/OSZICAR")
    energies = {f[2:-8]: load_oszicar_energy(f) for f in oszis}
    for p in energies:
        unwinded.loc[unwinded["path"] == p, "energy"] = energies[p]
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


def unfold_orbits(unwinded_slice):
    orbits = data_slice_orbits(unwinded_slice)

    for o in orbits:
        equivs=unwinded_slice.loc[unwinded_slice["equivalent_shifts"].astype(str)
                           == str(o), "energy"]
        not_nan=equivs[equivs.notnull()]
        assert(len(not_nan)>0)

        e = float(not_nan)
        unwinded_slice.loc[unwinded_slice["equivalent_shifts"].astype(str)
                           == str(o), "energy"] = e
    return unwinded_slice


def uber(delta, gamma, lamb, shift):
    return -2 * gamma * (1 + delta / lamb + 0) * np.exp(-delta / lamb) + shift


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

def plot_shifted_uber_fit(ax,unwinded_slice,gamma):
    ax.scatter(unwinded_slice["cleavage"], unwinded_slice["energy"])

    unwinded_slice = unwinded_slice.sort_values(["cleavage"])
    sigma=np.ones(len(unwinded_slice))
    sigma[0::10]=0.001

    partial=lambda d,g,l:uber(d,g,l,0.0)

    popt, pcov = curve_fit(partial, unwinded_slice["cleavage"], unwinded_slice["energy"],sigma=sigma,absolute_sigma=True)
    deltas = np.linspace(
        min(unwinded_slice["cleavage"]), max(unwinded_slice["cleavage"]), 1000)
    fitenergy = partial(deltas, *popt)
    ax.set_xlabel("$\delta$ $\mathrm{[\AA]}$")
    ax.set_ylabel("Energy [eV]")
    ax.plot(deltas, fitenergy, "g--")

    return ax

def plot_gamma_surface_heatmap(ax,unwinded_slice,avec,bvec):
    X, Y, Z = bin_data_slice(unwinded_slice, avec, bvec)

    pc=ax.pcolormesh(X, Y, Z)
    ax.set_aspect('equal')

    im_ratio=(max(Y.ravel())-min(Y.ravel()))/(max(X.ravel())-min(X.ravel()))
    cbar=plt.colorbar(pc,ax=ax,fraction=0.047*im_ratio,pad=0.04)
    cbar.set_label("Energy [eV]")

    ax.set_xlabel("x $\mathrm{[\AA]}$")
    ax.set_ylabel("y $\mathrm{[\AA]}$")
    return ax

def main():
    records = load_record("./Al-FCC.chain/record.json")
    # print(json.dumps(records["shift-cleave"],indent=4,sort_keys=True))

    unwinded = unwind_data(records)
    unwinded = load_energies(unwinded)

    avec, bvec = load_slab_plane_vectors("./Al-FCC.chain/slab.vasp")

    unfolded=[unfold_orbits(unwinded.loc[unwinded["cleavage"]==cleave].copy()) for cleave in set(unwinded["cleavage"])]
    unwinded=pd.concat(unfolded)

    unwinded = shift_energies_to_surface_energy(unwinded)
    gamma=surface_energy(unwinded)

    nocleave = unwinded.loc[unwinded["cleavage"] == 0.0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    Z1 = nocleave["energy"]
    plot_gamma_surface_heatmap(ax,nocleave,avec,bvec)
    plt.show()
    # plt.savefig("./figs/gamma0.0000.pdf",bbox_inches='tight',pad_inches=0)

    # A,B=fractional_coordinates(nocleave)
    # E=nocleave["energy"]

    # out={}
    # out["a_frac"]=list(A)
    # out["b_frac"]=list(B)
    # out["values"]=list(E)

    # with open('dump.json', 'w') as f:
    #     json.dump(out, f)

    # exit()

    # print(Z1 - Z2)
    # print("max error {} meV".format(max((Z1 - Z2).ravel()) * 1000))

    noshift = unwinded.loc[(unwinded["a_index"] == 0) &
                           (unwinded["b_index"] == 0)]

    id_orbits=data_slice_orbits(nocleave)
    prototypes=[o[0] for o in id_orbits]
    # unqa=[int(nocleave.loc[nocleave["id"]==p]["a_index"]) for p in prototypes]
    # unqb=[int(nocleave.loc[nocleave["id"]==p]["b_index"]) for p in prototypes]

    for a,b in prototypes:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        shiftdata=unwinded.loc[(unwinded["a_index"]==a)&(unwinded["b_index"]==b)]
        ax=plot_shifted_uber_fit(ax,shiftdata,gamma)
        plt.savefig("./figs/uber{}.{}.pdf".format(a,b),bbox_inches='tight',pad_inches=0)



if __name__ == "__main__":
    main()
