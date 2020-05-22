import json
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

    aline=content[2]
    bline=content[3]

    a=np.array([float(x) for x in aline.split()])
    b=np.array([float(x) for x in bline.split()])

    return a,b

def load_record():
    with open("./Al-FCC.chain/record.json") as f:
        records = json.load(f)
    return records

def scatter_equivalent_shifts():
    avec,bvec=load_slab_plane_vectors("./Al-FCC.chain/slab.vasp")

    records=load_record()

    shift_records=records["shift"]

    aix=[shift_records[s]["a_index"] for s in shift_records]
    bix=[shift_records[s]["b_index"] for s in shift_records]

    acount=len(set(aix))
    bcount=len(set(bix))

    eqs=[str(shift_records[s]["equivalent_structures"]) for s in shift_records]

    unq=set(eqs)
    colors=np.arange(len(unq))/len(unq)
    cmap={u:c for u,c in zip(list(unq),colors)}

    fig=plt.figure()
    ax=fig.add_subplot(111)

    pts=np.array([avec*a/acount+bvec*b/bcount for a,b in zip(aix,bix)])
    c=[cmap[eq] for eq in eqs]

    plt.scatter(pts[:,0],pts[:,1],c=c)

    ax.set_aspect('equal')
    plt.show()

def load_oszicar_energy(oszicar_path):
    with open(oszicar_path) as f:
        content = f.readlines()

    return float(content[-1].split()[2])

def unwind_data(records):
    paths=[p for p in records["shift-cleave"]]
    avals=[records["shift-cleave"][p]["a_index"] for p in records["shift-cleave"]]
    bvals=[records["shift-cleave"][p]["b_index"] for p in records["shift-cleave"]]
    angles=[records["shift-cleave"][p]["angle"] for p in records["shift-cleave"]]
    ids=[records["shift-cleave"][p]["id"] for p in records["shift-cleave"]]
    cleaves=[records["shift-cleave"][p]["cleavage"] for p in records["shift-cleave"]]
    equivs=[records["shift-cleave"][p]["equivalent_structures"] for p in records["shift-cleave"]]

    unwinded=pd.DataFrame.from_dict({"path":paths,"a_index":avals,"b_index":bvals,"angle":angles,"id":ids,"cleavage":cleaves,"equivalent_structures":equivs})

    unwinded=unwinded.sort_values(["a_index","b_index"])
    return unwinded

def bin_data_slice(unwinded,avec,bvec):
    L=np.array([avec[0:2],bvec[0:2]]).T

    amax,bmax=max(unwinded["a_index"])+1,max(unwinded["b_index"])+1
    A=np.array(unwinded["a_index"])/amax
    B=np.array(unwinded["b_index"])/bmax
    X,Y=L.dot(np.array([A,B]))

    X=X.reshape(amax,bmax)/amax
    Y=Y.reshape(amax,bmax)/bmax
    Z=np.array(unwinded["energy"]).reshape(X.shape)

    return X,Y,Z

def data_slice_orbits(unwinded):
    id_orbits=[]
    for s in unwinded["equivalent_structures"]:
        if s not in id_orbits:
            id_orbits.append(s)
    return id_orbits

def main():
    records=load_record()
    # print(json.dumps(records["shift-cleave"],indent=4,sort_keys=True))

    unwinded=unwind_data(records)
    unwinded["energy"]=0.0;

    oszis=glob.glob("./Al-FCC.chain/shift__*/cleave__*/OSZICAR")
    energies={f[2:-8]:load_oszicar_energy(f) for f in oszis}
    for p in energies:
        unwinded.loc[unwinded["path"]==p,"energy"]=energies[p]


    avec,bvec=load_slab_plane_vectors("./Al-FCC.chain/slab.vasp")

    nocleave=unwinded.loc[unwinded["cleavage"]==0.0]
    X,Y,Z=bin_data_slice(nocleave,avec,bvec)
    X1,Y1,Z1=X,Y,Z

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    ax.plot_surface(X,Y,Z)
    ax=set_axes_equal(ax)
    plt.show()


    orbits=data_slice_orbits(nocleave)

    for o in orbits:
        assert(len(nocleave.loc[nocleave["id"]==o[0]])==1)
        e=float(nocleave.loc[nocleave["id"]==o[0]]["energy"])
        nocleave.loc[nocleave["equivalent_structures"].astype(str)==str(o),"energy"]=e
    
    X,Y,Z=bin_data_slice(nocleave,avec,bvec)
    X2,Y2,Z2=X,Y,Z

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    ax.plot_surface(X,Y,Z)
    ax=set_axes_equal(ax)
    plt.show()

    print(Z1-Z2)
    print("max error {}".format(max((Z1-Z2).ravel())*1000))

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    ax.plot_surface(X2,Y2,Z1-Z2)
    ax=set_axes_equal(ax)
    plt.show()



if __name__=="__main__":
    main()
