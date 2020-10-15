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

def load_record(record_file):
    with open(record_file) as f:
        records = json.load(f)
    return records

def write_json(data,target_file):
    with open(target_file,'w') as f:
        json.dump(data,f)
    return

def in_plane_vectors(record):
    a,b=record["grid"]
    au,bu=record["shift_units"]
    avec,bvec=np.array(au)*a,np.array(bu)*b

    return avec,bvec

def tabulize_record(record):
    df=pd.read_json(json.dumps(record["ids"]), orient="index")
    return df


def scatter_equivalent_shifts(ax,record,**kwargs):
    unwinded=tabulize_record(record)

    #Find shifts for any cleave entry (if you just ran "shift" you don't need this step)
    nocleave=unwinded.loc[unwinded["cleavage"]==unwinded["cleavage"].iloc[0]]

    #Get in plane shift vectors
    xy=np.array(nocleave["shift"].to_list())
    ab=np.array(nocleave["grid_point"].to_list())
    #We'll use their equivalent group value to color the plot
    c=nocleave["orbit"].to_numpy()

    #This is all we need, we save it to a single array
    points=np.hstack((ab,xy,c[:,np.newaxis]))

    #We'll duplicate the "walls" along a and b onto the next periodic image
    avec,bvec=in_plane_vectors(record)
    a,b=record["grid"]

    #Get the walls and repeat them on the next boundary
    aarg=points[:,1].argsort()[0:a]
    awall=points[aarg]
    awall[:,[2,3]]+=bvec
    points=np.vstack((points,awall))

    barg=points[:,0].argsort()[0:b+1]
    bwall=points[barg]
    bwall[:,[2,3]]+=avec
    points=np.vstack((points,bwall))

    #Scatter the grid points, coloring symmetrically equivalent ones with the same value
    ax.scatter(points[:,2],points[:,3],c=points[:,4],**kwargs) 

    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

def transform_energy_to_surface_energy(E,avec,bvec):
    area=np.linalg.norm(np.cross(avec,bvec))
    ev_to_joules=1.602176565e-19
    angstrom_to_meter=1e-10
    return ev_to_joules*E/(area*angstrom_to_meter*angstrom_to_meter)

def read_oszicar_energy(oszicar_path):
    with open(oszicar_path) as f:
        content = f.readlines()

    return float(content[-1].split()[2])

def surface_energy(unwinded):
    eqe = equilibrium_energy(unwinded)
    maxcleave = max(unwinded["cleavage"])

    surface_energies = (unwinded.loc[unwinded["cleavage"] ==
                                    maxcleave]["energy"] - eqe)*0.5
    print(surface_energies)
    exit()
    if max(surface_energies)-min(surface_energies)>0.001:
        raise ValueError("Surface energies don't match at different shifts. Did you separate enough?")

    return np.mean(surface_energies)

def equilibrium_energy(unwinded):
    print(unwinded["grid_point"])
    return float(unwinded.loc[(unwinded["grid_point"].apply(lambda x: x==[0,0])) & (unwinded["cleavage"] == 0.0)]["energy"])

def shift_energies_to_surface_energy(unwinded):
    eqe = equilibrium_energy(unwinded)
    gamma=surface_energy(unwinded)
    unwinded["energy"]-=(2*gamma+eqe)
    return unwinded


def uber(d, d0, gamma, lamb, shift):
    delta=d-d0
    return -2 * gamma * (1 + delta / lamb + 0) * np.exp(-delta / lamb) + shift


def uber_fit(unwinded_slice):
    unwinded_slice = unwinded_slice.sort_values(["cleavage"])
    sigma=np.ones(len(unwinded_slice))

    popt, pcov = curve_fit(uber, unwinded_slice["cleavage"], unwinded_slice["energy"],sigma=sigma,absolute_sigma=True)

    return popt,pcov

def plot_shifted_uber_fit(ax,unwinded_slice):
    unwinded_slice=shift_energies_to_surface_energy(unwinded_slice)
    popt, pcov = uber_fit(unwinded_slice)

    deltas = np.linspace(
        min(unwinded_slice["cleavage"]), max(unwinded_slice["cleavage"]), 1000)
    fitenergy = uber(deltas, *popt)
    
    ax.plot(deltas, fitenergy-popt[3], "g--", color='k',zorder=1)
    ax.scatter(unwinded_slice["cleavage"], unwinded_slice["energy"]-popt[3],zorder=2)

    return ax

def surface_data_reshape(X,Y,Z,amax,bmax):
    X = X.reshape(amax, bmax)
    Y = Y.reshape(amax, bmax)
    Z = Z.reshape(X.shape)

    return X,Y,Z

def plot_surface_heatmap(ax,X,Y,Z):
    pc=ax.pcolormesh(X, Y, Z)
    ax.set_aspect('equal')
    im_ratio=(max(Y.ravel())-min(Y.ravel()))/(max(X.ravel())-min(X.ravel()))
    cbar=plt.colorbar(pc,ax=ax,fraction=0.047*im_ratio,pad=0.04)

    return ax,cbar

def insert_vasp_energy(record):
    for orbit in record["equivalents"]:
        p=orbit[0]
        d=record["ids"][p]["directory"]
        energy=read_oszicar_energy(os.path.join(d,"OSZICAR")) 

        for i in orbit:
            record["ids"][i]["vasp"]=energy

    return record

def main():
    record=load_record("./record.json")
    record=insert_vasp_energy(record)
    unwinded=tabulize_record(record)

    energies=unwinded["vasp"]
    
    avec,bvec=in_plane_vectors(record)
    energies=[transform_energy_to_surface_energy(e,avec,bvec) for e in energies]
    unwinded["energy"]=energies

    fig=plt.figure()
    ax=fig.add_subplot(111)

    xy=np.array(unwinded["shift"].to_list())
    a,b=record["grid"]
    X,Y,Z=surface_data_reshape(xy[:,0],xy[:,1],unwinded["energy"].to_numpy(),a,b)
    Z=Z-(min(Z.ravel()))

    ax,cbar=plot_surface_heatmap(ax,X,Y,Z)

    ax.set_xlabel("x $\mathrm{[\AA]}$")
    ax.set_ylabel("y $\mathrm{[\AA]}$")
    cbar.set_label("Energy $\mathrm{[J/m^2]}$")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
