import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
import pandas as pd

font = {'family' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)


def surface(x,y):
    return 0+1.983616e-08*np.sin(1.968300e+00*x-3.409195e+00*y)-3.608439e-08*np.sin(3.936600e+00*x-6.818391e+00*y)-1.111111e-07*np.sin(5.904900e+00*x-1.022759e+01*y)+1.556357e-07*np.cos(9.841499e+00*x-7.954789e+00*y)-1.726139e-07*np.sin(9.841499e+00*x-1.704598e+01*y)-2.766470e-07*np.sin(7.873199e+00*x-1.363678e+01*y)+3.721421e-07*np.cos(1.968300e+00*x-1.250038e+01*y)+7.216878e-07*np.sin(7.873199e+00*x-1.818238e+01*y)+7.216878e-07*np.sin(7.873199e+00*x+9.091188e+00*y)+1.106588e-06*np.sin(1.180980e+01*x+2.272797e+00*y)-1.106588e-06*np.sin(1.180980e+01*x-1.590958e+01*y)+1.315616e-06*np.cos(9.841499e+00*x-1.022759e+01*y)+1.386423e-06*np.cos(3.936600e+00*x-1.363678e+01*y)+2.048611e-06*np.cos(7.873199e+00*x-1.363678e+01*y)-2.701389e-06*np.cos(3.936600e+00*x-1.136398e+01*y)-2.854167e-06*np.cos(7.873199e+00*x-9.091188e+00*y)-3.548294e-06*np.cos(7.873199e+00*x+2.272797e+00*y)-3.582518e-06*np.cos(5.904900e+00*x+5.681992e+00*y)+4.822805e-06*np.cos(7.873199e+00*x-6.818391e+00*y)+5.025058e-06*np.cos(1.968300e+00*x-1.022759e+01*y)-5.184124e-06*np.sin(7.873199e+00*x+4.545594e+00*y)+7.930556e-06*np.cos(1.180980e+01*x-1.136398e+01*y)+7.930556e-06*np.cos(1.180980e+01*x-2.272797e+00*y)+8.013889e-06*np.cos(3.936600e+00*x+1.136398e+01*y)+8.013889e-06*np.cos(3.936600e+00*x-1.590958e+01*y)+8.151588e-06*np.sin(9.841499e+00*x-1.022759e+01*y)-8.277218e-06*np.cos(5.904900e+00*x+7.954789e+00*y)-8.283208e-06*np.cos(9.841499e+00*x-1.477318e+01*y)-8.293011e-06*np.cos(7.873199e+00*x-1.590958e+01*y)-8.299835e-06*np.cos(9.841499e+00*x-1.931877e+01*y)-8.299835e-06*np.cos(9.841499e+00*x+7.954789e+00*y)-8.349941e-06*np.cos(1.180980e+01*x-1.818238e+01*y)-8.349941e-06*np.cos(1.180980e+01*x+4.545594e+00*y)-8.381618e-06*np.sin(3.936600e+00*x-1.363678e+01*y)-8.391053e-06*np.cos(9.841499e+00*x+1.136398e+00*y)+9.011452e-06*np.cos(9.841499e+00*x-1.136398e+00*y)+9.154428e-06*np.cos(3.936600e+00*x+9.091188e+00*y)+9.182453e-06*np.sin(7.873199e+00*x+6.818391e+00*y)+9.343039e-06*np.sin(9.841499e+00*x+3.409195e+00*y)+1.014881e-05*np.sin(7.873199e+00*x-1.136398e+01*y)-1.042559e-05*np.sin(5.904900e+00*x-1.250038e+01*y)+1.053201e-05*np.sin(7.873199e+00*x-1.590958e+01*y)-1.089721e-05*np.sin(9.841499e+00*x-1.477318e+01*y)+1.138030e-05*np.sin(9.841499e+00*x-1.250038e+01*y)-1.180514e-05*np.sin(5.904900e+00*x-1.477318e+01*y)+1.191988e-05*np.sin(7.873199e+00*x-9.091188e+00*y)-1.208827e-05*np.sin(3.936600e+00*x-1.136398e+01*y)+1.368121e-05*np.cos(5.904900e+00*x+1.136398e+00*y)-1.376663e-05*np.sin(9.841499e+00*x-1.136398e+00*y)+1.380893e-05*np.cos(3.936600e+00*x+4.545594e+00*y)-1.401427e-05*np.sin(3.936600e+00*x+9.091188e+00*y)-1.451389e-05*np.cos(1.180980e+01*x+2.272797e+00*y)-1.451389e-05*np.cos(1.180980e+01*x-1.590958e+01*y)-1.454167e-05*np.cos(7.873199e+00*x-1.818238e+01*y)-1.454167e-05*np.cos(7.873199e+00*x+9.091188e+00*y)+1.496299e-05*np.sin(1.180980e+01*x-1.136398e+01*y)-1.496299e-05*np.sin(1.180980e+01*x-2.272797e+00*y)-1.525167e-05*np.sin(3.936600e+00*x+1.136398e+01*y)-1.525167e-05*np.sin(3.936600e+00*x-1.590958e+01*y)+1.564720e-05*np.cos(9.841499e+00*x+5.681992e+00*y)-1.581699e-05*np.sin(7.873199e+00*x-2.554503e-08*y)-1.600944e-05*np.sin(3.936600e+00*x+6.818391e+00*y)+1.721528e-05*np.cos(7.873199e+00*x+4.545594e+00*y)-1.852415e-05*np.cos(9.841499e+00*x-1.704598e+01*y)+1.888889e-05*np.cos(1.180980e+01*x-2.045517e+01*y)+1.888889e-05*np.cos(1.180980e+01*x+6.818391e+00*y)-1.981944e-05*np.cos(5.904900e+00*x-1.704598e+01*y)-1.981944e-05*np.cos(5.904900e+00*x+1.022759e+01*y)-1.993056e-05*np.cos(1.180980e+01*x-1.363678e+01*y)-1.993056e-05*np.cos(1.180980e+01*x-3.831754e-08*y)+2.031944e-05*np.sin(1.180980e+01*x-1.363678e+01*y)-2.031944e-05*np.sin(1.180980e+01*x-3.831754e-08*y)-2.056944e-05*np.sin(5.904900e+00*x-1.704598e+01*y)-2.056944e-05*np.sin(5.904900e+00*x+1.022759e+01*y)+2.203213e-05*np.sin(1.180980e+01*x-1.818238e+01*y)-2.203213e-05*np.sin(1.180980e+01*x+4.545594e+00*y)-2.215704e-05*np.sin(9.841499e+00*x+1.136398e+00*y)-2.224691e-05*np.sin(9.841499e+00*x-1.931877e+01*y)-2.224691e-05*np.sin(9.841499e+00*x+7.954789e+00*y)-2.240505e-05*np.sin(5.904900e+00*x+7.954789e+00*y)-2.707010e-05*np.cos(5.904900e+00*x-7.954789e+00*y)-2.708671e-05*np.cos(3.936600e+00*x-9.091188e+00*y)-2.733061e-05*np.cos(9.841499e+00*x+3.409195e+00*y)-2.736447e-05*np.cos(7.873199e+00*x+6.818391e+00*y)+2.869563e-05*np.sin(1.968300e+00*x+3.409195e+00*y)+2.875727e-05*np.sin(3.936600e+00*x-1.277251e-08*y)-2.888970e-05*np.cos(5.904900e+00*x-1.250038e+01*y)-2.893617e-05*np.cos(5.904900e+00*x-1.477318e+01*y)-2.895171e-05*np.cos(7.873199e+00*x-1.136398e+01*y)-2.898428e-05*np.cos(9.841499e+00*x-1.250038e+01*y)+2.905556e-05*np.cos(5.904900e+00*x+3.409195e+00*y)+2.932216e-05*np.cos(1.180980e+01*x-9.091188e+00*y)+2.932216e-05*np.cos(1.180980e+01*x-4.545594e+00*y)+2.966095e-05*np.cos(1.968300e+00*x-1.477318e+01*y)+2.966095e-05*np.cos(1.968300e+00*x+1.250038e+01*y)+3.086806e-05*np.cos(3.936600e+00*x-6.818391e+00*y)+3.087482e-05*np.sin(5.904900e+00*x-7.954789e+00*y)-3.094699e-05*np.sin(3.936600e+00*x-9.091188e+00*y)+3.277083e-05*np.cos(7.873199e+00*x-2.554503e-08*y)+3.284028e-05*np.cos(3.936600e+00*x+6.818391e+00*y)+4.283079e-05*np.sin(3.936600e+00*x+4.545594e+00*y)+4.290295e-05*np.sin(5.904900e+00*x+1.136398e+00*y)-4.438194e-05*np.cos(3.936600e+00*x+2.272797e+00*y)-4.986111e-05*np.cos(5.904900e+00*x-1.022759e+01*y)+5.115082e-05*np.sin(9.841499e+00*x+5.681992e+00*y)-5.617782e-05*np.sin(1.968300e+00*x-5.681992e+00*y)+5.630518e-05*np.sin(3.936600e+00*x-4.545594e+00*y)-5.696230e-05*np.sin(7.873199e+00*x+2.272797e+00*y)-5.704663e-05*np.sin(5.904900e+00*x+5.681992e+00*y)+7.051303e-05*np.cos(9.841499e+00*x-3.409195e+00*y)+7.068142e-05*np.cos(1.968300e+00*x+1.022759e+01*y)-8.079217e-05*np.sin(9.841499e+00*x-5.681992e+00*y)+8.083011e-05*np.sin(0.000000e+00*x-1.136398e+01*y)-8.157127e-05*np.cos(1.968300e+00*x-7.954789e+00*y)-8.161938e-05*np.cos(5.904900e+00*x-5.681992e+00*y)-8.573188e-05*np.sin(7.873199e+00*x-2.272797e+00*y)-8.588057e-05*np.sin(1.968300e+00*x+7.954789e+00*y)-9.856717e-05*np.sin(5.904900e+00*x-1.136398e+00*y)-9.870407e-05*np.sin(1.968300e+00*x+5.681992e+00*y)+9.898085e-05*np.sin(9.841499e+00*x-7.954789e+00*y)-9.905029e-05*np.sin(1.968300e+00*x-1.250038e+01*y)-1.013587e-04*np.cos(1.968300e+00*x-5.681992e+00*y)-1.014322e-04*np.cos(3.936600e+00*x-4.545594e+00*y)-1.105181e-04*np.sin(9.841499e+00*x-3.409195e+00*y)-1.106624e-04*np.sin(1.968300e+00*x+1.022759e+01*y)+1.107380e-04*np.sin(7.873199e+00*x-6.818391e+00*y)-1.107458e-04*np.sin(1.968300e+00*x-1.022759e+01*y)-1.113290e-04*np.sin(1.180980e+01*x-4.545594e+00*y)+1.113290e-04*np.sin(1.180980e+01*x-9.091188e+00*y)-1.114475e-04*np.sin(1.968300e+00*x-1.477318e+01*y)-1.114475e-04*np.sin(1.968300e+00*x+1.250038e+01*y)+1.231263e-04*np.cos(7.873199e+00*x-2.272797e+00*y)+1.233527e-04*np.cos(1.968300e+00*x+7.954789e+00*y)-1.343525e-04*np.sin(1.968300e+00*x-7.954789e+00*y)+1.343607e-04*np.sin(5.904900e+00*x-5.681992e+00*y)-1.370556e-04*np.sin(5.904900e+00*x+3.409195e+00*y)-1.421244e-04*np.sin(7.873199e+00*x-4.545594e+00*y)+1.423168e-04*np.sin(0.000000e+00*x-9.091188e+00*y)+1.808050e-04*np.cos(5.904900e+00*x-1.136398e+00*y)+1.808355e-04*np.cos(1.968300e+00*x+5.681992e+00*y)-1.882083e-04*np.cos(0.000000e+00*x-6.818391e+00*y)-1.883194e-04*np.cos(5.904900e+00*x-3.409195e+00*y)-2.687500e-04*np.cos(0.000000e+00*x-1.363678e+01*y)-2.689722e-04*np.cos(1.180980e+01*x-6.818391e+00*y)-2.839446e-04*np.cos(0.000000e+00*x-1.136398e+01*y)-2.841824e-04*np.cos(9.841499e+00*x-5.681992e+00*y)-3.074861e-04*np.sin(5.904900e+00*x-3.409195e+00*y)+3.075139e-04*np.sin(0.000000e+00*x-6.818391e+00*y)-3.488750e-04*np.cos(0.000000e+00*x-4.545594e+00*y)-3.490139e-04*np.cos(3.936600e+00*x-2.272797e+00*y)-3.760139e-04*np.cos(0.000000e+00*x-9.091188e+00*y)-3.761806e-04*np.cos(7.873199e+00*x-4.545594e+00*y)+4.939833e-04*np.sin(3.936600e+00*x+2.272797e+00*y)-7.472837e-04*np.sin(3.936600e+00*x-2.272797e+00*y)+7.474280e-04*np.sin(0.000000e+00*x-4.545594e+00*y)+9.441353e-04*np.cos(1.968300e+00*x-3.409195e+00*y)+1.756627e-03*np.cos(3.936600e+00*x-1.277251e-08*y)+1.756629e-03*np.cos(1.968300e+00*x+3.409195e+00*y)-5.233629e-02*np.cos(1.968300e+00*x+1.136398e+00*y)-5.335879e-02*np.cos(1.968300e+00*x-1.136398e+00*y)-5.335889e-02*np.cos(0.000000e+00*x-2.272797e+00*y)+8.207326e-02*np.sin(1.968300e+00*x-1.136398e+00*y)-8.207327e-02*np.sin(0.000000e+00*x-2.272797e+00*y)-8.337125e-02*np.sin(1.968300e+00*x+1.136398e+00*y)-1.261405e+01*np.cos(0.000000e+00*x+0.000000e+00*y)

def plot_surface_heatmap(ax,X,Y,Z):
    pc=ax.pcolormesh(X, Y, Z)
    ax.set_aspect('equal')
    im_ratio=(max(Y.ravel())-min(Y.ravel()))/(max(X.ravel())-min(X.ravel()))
    cbar=plt.colorbar(pc,ax=ax,fraction=0.047*im_ratio,pad=0.04)

    return ax,cbar

def in_plane_vectors(record):
    a,b=record["grid"]
    au,bu=record["shift_units"]
    avec,bvec=np.array(au)*a,np.array(bu)*b

    return avec,bvec

def transform_energy_to_surface_energy(E,avec,bvec):
    area=np.linalg.norm(np.cross(avec,bvec))
    ev_to_joules=1.602176565e-19
    angstrom_to_meter=1e-10
    return ev_to_joules*E/(area*angstrom_to_meter*angstrom_to_meter)

def load_record(record_file):
    with open(record_file) as f:
        records = json.load(f)
    return records

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

def main():
    record=load_record("./modified_record.json")
    avec,bvec=in_plane_vectors(record)
    
    fig=plt.figure()
    ax=fig.add_subplot(111)

    x=np.arange(0,7,0.01)
    y=np.arange(0,7,0.01)
    X,Y=np.meshgrid(x,y)
    Z=transform_energy_to_surface_energy(surface(X,Y),avec,bvec)/2.0
    Z=Z-min(Z.ravel())

    ax,cbar=plot_surface_heatmap(ax,X,Y,Z)
    ax.set_xlabel("x $\mathrm{[\AA]}$")
    ax.set_ylabel("y $\mathrm{[\AA]}$")
    cbar.set_label("Energy $\mathrm{[J/m^2]}$")

    # scatter_equivalent_shifts(ax,record,s=85,cmap=matplotlib.cm.get_cmap("twilight"))
    plt.tight_layout()
    plt.show()

if __name__=="__main__":
    main()
