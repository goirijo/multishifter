import json
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

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
    c=nocleave["group"].to_numpy()

    #This is all we need, we save it to a single array
    points=np.hstack((ab,xy,c[:,np.newaxis]))

    #We'll duplicate the "walls" along a and b onto the next periodic image
    a,b=record["grid"]
    au,bu=record["shift_units"]
    avec,bvec=np.array(au)*a,np.array(bu)*b

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

def main():
    record=load_record("./record.json")

    fig=plt.figure()
    ax=fig.add_subplot(111)

    scatter_equivalent_shifts(ax,record,s=85,cmap=matplotlib.cm.get_cmap("twilight"))
    plt.show()

if __name__ == "__main__":
    main()
