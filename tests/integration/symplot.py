import json
import numpy as np
import matplotlib.pyplot as plt

def load_slab_plane_vectors(slabfile):
    with open(slabfile) as f:
        content = f.readlines()

    aline=content[2]
    bline=content[3]

    a=np.array([float(x) for x in aline.split()])
    b=np.array([float(x) for x in bline.split()])

    return a,b

def main():
    avec,bvec=load_slab_plane_vectors("./slab.vasp")

    with open("./record.json") as f:
        records = json.load(f)

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


if __name__=="__main__":
    main()
