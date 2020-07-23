from reciprocals import *
from scipy.spatial import Voronoi
import glob
import os
import sys

def load_lattice(filename):
    with open(filename,'r') as filestream:
        data=filestream.read().replace('\n',' ')
    entries=data.split()
    entries=[float(e) for e in entries]
    return Lattice([entries[0],entries[3]],[entries[1],entries[4]])

def voronoi_plot(ax,vor,**kwargs):
    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            ax.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], **kwargs)
    return ax

def plot_voronoi_lattice(ax,L,arad,brad,**kwargs):
    col_matrix=L.column_vector_matrix()

    points=[]
    for i in range(-arad,arad):
        for j in range(-brad,brad):
            frac=np.array([i,j])
            cart=col_matrix.dot(frac)
            points.append(cart)
    vor=Voronoi(points)
    voronoi_plot(ax,vor,**kwargs)
    return ax

def plot_lattice_overlay(ax, L, Lt):
    plot_periodic_lattice_cells(ax,L,8,8,c='red')
    plot_lattice_unit_vectors(ax,L,width=0.1,color="red")
    plot_periodic_lattice_cells(ax,Lt,8,8,c='green')
    plot_lattice_unit_vectors(ax,Lt,width=0.1,color="green")
    # ax.set_ylim([-2,4])
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])
    return ax

def plot_reciprocal_voronoi(ax, K, Kt, G, Gz, V, Vcolor):
    plot_lattice_unit_vectors(ax,K,width=0.1,color="red")
    plot_lattice_unit_vectors(ax,Kt,width=0.1,color="green")
    plot_lattice_unit_vectors(ax,G,width=0.1,color="gray")
    plot_lattice_unit_vectors(ax,Gz,width=0.1,color="k")
    plot_periodic_lattice_cells(ax,V,8,8,c=Vcolor)
    plot_voronoi_lattice(ax,V,8,8,c='k')
    # ax.set_ylim([-3,3])
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])

    return ax

def plot_moire(ax,L,Lt,M):
    plot_periodic_lattice_cells(ax,L,150,150,c='red')
    plot_lattice_unit_vectors(ax,L,width=0.1,color="red")
    plot_periodic_lattice_cells(ax,Lt,150,150,c='green')
    plot_lattice_unit_vectors(ax,Lt,width=0.1,color="green")
    plot_lattice_unit_vectors(ax,M,width=0.5,color="k")
    plot_periodic_lattice_cells(ax,M,6,6,lw=2,c='k')
    ax.set_ylim([-50,50])
    ax.set_xlim([-50,50])

    return ax

def main1():
    dirs=glob.glob("./frames/1???")

    Ms=[load_lattice(os.path.join(d,"M.txt")) for d in dirs]
    Mts=[load_lattice(os.path.join(d,"Mt.txt")) for d in dirs]
    Mds=[load_lattice(os.path.join(d,"Md.txt")) for d in dirs]
    
    def deg(d):
        f=open(os.path.join(d,"degrees.txt"),'r')
        ls=f.readlines()
        return[float(ls[0]),bool(int(ls[1])),bool(int(ls[2])),bool(int(ls[3])),bool(int(ls[4]))]
    
    extras=[deg(d) for d in dirs]

    ds=[e[0] for e in extras]

    GzainKt=[e[1] for e in extras]
    GzbinKt=[e[2] for e in extras]
    GztainK=[e[3] for e in extras]
    GztbinK=[e[4] for e in extras]

    # Msvol=[abs(np.linalg.det(M.column_vector_matrix())) for M in Ms]
    # Mtsvol=[abs(np.linalg.det(M.column_vector_matrix())) for M in Mts]
    Msvol=[np.linalg.det(M.column_vector_matrix()) for M in Ms]
    Mtsvol=[np.linalg.det(M.column_vector_matrix()) for M in Mts]
    Mdsvol=[np.linalg.det(M.column_vector_matrix()) for M in Mds]

    Msvol=np.array(Msvol)
    Mtsvol=np.array(Mtsvol)

    Mcolor=["red" if a and b else "orange" for a,b in zip(GzainKt,GzbinKt) ]
    Mtcolor=["green" if a and b else "blue" for a,b in zip(GztainK,GztbinK) ]
    diffcolor=["gray" if a and b and c and d else "black" for _,a,b,c,d in extras ]

    fig=plt.figure(0)
    ax=fig.add_subplot('311')
    ax.scatter(ds,Msvol,c=Mcolor)
    ax.set_ylim([-800,800])

    ax=fig.add_subplot('312')
    ax.scatter(ds,Mtsvol,c=Mtcolor)
    ax.set_ylim([-800,800])

    ax=fig.add_subplot('313')
    ax.scatter(ds,abs(Msvol-Mtsvol),c=diffcolor)
    ax.set_ylim([-800,800])
    plt.show()


def main2():
    L=load_lattice("./frames/{}/L.txt".format(sys.argv[1]))
    Lt=load_lattice("./frames/{}/Lt.txt".format(sys.argv[1]))
    K=load_lattice("./frames/{}/K.txt".format(sys.argv[1]))
    Kt=load_lattice("./frames/{}/Kt.txt".format(sys.argv[1]))
    G=load_lattice("./frames/{}/G.txt".format(sys.argv[1]))
    Gz=load_lattice("./frames/{}/Gz.txt".format(sys.argv[1]))
    Gzt=load_lattice("./frames/{}/Gzt.txt".format(sys.argv[1]))
    Gzx=load_lattice("./frames/{}/Gzx.txt".format(sys.argv[1]))
    Gzx2=load_lattice("./frames/{}/Gzx2.txt".format(sys.argv[1]))
    Gd=load_lattice("./frames/{}/Gd.txt".format(sys.argv[1]))
    M=load_lattice("./frames/{}/M.txt".format(sys.argv[1]))
    Mt=load_lattice("./frames/{}/Mt.txt".format(sys.argv[1]))
    Mx=load_lattice("./frames/{}/Mx.txt".format(sys.argv[1]))
    Mx2=load_lattice("./frames/{}/Mx2.txt".format(sys.argv[1]))
    Md=load_lattice("./frames/{}/Md.txt".format(sys.argv[1]))

    fig=plt.figure(0,figsize=[15,10])

    ax=fig.add_subplot('231')
    ax.set_aspect("equal")
    plot_lattice_overlay(ax,L,Lt)

    ax=fig.add_subplot('232')
    ax.set_aspect("equal")
    plot_reciprocal_voronoi(ax,K,Kt,G,Gz,K,"red")

    ax=fig.add_subplot('233')
    ax.set_aspect("equal")
    plot_reciprocal_voronoi(ax,K,Kt,G,Gzt,Kt,"green")

    ax=fig.add_subplot('234')
    ax.set_aspect("equal")
    plot_voronoi_lattice(ax,K,8,8,c='red')
    plot_voronoi_lattice(ax,Kt,8,8,c='green')
    plot_lattice_unit_vectors(ax,G,width=0.1,color="gray")
    plot_lattice_unit_vectors(ax,Gz,width=0.1,color="red")
    plot_lattice_unit_vectors(ax,Gzt,width=0.1,color="green")
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])

    ax=fig.add_subplot('235')
    ax.set_aspect("equal")
    plot_moire(ax,L,Lt,M)

    ax=fig.add_subplot('236')
    ax.set_aspect("equal")
    plot_moire(ax,L,Lt,Mt)

    print(np.linalg.inv(M.column_vector_matrix()).dot(Mt.column_vector_matrix()))

    fig.tight_layout()
    plt.show()
    exit()
    print(np.linalg.det(M.column_vector_matrix()),np.linalg.det(Mt.column_vector_matrix()))
    plt.savefig("./vorimation/img{}.png".format(sys.argv[1]),bbox_inches='tight',dpi=100)


if __name__=="__main__":
    main1()
    main2()
