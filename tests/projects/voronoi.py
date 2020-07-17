from reciprocals import *
from scipy.spatial import Voronoi

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

def main():
    L=load_lattice("./L.txt")
    Lt=load_lattice("./Lt.txt")
    K=load_lattice("./K.txt")
    Kt=load_lattice("./Kt.txt")
    G=load_lattice("./G.txt")
    Gz=load_lattice("./Gz.txt")
    Gzt=load_lattice("./Gzt.txt")
    M=load_lattice("./M.txt")
    Mt=load_lattice("./Mt.txt")

    fig=plt.figure(0,figsize=[15,15])
    ax=fig.add_subplot('231')
    ax.set_aspect("equal")

    plot_periodic_lattice_cells(ax,L,8,8,c='red')
    plot_lattice_unit_vectors(ax,L,width=0.1,color="red")
    plot_periodic_lattice_cells(ax,Lt,8,8,c='green')
    plot_lattice_unit_vectors(ax,Lt,width=0.1,color="green")
    # ax.set_ylim([-2,4])
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])

    ax=fig.add_subplot('232')
    ax.set_aspect("equal")

    plot_lattice_unit_vectors(ax,K,width=0.1,color="red")
    plot_lattice_unit_vectors(ax,Kt,width=0.1,color="green")

    plot_lattice_unit_vectors(ax,G,width=0.1,color="gray")
    plot_lattice_unit_vectors(ax,Gz,width=0.1,color="k")
    plot_periodic_lattice_cells(ax,K,8,8,c='red')
    plot_voronoi_lattice(ax,K,8,8,c='k')
    # ax.set_ylim([-3,3])
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])

    ax=fig.add_subplot('233')
    ax.set_aspect("equal")

    plot_lattice_unit_vectors(ax,K,width=0.1,color="red")
    plot_lattice_unit_vectors(ax,Kt,width=0.1,color="green")

    plot_lattice_unit_vectors(ax,G,width=0.1,color="gray")
    plot_lattice_unit_vectors(ax,Gzt,width=0.1,color="k")
    plot_periodic_lattice_cells(ax,Kt,8,8,c='green')
    plot_voronoi_lattice(ax,Kt,8,8,c='k')
    ax.set_ylim([-6,6])
    ax.set_xlim([-6,6])

    ax=fig.add_subplot('235')
    ax.set_aspect("equal")

    plot_periodic_lattice_cells(ax,L,150,150,c='red')
    plot_lattice_unit_vectors(ax,L,width=0.1,color="red")
    plot_periodic_lattice_cells(ax,Lt,150,150,c='green')
    plot_lattice_unit_vectors(ax,Lt,width=0.1,color="green")
    plot_lattice_unit_vectors(ax,M,width=0.5,color="k")
    plot_periodic_lattice_cells(ax,M,6,6,lw=2,c='k')
    ax.set_ylim([-50,50])
    ax.set_xlim([-50,50])

    ax=fig.add_subplot('236')
    ax.set_aspect("equal")

    plot_periodic_lattice_cells(ax,L,150,150,c='red')
    plot_lattice_unit_vectors(ax,L,width=0.1,color="red")
    plot_periodic_lattice_cells(ax,Lt,150,150,c='green')
    plot_lattice_unit_vectors(ax,Lt,width=0.1,color="green")
    plot_lattice_unit_vectors(ax,Mt,width=0.5,color="k")
    plot_periodic_lattice_cells(ax,Mt,6,6,lw=2,c='k')
    ax.set_ylim([-50,50])
    ax.set_xlim([-50,50])

    print(np.linalg.det(M.column_vector_matrix()),np.linalg.det(Mt.column_vector_matrix()))
    plt.savefig("./tmp.png",bbox_inches='tight',dpi=100)


if __name__=="__main__":
    main()
