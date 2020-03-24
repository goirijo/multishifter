import math
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines


def reciprocal_vectors(a,b):
    """Construct a reciprocal lattice from the given
    lattice vectors

    :a: 2x1 lattice vector
    :b: 2x1 lattice vector
    :returns: vector,vector

    """
    real=np.array([a,b]).T
    recip=2*math.pi*np.linalg.inv(real).T

    astar=recip[:,0]
    bstar=recip[:,1]

    return astar,bstar

def make_rotation_matrix(angle):
    """Given an angle, get a 2D rotation matrix

    Parameters
    ----------
    angle : float

    Returns
    -------
    2x2 matrix

    """
    rad=np.radians(angle)
    c, s = np.cos(rad), np.sin(rad)
    return np.array(((c, -s), (s, c)))

class Lattice(object):

    """Docstring for Lattice. """

    def __init__(self,a,b):
        """2D lattice

        Parameters
        ----------
        a : 2x1 vector
        b : 2x1 vector


        """
        self._a = np.array(a)
        self._b = np.array(b)

        self._column_matrix=np.array([a,b]).T

    def a(self):
        """Get the first vector of the lattice
        Returns
        -------
        2x1 vector

        """
        return self._a

    def b(self):
        """Get the second vector of the lattice
        Returns
        -------
        2x1 vector

        """
        return self._b

    def column_vector_matrix(self):
        """Get both vectors in a matrix
        Returns
        -------
        2x2 matrix

        """
        return self._column_matrix

    def reciprocal(self):
        """Construct a new lattice that is reciprocal to self
        Returns
        -------
        Lattice

        """
        astar,bstar=reciprocal_vectors(self.a(),self.b())
        return Lattice(astar,bstar)

    def transform(self, T):
        """Apply a transformation matrix to the lattice vectors
        by multiplying on the left of the vector matrix, with
        lattice vectors as columns

        Parameters
        ----------
        T : 2x2 matrix

        Returns
        -------
        Lattice

        """
        transformed_column_matrix=T.dot(self.column_vector_matrix())
        return Lattice(transformed_column_matrix[:,0],transformed_column_matrix[:,1])


def scatter_lattice_points(ax,lattice: Lattice,arad,brad,**kwargs):
    """Scatter a bunch of lattice points around the origin to visualize
    the grid

    Parameters
    ----------
    ax : pyplot axes
    lattice : Lattice you want to visualize
    arad : maximum range along a from origin
    brad : maximum range along a from origin
    kwargs : plotting options

    Returns
    -------
    pyplot axes

    """
    col_matrix=lattice.column_vector_matrix()
    for i in range(-arad,arad):
        for j in range(-brad,brad):
            frac=np.array([i,j])
            cart=col_matrix.dot(frac)

            ax.scatter(cart[0],cart[1],**kwargs)

def plot_periodic_lattice_cells(ax,lattice: Lattice, arad, brad, **kwargs):
    """TODO: Docstring for plot_lattice_cells.

    Parameters
    ----------
    ax : pyplot axes
    lattice : Lattice you want to visualize
    arad : maximum range along a from origin
    brad : maximum range along a from origin
    kwargs : plotting options

    Returns
    -------
    pyplot axes

    """
    col_matrix=lattice.column_vector_matrix()

    def plot_line_from_fractional_points(lattice,frac0,frac1):
        cart0=col_matrix.dot(frac0)
        cart1=col_matrix.dot(frac1)
        cart=np.array([cart0,cart1])
        ax.plot(cart[:,0],cart[:,1],**kwargs)
        return

    for i in range(-arad,arad):
        frac0=np.array([i,brad])
        frac1=np.array([i,-brad])
        plot_line_from_fractional_points(lattice,frac0,frac1)

    for j in range(-brad,brad):
        frac0=np.array([arad,j])
        frac1=np.array([-arad,j])
        plot_line_from_fractional_points(lattice,frac0,frac1)
    
    return ax

def plot_lattice_unit_vectors(ax,lattice:Lattice,offset=[0,0],**kwargs):
    """Plot a single lattice cell at the offset with arrows

    Parameters
    ----------
    ax : pyplot axes
    lattice : Lattice to visualize
    **kwargs : plotting options
    offset : where to place the origin, optional

    Returns
    -------
    pyplot axes

    """
    x,y=offset
    a_x,a_y=lattice.a()
    b_x,b_y=lattice.b()

    ax.add_patch(patches.FancyArrow(x,y,a_x,a_y,length_includes_head=True,**kwargs))
    ax.add_patch(patches.FancyArrow(x,y,b_x,b_y,length_includes_head=True,**kwargs))

    return ax

def make_hexagonal_lattice():
    a=np.array([0,3])
    R=make_rotation_matrix(60)
    b=R.dot(a)
    return Lattice(a,b)

def make_square_lattice():
    a=[3,0]
    b=[0,3]
    return Lattice(a,b)

def _main():
    fig=plt.figure()
    ax=fig.add_subplot('111')

    x=np.arange(0,500,0.1)

    kbig=0.9
    ybig=np.sin(kbig*x)

    ksmall=0.85
    ysmall=np.sin(ksmall*x)
    
    kmoire=kbig-ksmall
    ymoire=np.sin(kmoire*x)

    ax.scatter(x,ybig,c='r')
    ax.scatter(x,ysmall,c='g')
    ax.scatter(x,ybig+ysmall,c='b')
    ax.scatter(x,ymoire,c='k')

    plt.show()


def main():
    a=[1,1]
    b=[-2,2]
    real_l=Lattice(a,b)
    real_l=make_hexagonal_lattice()
    real_l=make_square_lattice()

    fig=plt.figure(0)
    ax=fig.add_subplot('111')
    ax.set_aspect("equal")

    arad=150
    brad=150
    point_size=5

    # scatter_lattice_points(ax,real_l,arad,brad,c='red',s=point_size)
    plot_periodic_lattice_cells(ax,real_l,arad,brad,c='red')
    plot_lattice_unit_vectors(ax,real_l,width=0.1,color='red')

    rot_m=make_rotation_matrix(180//3)
    rot_m=make_rotation_matrix(2)
    # rot_m=np.eye(2)*0.8
    rot_l=real_l.transform(rot_m)

    # recip_l=real_l.reciprocal()
    # plot_periodic_lattice_cells(ax,recip_l,arad,brad,c='green')
    # plot_lattice_unit_vectors(ax,recip_l,width=0.1,color='green')

    # plt.show()
    # return

    # scatter_lattice_points(ax,rot_l,arad,brad,c='green',s=point_size)
    plot_periodic_lattice_cells(ax,rot_l,arad,brad,c='green')
    plot_lattice_unit_vectors(ax,rot_l,width=0.1,color='green')

    ax.set_aspect("equal")

    R=real_l.column_vector_matrix()
    M=rot_m
    R_prime=M.dot(R)

    print(R_prime-rot_l.column_vector_matrix())

    G=real_l.reciprocal().column_vector_matrix()
    G_prime=rot_l.reciprocal().column_vector_matrix()
    G_moire=G_prime-G
    K=np.linalg.inv(G).dot(G_prime)

    print(M.dot(K.T))

    print("M\n",M)
    print("K\n",K)
    print("M.T^-1\n",np.linalg.inv(M.T))

    I=np.eye(2)
    P=np.linalg.inv((K-I).T)
    print("P\n",P)
    print(P-np.linalg.inv((K-I).T))


    print("R\n",R)
    print("R'\n",rot_l.column_vector_matrix())
    print("M\n",rot_m)
    print("P\n",P)

    moire_l=Lattice(G_moire[:,0],G_moire[:,1]).reciprocal()
    R_moire=moire_l.column_vector_matrix()
    T=np.linalg.inv(R).dot(R_moire)
    T_prime=np.linalg.inv(R_prime).dot(R_moire);

    print("T\n",T)
    print("T_prime\n",T_prime)

    T=np.rint(T)
    T_prime=np.rint(T_prime)

    print("T\n",T)
    print("T_prime\n",T_prime)
    R_moire=R.dot(T)
    moire_l=Lattice(R_moire[:,0],R_moire[:,1])

    arad=3
    brad=3

    plot_periodic_lattice_cells(ax,moire_l,arad,brad,lw=3,c='k')
    plot_lattice_unit_vectors(ax,moire_l,width=0.5,color='black')

    #############################

    # fig=plt.figure(1)
    # ax=fig.add_subplot('111')
    # ax.set_aspect("equal")

    RR=R.dot(T)
    RR_prime=R_prime.dot(T_prime)


    Z=(RR+RR_prime)//2
    Z_l=Lattice(Z[:,0],Z[:,1])
    plot_lattice_unit_vectors(ax,Z_l,width=0.5,color="gray")
    plot_periodic_lattice_cells(ax,Z_l,arad,brad,lw=3,c='gray')

    # RR_l=Lattice(RR[:,0],RR[:,1])
    # plot_lattice_unit_vectors(ax,RR_l,width=0.5,color="red")
    # plot_periodic_lattice_cells(ax,RR_l,arad,brad,lw=3,c='k')

    # RR_prime_l=Lattice(RR_prime[:,0],RR_prime[0:,1])
    # plot_lattice_unit_vectors(ax,RR_prime_l,width=0.5,color='green')
    # plot_periodic_lattice_cells(ax,RR_prime_l,arad,brad,lw=3,c='green')

    fig=plt.figure(1)
    ax=fig.add_subplot('111')
    ax.set_aspect("equal")

    arad=150
    brad=150

    T=Z.dot(np.linalg.inv(R))
    T_prime=Z.dot(np.linalg.inv(R_prime))

    T=np.rint(T)
    T_prime=np.rint(T_prime)

    R=np.linalg.inv(R)
    R_prime=np.linalg.inv(T).dot(Z)

    R_l=Lattice(R[:,0],R[:,1])
    R_prime_l=Lattice(R[:,0],R[:,1])

    plot_periodic_lattice_cells(ax,R_l,arad,brad,c='green')
    plot_lattice_unit_vectors(ax,R_l,width=0.1,color='green')
    plot_periodic_lattice_cells(ax,R_prime_l,arad,brad,c='red')
    plot_lattice_unit_vectors(ax,R_prime_l,width=0.1,color='red')
    
    plt.show()
    

if __name__=="__main__":
    main()
