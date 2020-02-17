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

def make_hcp_lattice():
    a=np.array([0,3])
    R=make_rotation_matrix(60)
    b=R.dot(a)
    return Lattice(a,b)

def main():
    a=[1,1]
    b=[-2,2]
    real_l=Lattice(a,b)
    real_l=make_hcp_lattice()

    fig=plt.figure()
    ax=fig.add_subplot('111')
    ax.set_aspect("equal")

    arad=140
    brad=140
    oint_size=5

    # scatter_lattice_points(ax,real_l,arad,brad,c='red',s=point_size)
    plot_periodic_lattice_cells(ax,real_l,arad,brad,c='red')
    plot_lattice_unit_vectors(ax,real_l,width=0.1,color='red')

    rot_m=make_rotation_matrix(180//3)
    rot_m=make_rotation_matrix(1)
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

    print("T\n",T)
    print("T\n",np.round(T))
    R_moire=R.dot(T)
    moire_l=Lattice(R_moire[:,0],R_moire[:,1])

    plot_periodic_lattice_cells(ax,moire_l,3,3,c='k')
    plot_lattice_unit_vectors(ax,moire_l,width=0.1,color='black')


    plt.show()
    

if __name__=="__main__":
    main()
