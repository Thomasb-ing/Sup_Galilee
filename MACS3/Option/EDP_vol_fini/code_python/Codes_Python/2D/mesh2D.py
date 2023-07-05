# -*- coding: utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt
import data as d
mesh_dir = "./meshes/"

_X = 0
_Y = 1
_VOL = 2
_V1 = 3
_V2 = 4
_V3 = 5
_V4 = 6

_DEB = 2
_FIN = 3
_K   = 4
_L   = 5
_DKsigma = 6
_DLsigma = 7
_DKL     = 8
_MES     = 9
_LABEL   = 16

_CELL = 0

class mesh2D:
    def __init__(self,name,nb_cells,nb_edges,nb_vertices,cells,edges,vertices,dim):
        self.name        = name
        self.nb_cells    = nb_cells
        self.nb_vertices = nb_vertices
        self.nb_edges    = nb_edges
        self.cells       = cells
        self.edges       = edges
        self.vertices    = vertices
        self.dim         = dim

        
def read_mesh(filename):
    # lecture des coordonnees des sommets
    f=open(mesh_dir+filename+"_sommets","r")    
    for ligne in f:
        if ligne.startswith(" nbs"):
            chaine = ligne.split()
            nbs = int(chaine[1])
            # print nbs
            break
    vertices = np.zeros((nbs,2))
    i=0
    for ligne in f:
        chaine = ligne.split()
        vertices[i,0] = float(chaine[0])
        vertices[i,1] = float(chaine[1])
        i=i+1
        print(chaine[1])
    f.close()
     # lecture des donnees sur les aretes
    f=open(mesh_dir+filename+"_aretes","r")
    for ligne in f:
        if ligne.startswith(" nba"):
            chaine = ligne.split()
            nba = int(chaine[1])
            #print nba
            break
    edges = np.zeros((nba,17))
    i=0
    for ligne in f:
        chaine = ligne.split()
        for j in range(0,17):
            edges[i][j] = float(chaine[j])
        i = i+1
    f.close()
    # lecture des centres des mailles
    f=open(mesh_dir+filename+"_centres","r")    
    for ligne in f:
        if ligne.startswith(" nbt"):
            chaine = ligne.split()
            nbc = int(chaine[1])
            # print nbc
            break
    cells = np.zeros((nbc,2))
    i=0
    for ligne in f:
        chaine = ligne.split()
        cells[i][0] = float(chaine[0])
        cells[i][1] = float(chaine[1])
        i = i+1
    f.close()
    # initialisation de la structure maillage
    return mesh2D(filename,nbc,nba,nbs,cells,edges,vertices,[1,1])
#
# Construction d'un maillage cartésien régulier 2D
#
# Données d'entrée
# Lx    : longueur du domaine dans la direction x
# Ly    : longueur du domaine dans la direction y
# Nx    : nombre d'intervalles dans la direction x
# Ny    : nombre d'intervalles dans la direction y
#
# Données de sortie :
# mesh2D  : la structure du maillage 2D
#    
def cartesian2D_mesh(Lx,Ly,Nx,Ny):
    #
    name = 'cartesian_grid'+str(Nx)+'_'+str(Ny)
    dx = Lx/float(Nx)
    dy = Ly/float(Ny)
    #print('dx=',dx)
    #print('dy=',dy)
    # initialisation du nombre de sommets, d'aretes et de cellules
    nbs = (Nx+1)*(Ny+1)
    nbe = (Nx+1)*Ny+Nx*(Ny+1)
    nbc = (Nx*Ny)
    
    x = np.zeros(Nx+1)
    y = np.zeros(Ny+1)
    for i in range(Nx+1):
        x[i] = i*dx
    for j in range(Ny+1):
        y[j] = j*dy

    vertices = np.zeros((nbs,2))
    for j in range(Ny+1):  
        for i in range(Nx+1):
            v = i + j*(Nx+1)
            vertices[v][_X] = x[i]
            vertices[v][_Y] = y[j]
            #print('Vertice(',v,')= ',vertices[v])
            
    edges = np.zeros((nbe,17))
    e=0
    # initialisation des aretes verticales
    mesure = dy
    for i in range(Nx+1):
        for j in range(Ny):
            k = i-1+j*Nx
            l = i+j*Nx
            s = i+j*(Nx+1)
            dk_sigma = dx/2.0
            dl_sigma = dk_sigma
            # bord gauche
            if(i==0):
                k=l
                l=-1
                dl_sigma = 0.0
            # bord droit
            if(i==Nx):
                l=-1
                dl_sigma = 0.0
                
            edges[e][_X] = x[i]
            edges[e][_Y] = (y[j]+y[j+1])/2.0
            edges[e][_DEB] = s
            edges[e][_FIN] = s + Nx + 1
            edges[e][_K]   = k
            edges[e][_L]   = l
            edges[e][_DKsigma] = dk_sigma
            edges[e][_DLsigma] = dl_sigma
            edges[e][_MES]     = mesure
            edges[e][_LABEL]   = -1
            #print("edge(",e,")= ",edges[e])
            e = e+1
     # initialisation des aretes horizontales
    mesure = dx
    for i in range(Nx):
        for j in range(Ny+1):
            k = i+(j-1)*Nx
            l = i+j*Nx
            s = i+j*(Nx+1)
            dk_sigma = dy/2.0
            dl_sigma = dk_sigma
            # bord bas
            if(j==0):
                k=l
                l=-1
                dl_sigma = 0.0
            # bord haut
            if(j==Ny):
                l=-1
                dl_sigma = 0.0
                
            edges[e][_X] = (x[i]+x[i+1])/2.0
            edges[e][_Y] = y[j]
            edges[e][_DEB] = s
            edges[e][_FIN] = s + 1
            edges[e][_K]   = k
            edges[e][_L]   = l
            edges[e][_DKsigma] = dk_sigma
            edges[e][_DLsigma] = dl_sigma
            edges[e][_MES]     = mesure
            edges[e][_LABEL]   = -1
            #print("edge(",e,")= ",edges[e])
            e = e + 1
    # initialisation des cellules
    cells = np.zeros((nbc,7))
    volume = dx*dy
    for i in range(Nx):
        for j in range(Ny):
            c = i+j*Nx
            cells[c][_X] = (x[i]+x[i+1])/2.0
            cells[c][_Y] = (y[j]+y[j+1])/2.0
            cells[c][_VOL] = volume
            cells[c][_V1] = i+j*(Nx+1)
            cells[c][_V2] = i+j*(Nx+1) + 1
            cells[c][_V3] = i+(j+1)*(Nx+1)
            cells[c][_V4] = i+(j+1)*(Nx+1) + 1
            #print("cell(",c,")= ", cells[c])
                                          
    return mesh2D(name,nbc,nbe,nbs,cells,edges,vertices,[Lx,Ly])
            
def plot_mesh(m):
    for e in m.edges:
        #i1 = int(e[_DEB])-1
        #i2 = int(e[_FIN])-1
        i1 = int(e[_DEB])
        i2 = int(e[_FIN])
        lx = [m.vertices[i1][_X],m.vertices[i2][_X]]
        ly = [m.vertices[i1][_Y],m.vertices[i2][_Y]]
        plt.plot(lx,ly,'b')
    plt.show()

def plot_val(m,val):
    Lx = m.dim[_X]
    Ly = m.dim[_Y]
    xlist=np.linspace(0.0,Lx,100)
    ylist=np.linspace(0.0,Ly,100)
    X,Y = np.meshgrid(xlist,ylist)

    plt.figure()
    cp=plt.contourf(X,Y,val)
    plt.colorbar(cp)
    plt.title("pressure")
    plt.show()
#
# Affecte le type de condition aux limites sur les arêtes de bord
# du maillage:
#
#        Conditions de Neumann homogenes  = -1
#        Conditions de Neumann non homogenes  = -2
#        Conditions de Dirichlet = -3
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
#
#
def setLabelBoundaryConditions(mesh,data):
    for e in range(mesh.nb_edges):
        l = mesh.edges[e][_L]
        if(l<0):
            # Label par defaut : condition de flux nul
            mesh.edges[e][_LABEL] = -1
            x = mesh.edges[e][_X]
            y = mesh.edges[e][_Y]
            xmax = mesh.dim[_X]
            ymax = mesh.dim[_Y]
            valD=data.dirichlet(x,y,xmax,ymax)
            valN=data.neuman(x,y,xmax,ymax)
            if(valD != d.NOVAL):
                mesh.edges[e][_LABEL] = -3
            elif(valN != d.NOVAL):
                mesh.edges[e][_LABEL] = -2
            

                         
                         
             
