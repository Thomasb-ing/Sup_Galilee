# -*- coding: utf-8 -*-
import numpy as np
NOVAL = -123456789.0
_X=0
_Y=1
_R=2
_P=3
#
class data:
    def __init__(self,nom,perm,poros,swi,sor,muw,muo,Lx,Ly,Nx,Ny,Tend,ndtMax,dirichlet,neuman,p_init,s_init):
        self.name = nom
        self.perm = perm
        self.phi  = poros
        self.swi = swi
        self.sor = sor
        self.muw = muw
        self.muo = muo
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = Nx
        self.Ny = Ny
        self.Tend = Tend
        self.ndt_max = ndtMax
        self.dirichlet = dirichlet
        self.neuman   = neuman
        self.p_init   = p_init
        self.s_init   = s_init

# permeabilite [m2]
def perm_homogene(x,y,xmax,ymax):
    perm = 1.0E-12
    return perm

# porosité du domaine
phi0 = 0.1
# saturation irreductible en eau
swi  = 0.1
# saturation residuelle en huile
sor  = 0.3
# viscosite de l'eau [Pa.s]
muw = 1.0E-3
# viscosite de l'huile [Pa.s]
muo = 5.0E-3
# dimensions du domaine en [m]
Lx1 = 1000
Ly1 = 500
# le pas de discretisation dans chaque direction
Nx1 = 100
Ny1 = 10
# Temps  final (secondes [s])
Tend = 3600*24*365*10 # 10 ans en seconde
# nombre max de pas de temps
ndtMax = 1000
# conditions aux limites de type Dirichlet (pression [Pa])
def dirichlet1(x,y,xmax,ymax):
    if(x==0.0):
        return 110.0E05
    elif(x==xmax):
        return 95.0E05
    else:
        return NOVAL
# conditions aux limites de type Neumann homogene    
def neuman_homogene(x,y,xmax,ymax):
    return NOVAL

# conditions initiales en pression
def p_init(x,y,xmax,ymax):
    p=95.0E05
    return p

# conditions initiales en saturation d'huile
def s_init(x,y,xmax,ymax):
    s=1.0-swi
    return s

# Construction de la liste des tests
liste_test=[]
#
# Test 1
#
liste_test.append(data("test1: homogène + CL Dirichlet",perm_homogene,
                       phi0,swi,sor,muw,muo,Lx1,Ly1,Nx1,Ny1,Tend,ndtMax,
                       dirichlet1,neuman_homogene,p_init,s_init))
#
# Test 2
#

def perm_heterogene(x,y,xmax,ymax):
    if(x<=xmax/2.0):
        perm=1.0E-12
    else:
        perm=10.0E-12
    return perm
#
liste_test.append(data("test2: hétérogène + CL Dirichlet",perm_heterogene,
                       phi0,swi,sor,muw,muo,Lx1,Ly1,Nx1,Ny1,Tend,ndtMax,
                       dirichlet1,neuman_homogene,p_init,s_init))

#
# Test 3
#    
# condition aux limites de type Dirichlet
def dirichlet3(x,y,xmax,ymax):
    if(x==0.0 and y<=50.0) or (y==0.0 and x<=10.0):
        return 110.0E05
    elif(x==xmax and y>=ymax-50.0) or (x>=xmax-10.0 and y==ymax):
        return 95.0E05
    else:
        return NOVAL
#    
liste_test.append(data("test3: homogène + CL Dirichlet aux deux coins opposés",perm_homogene,
                        phi0,swi,sor,muw,muo,Lx1,Ly1,Nx1,Ny1,Tend,ndtMax,
                        dirichlet3,neuman_homogene,p_init,s_init))
