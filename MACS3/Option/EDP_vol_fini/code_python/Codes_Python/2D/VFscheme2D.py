# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
from scipy.sparse.linalg import spsolve
import mesh2D as m
import utils as u
import data as d
#
# Initialisation des pressions et des saturations d'huile 
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
#
# Données de sortie :
# p     : les pressions initiales [Pa]
# s     : les saturations d'huile initiales [-]
#
def initUnknowns(mesh,data):
    p=np.zeros(mesh.nb_cells)
    s=np.zeros(mesh.nb_cells)
    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]
    for c in range(mesh.nb_cells):
        x=mesh.cells[c][m._X]
        y=mesh.cells[c][m._Y]
        p[c] = data.p_init(x,y,Lx,Ly)
        s[c] = data.s_init(x,y,Lx,Ly)
    return p,s
#
# Initialisation des transmissivités
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
#
# Données de sortie : Tkl : les transmissivités sur les arêtes
# internes et sur les arêtes de bord
#    
def computeTransmissivity(mesh,data):

    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]

    # initialisation des permeabilites dans les mailles
    perm = np.zeros(mesh.nb_cells)
    for c in range(mesh.nb_cells):
        perm[c] = data.perm(mesh.cells[c][m._X],mesh.cells[c][m._Y],Lx,Ly)
        
    # initialisation des transmissivites sur les aretes
    Tkl = np.zeros(mesh.nb_edges)
    for e in range(mesh.nb_edges):
        dk_sigma = mesh.edges[e][m._DKsigma]
        dl_sigma = mesh.edges[e][m._DLsigma]
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        mesure = mesh.edges[e][m._MES]
        if(l>=0):
            # moyenne harmonique des permeabilites
            Tkl[e] = 0.0 # A COMPLETER
        else:
            # arete de bord
            Tkl[e] = 0.0 # A COMPLETER

    return Tkl
#
# Calcul des mobilités totales
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
# p     : les pressions [Pa]
# s     : les saturations d'huile [-]
#
# Données de sortie :
# mobT  : les mobilités totales sur les arêtes internes et 
#         sur les arêtes de bord
#    

def computeTotalMobility(mesh,data,p,s):
    # initialisation des mobilites totales sur les aretes
    mobT = np.zeros(mesh.nb_edges)
    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]
    for e in range(mesh.nb_edges):
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        if(l>=0):
            # arete interne : decentrage amont
            if(p[k]>=p[l]):
                mobT[e] = u.mw(1.0-s[k],data.swi,data.muw)+u.mo(s[k],data.sor,data.muo)
            else:
                mobT[e] = u.mw(1.0-s[l],data.swi,data.muw)+u.mo(s[l],data.sor,data.muo)
        else:
            # arete de bord a pression imopsee:
            if(mesh.edges[e][m._LABEL]==-3):
                x = mesh.edges[e][m._X]
                y = mesh.edges[e][m._Y]
                pbound = data.dirichlet(x,y,Lx,Ly)
                if(p[k]>=pbound):
                    # bord sortant
                    mobT[e] = u.mw(1.0-s[k],data.swi,data.muw)+u.mo(s[k],data.sor,data.muo)
                else:
                    # bord entrant: il ne rentre que de l'eau mobT = mw
                    mobT[e] = u.mw(1.0,data.swi,data.muw)
        
    return mobT
#
# Calcul des pressions par un schéma volumes finis
#
# Données d'entrée 
# mesh : le maillage 
# data : les données du problème
# Tkl  : les transmissivités sur les arêtes 
# mobT : les mobilités totales sur les arêtes
#
# Données de sortie :
# p     : les nouvelles pressions
#    
def compute_pressure(mesh,data,Tkl,mobT):
    #
    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]
    nb_cells = mesh.nb_cells
    A=sp.sparse.dok_matrix((nb_cells,nb_cells))
    b = np.zeros(nb_cells)
        
    for e in range(mesh.nb_edges): 
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        # arete interne
        if(l>=0):
            coeff = 0.0 # A COMPLETER
            A[k,k] = A[k,k] + coeff
            A[l,l] = A[l,l] + coeff
            A[k,l] = A[k,l] - coeff
            A[l,k] = A[l,k] - coeff
        else:
            #aretes de bord
            if(mesh.edges[e][m._LABEL]==-3):
                x = mesh.edges[e][m._X]
                y = mesh.edges[e][m._Y]
                pbound = data.dirichlet(x,y,Lx,Ly)
                coeff = 0.0 # A COMPLETER
                A[k,k] = A[k,k] + coeff
                b[k]   = b[k] + coeff * pbound
            elif(mesh.edges[e][m._LABEL]==-2):
                x = mesh.edges[e][m._X]
                y = mesh.edges[e][m._Y]
                flux = data.neumann(x,y,Lx,Ly)
                b[k] =  b[k] + flux

    # resolution du systeme lineaire
    #print A.todense()
    p = spsolve(A.tocsr(),b)
        
    return p
#
# Calcul des flux
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
# p     : les pressions [Pa]
# Tkl   : les transmissivités sur les arêtes 
# mobT  : les mobilités totales sur les arêtes
#
# Données de sortie :
# fluxT  : les flux sur les arêtes internes et 
#          sur les arêtes de bord
#    
def computeTotalFlux(mesh,data,p,Tkl,mobT):
    # initialisation des flux totaux sur les aretes
    fluxT = np.zeros(mesh.nb_edges)
    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]
    for e in range(mesh.nb_edges):
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        if(l>=0):
            # arete interne : fluxT = Tkl*mobT*(pk-pl)
            fluxT[e] = Tkl[e]*mobT[e]*(p[k]-p[l])
        else:
            # arete de bord a pression imopsee: fluxT = Tkl*mobT*(pk-pbound)
            if(mesh.edges[e][m._LABEL]==-3):
                x = mesh.edges[e][m._X]
                y = mesh.edges[e][m._Y]
                pbound = data.dirichlet(x,y,Lx,Ly)
                fluxT[e] = Tkl[e]*mobT[e]*(p[k]-pbound)

    return fluxT
#
# Calcul des saturations par un schéma volumes finis explicite
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
# sn    : saturations au temps t(n)
# fluxT : les flux sur les arêtes internes et 
#         sur les arêtes de bord
# dt    : valeur du pas de temps
#
# Données de sortie :
# s     : les nouvelles saturations au temps t(n+1)
#    
def computeExplicitSaturation(mesh,data,sn,fluxT,dt):
    # initialisation des nouvelles saturations
    s = sn
    for e in range(mesh.nb_edges):
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        if(l>=0):
            # arete interne :
            if(fluxT[e] >= 0.0):
                coeff = 0.0 # A COMPLETER
            else:
                coeff = 0.0 # A COMPLETER
            # mise a jour des saturations dans les mailles k et l
            s[k] = s[k] - coeff / mesh.cells[k][m._VOL]
            s[l] = s[l] + coeff / mesh.cells[l][m._VOL]
        else:
            #aretes de bord
            if(mesh.edges[e][m._LABEL]!=-1):
                if(fluxT[e] >= 0.0):
                    coeff = 0.0 # A COMPLETER
                else:
                    coeff = 0.0 # A COMPLETER
                # mise a jour de la saturation dans la maille de bord
                s[k] = s[k] - coeff /  mesh.cells[k][m._VOL]
           
    return s
#
# Calcul de la condition CFL
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
# fluxT : les flux sur les arêtes
#
# Données de sortie :
# dt    : pas de temps dt=t(n+1)-t(n)
#
def computeCFL(mesh,data,fluxT):
    # calcul du max de la derivee du flux fractionnaire Fo
    # ce calcul est independant du temps
    if not hasattr(computeCFL, 'max_fop'):
        max_fop = u.max_dfodso(data.swi,data.sor,data.muw,data.muo)
        computeCFL.max_fop = max_fop
    max_fop = computeCFL.max_fop
    # initialisation de la somme des flux rentrant dans les cellules
    fluxT_in = np.zeros(mesh.nb_cells)
    for e in range(mesh.nb_edges):
        k = int(mesh.edges[e][m._K])
        l = int(mesh.edges[e][m._L])
        if(l>=0):
            # arete interne
            if(fluxT[e]<0.0):
                fluxT_in[k] = fluxT_in[k] - fluxT[e]
            else:
                fluxT_in[l] = fluxT_in[l] + fluxT[e]
        else:
            if(mesh.edges[e][m._LABEL]==-3):
                # arete de bord a pression imopsee: fluxT = Tkl*mobT*(pk-pbound)
                if(fluxT[e]<0.0):
                    fluxT_in[k] = fluxT_in[k] - fluxT[e]
    #
    max_flux_in = max(fluxT_in)
    min_vol = 123456789.0
    for c in range(mesh.nb_cells):
        min_vol = min(min_vol,mesh.cells[c][m._VOL])
    # calcul du pas de temps
    dt =1.0* (data.phi*min_vol)/(max_flux_in*max_fop)
    
    return dt
    