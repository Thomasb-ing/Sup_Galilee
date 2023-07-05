# -*- coding: utf-8 -*-
import numpy as np
import mesh2D as m
#
# Calcul de la mobilite de la phase huile : 
#
#           Mo = kro(So)/muo
#
# ou Kro(So) est la permeabilite relative de l'huile : 
#
#           kro(So) = (so - sor)^2 / (1-sor)^2
#
# Donnees d'entree
# So    : saturation de la phase huile [-]
# Sor   : saturation residuelle de la phase huile [-]
# muo   : viscosite de la phase huile [Pa.s]
#
# Donnees de sortie :
# moboil    : mobilite de l'huile  [Pa^-1.s^-1]
#
def mo(so,sor,muo):
    if so<sor :
        return(0.0)
    else:
        return((so-sor)**2.0/(1.0-sor)**2/muo)
#
# Calcul de la mobilite de la phase eau : 
#
#           Mw = krw(Sw)/muw
#
# ou Krw(Sw) est la permeabilite relative de l'eau: 
#
#           krw(Sw) = (sw - swi)^2 / (1-swi)^2
#
# Donnees d'entree
# Sw    : saturation de la phase eau [-]
# Swi   : saturation irreductible de la phase eau [-]
# muw   : viscosite de la phase eau [Pa.s]
#
# Donnees de sortie :
# mobwater  : mobilite de l'eau  [Pa^-1.s^-1]
#
def mw(sw,swi,muw):
    if sw<swi : 
        return(0.0)
    else:
        return((sw-swi)**2.0/(1.0-swi)**2/muw)
#
# Calcul de la derivee de la mobilite de la phase huile : 
#
#           dModSo = dkrodSo(So)/muo
#
# ou dKrodSo(So) est la derivee de la permeabilite relative de l'huile : 
#
#           dkrodSo(So) = 2*(so - sor) / (1-sor)^2
#
# Donnees d'entree
# So    : saturation de la phase huile [-]
# Sor   : saturation residuelle de la phase huile [-]
# muo   : viscosite de la phase huile [Pa.s]
#
# Donnees de sortie :
# dModSo : derivee de la mobilite de l'huile  [Pa^-1.s^-1]
#
#
def dmodso(so,sor,muo):
    if so<sor : 
        return(0.0)
    else:
        return(2.0*(so-sor)/(1.0-sor)**2.0/muo)
#
# Calcul de la derivee de la mobilite de la phase eau : 
#
#           dMwdSw = dkrwdSw(Sw)/muw
#
# ou dKrwdSw(Sw) est la derivee de la permeabilite relative de l'eau : 
#
#           dkrwdSw(Sw) = 2*(sw - swi) / (1-swi)^2
#
# Donnees d'entree
# Sw    : saturation de la phase eau [-]
# Swi   : saturation irreductible de la phase eau [-]
# muw   : viscosite de la phase eau[Pa.s]
#
# Donnees de sortie :
# dMwdSw : derivee de la mobilite de l'eau  [Pa^-1.s^-1]
#
#
def dmwdsw(sw,swi,muw):
    if sw<swi:
        return(0.0)
    else:
        return(2.0*(sw-swi)/(1.0-swi)**2.0/muw)
#
# Calcul du flux fractionnaire : 
#
#           Fo(so)= Mo(So)/(Mo(So)+Mw(1-So))
#
# Donnees d'entree
# So    : saturation de la phase huile [-]
# Swi   : saturation irreductible de la phase eau [-]
# Sor   : saturation residuelle de la phase huile [-]
# muw   : viscosite de la phase eau [Pa.s]
# muo   : viscosite de la phase huile [Pa.s]
#
# Donnees de sortie :
# fo : flux fractionnaire [-]
#
#
def fo(so,swi,sor,muw,muo):
    sw=1.0-so
    mobo = mo(so,sor,muo)
    mobw = mw(sw,swi,muw)
    return(mobo/(mobo+mobw))
#
# Calcul de la derivee du flux fractionnaire : 
#
#           dFodSo(so)= Mo(So)/(Mo(So)+Mw(1-So))
#
# Donnees d'entree
# So    : saturation de la phase huile [-]
# Swi   : saturation irreductible de la phase eau [-]
# Sor   : saturation residuelle de la phase huile [-]
# muw   : viscosite de la phase eau [Pa.s]
# muo   : viscosite de la phase huile [Pa.s]
#
# Donnees de sortie :
# dfodso : derivee du flux fractionnaire [-]
#
#
def dfodso(so,swi,sor,muw,muo):
    sw=1.0-so
    mobo = mo(so,sor,muo)
    mobw = mw(sw,swi,muw)
    dmobodso = dmodso(so,sor,muo)
    dmobwdsw = dmwdsw(sw,swi,muw)
    return((dmobodso*mobw+mobo*dmobwdsw)/(mobo+mobw)**2)

#
# Calcul de la valeur max de la derivee du flux fractionnaire : 
#
#           max_dFodSo(so)= Max(dfodso), 0<=so<=1
#
# Donnees d'entree
# Swi   : saturation irreductible de la phase eau [-]
# Sor   : saturation residuelle de la phase huile [-]
# muw   : viscosite de la phase eau [Pa.s]
# muo   : viscosite de la phase huile [Pa.s]
#
# Donnees de sortie :
# max_dfodso : max de la derivee du flux fractionnaire [-]
#
#
def max_dfodso(swi,sor,muw,muo):
    max_fop=-1.0
    for i in range(1000):
        so = i/999.0
        fop = dfodso(so,swi,sor,muw,muo)
        if(fop>max_fop):
            max_fop = fop
            
    return(max_fop)
#
# Ecriture des resultats dans un fichier au format VTK
# pour une visualisation avec Paraview
#
# Donnees d'entree
# n     : le numéro du pas de temps
# mesh  : le maillage
# p     : les pressions[Pa]
# sw    : les saturations d'eau [-]
# so    : les saturations d'huile [-]
#
def write_result_in_VTK_format(n,mesh,data,p,sw,so):
    # ouverture du fichier
    file = open("output/res_grid"+str(n)+".vtk","w")
    file.write("# vtk DataFile Version 3.98\n")
    file.write("Exemple RECTILINEAR_GRID\n")
    file.write("ASCII\n")
    file.write("DATASET UNSTRUCTURED_GRID\n")
    file.write("POINTS "+str(mesh.nb_vertices)+" float \n")
    for i in range(mesh.nb_vertices):
        file.write(str(mesh.vertices[i][m._X])+" ")
        file.write(str(mesh.vertices[i][m._Y]))
        file.write(" 0.0\n")
    file.write("CELLS "+str(mesh.nb_cells)+" "+str(mesh.nb_cells*5)+"\n")
    for i in range(mesh.nb_cells):
        file.write("4 "+ str(int(mesh.cells[i][m._V1])) + " ")
        file.write(str(int(mesh.cells[i][m._V2])) + " ")
        file.write(str(int(mesh.cells[i][m._V3])) + " ")
        file.write(str(int(mesh.cells[i][m._V4])) + "\n")
    file.write("CELL_TYPES "+str(mesh.nb_cells)+"\n")
    for i in range(mesh.nb_cells):
        file.write("8\n")
    file.write("CELL_DATA "+str(mesh.nb_cells)+"\n")
    # Valeurs des permeabilités
    file.write("SCALARS permeability float\n")
    file.write("LOOKUP_TABLE default\n")
    Lx = mesh.dim[m._X]
    Ly = mesh.dim[m._Y]
    for i in range(mesh.nb_cells):
        perm = data.perm(mesh.cells[i][m._X],mesh.cells[i][m._Y],Lx,Ly)
        file.write(str(perm)+"\n")
    # Valeurs des pressions
    file.write("SCALARS pressure float\n")
    file.write("LOOKUP_TABLE default\n")
    for i in range(mesh.nb_cells):
        file.write(str(p[i])+"\n")
    # Valeurs des saturations d'eau
    file.write("SCALARS water_saturation float\n")
    file.write("LOOKUP_TABLE default\n")
    for i in range(mesh.nb_cells):
        file.write(str(sw[i])+"\n")
    # Valeurs des saturations d'huile
    file.write("SCALARS oil_saturation float\n")
    file.write("LOOKUP_TABLE default\n")
    for i in range(mesh.nb_cells):
        file.write(str(so[i])+"\n")

    file.close()
                   
        

