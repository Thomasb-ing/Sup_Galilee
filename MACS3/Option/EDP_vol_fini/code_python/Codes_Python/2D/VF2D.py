#!/usr/bin/env python2.7
import numpy as np 
import mesh2D as m
import data as d
import VFscheme2D as vf
import impes_scheme as impes
import utils as u


print ("-----------------------------------")
print (" Liste des cas test:")
index=1
for test in d.liste_test:
    print (str(index)+'-'+test.name)
    index=index+1
    
choix=-1
while choix <=0 or choix > len(d.liste_test):
    choix = int(input("choisissez un test:"))

data = d.liste_test[choix-1]
print ("===== Execution du test "+data.name+" =====")
# construction du maillage
mesh = m.cartesian2D_mesh(data.Lx,data.Ly,data.Nx,data.Ny)
# affectation des conditions aux limites sur le maillage
m.setLabelBoundaryConditions(mesh,data)
# initialisation des inconnues
p,sn = vf.initUnknowns(mesh,data)
s=sn
# initiallisation des transmissivites sur les aretes
Tkl = vf.computeTransmissivity(mesh,data)
# initialisation du temps
t=0.0
dt=0.05 
ndt=0
# archivage de l'etat initial au format VTK 
u.write_result_in_VTK_format(ndt,mesh,data,p,1.0-sn,sn)
 # boucle en temps
while (t < data.Tend and ndt < data.ndt_max):
    # resolution d'un pas de temps par le schema IMPES
    p,s,dt = impes.solve_time_step(mesh,data,Tkl,sn,p,s,t,dt)
    #
    # mise a jour des inconnues pour le pas de temps suivant
    #
    t = t + dt
    ndt = ndt + 1
    sn = s
    #print ("===== Step= "+str(ndt)+", dt= "+str(dt)+", TIME= "+str(t))
    print ("===== Step= {:<5d}, dt= {:<10.5f}, time= {:10.5f}".format(ndt,dt,t))
    # archivage des resultats au format VTK
    if (ndt%10==0):
        print ("Archivage des resultats au temps {:10.5f}".format(t))
        u.write_result_in_VTK_format(ndt,mesh,data,p,1.0-sn,sn)
# archivage des resultats au format VTK au temps final 
print ("Archivage des resultats au temps {:10.5f}".format(t))
u.write_result_in_VTK_format(ndt,mesh,data,p,1.0-sn,sn)

