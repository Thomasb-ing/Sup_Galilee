# -*- coding: utf-8 -*-
import VFscheme2D as vf

#
# Calcul des pressions et des saturations à un pas de temps
# en utilisant le schéma IMPES: 
#
# Données d'entrée
# mesh  : le maillage
# data  : les données du problème
# Tkl   : les transmissivités
# sn    : les saturations d'huile au temps t(n) [-]
# p     : les pressions au temps t(n) [Pa]
# s     : les saturations d'huile au temps t(n+1) [-]
# t     : t(n) [s]
# dt    : le pas de temps t(n+1)-t(n) [s]
#
# Donnees de sortie :
# p     : les pressions au temps t(n+1) [Pa]
# s     : les saturations d'huile au temps t(n+1) [-]
# dt    : le pas de temps t(n+1)-t(n) [s]
#
def solve_time_step(mesh,data,Tkl,sn,p,s,t,dt):
    #
    # calcul des mobilites totales
    #
    mobT = vf.computeTotalMobility(mesh,data,p,sn)
    #
    # calcul des nouvelles pressions
    #
    p = vf.compute_pressure(mesh,data,Tkl,mobT)
    #
    # calcul des flux
    #
    fluxT = vf.computeTotalFlux(mesh,data,p,Tkl,mobT)
    #
    # calcul du pas de temps respectant la condition CFL
    #
    dt = vf.computeCFL(mesh,data,fluxT)
    #
    # correction du pas de temps en cas de depassement du temps final
    #
    if(t+dt > data.Tend):
        dt = data.Tend-t
    #print ("dt= "+str(dt))
    #
    # calcul des nouvelles saturations
    #
    s = vf.computeExplicitSaturation(mesh,data,sn,fluxT,dt)
    return p,s,dt


