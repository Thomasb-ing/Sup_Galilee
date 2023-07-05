# Installation_Test

Tests de l'installation des différents outils nécessaires au cours IN203.

## Pré-requis

### Utilisateurs sous environnement Windows et MSYS 2

Editer le fichier "Makefile" comme suit :

    # include Make.inc
    include Make_msys2.inc

### Utilisateurs sous environnement MacOS

Procéder à l'installation de homebrew, puis des paquets gcc et open-mpi :

    $ brew install gcc open-mpi

Éditer ensuite le fichier "Makefile" comme suit :

    # include Make.inc
    include Make_osx.inc

## Exécution du test

Tapez :

    make test

Si tout se passe bien, vous devriez avoir des sorties de divers programmes sur la console.

Si les messages vous semblent sortir *dans le désordre ou bien sont mélangés*, ne vous inquiétez pas,
**c'est un comportement tout à fait normal lorsqu'on lance des programmes en parallèle**. 

