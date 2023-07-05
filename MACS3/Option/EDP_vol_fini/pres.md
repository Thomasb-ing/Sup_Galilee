---
marp: true
author: Thomas Binet
title: Test
theme: dracula
paginate: true
# backgroundColor: 
#backgroundImage: url('https://marp.app/assets/hero-background.svg')
header: 'Présentation résultat projet extraction pétrole'
footer: 'T. Binet G. Polet'
math: katex
style: |
  video::-webkit-media-controls {
    will-change: transform;
  }

---

<!-- _header: ''-->
<!-- _footer: ''-->

![bg w:500](https://uploads-ssl.webflow.com/5dea07d5bb83abf6cdffcf8a/5e201b59924eb42cc03976e3_14685.svg)

### <!--fit><--> *Modélisation fluide diphasique incomprésible*

$$
\begin{align*}
\text{par } &BINET &Thomas\\
\text{et } &POLET &Guillaume
\end{align*}
$$





---
<!--_class: fmmmm-->
* Questions de cours
  * Calcul des transmissivités sur les arêtes internes
  * Calcul des mobilités totales aux temps $t^{(n)} \geq 0$.
  * Equation discrète en pression.
  * Equation discrète en saturation.
  * Ecrire et démontrer la condition CFL en tenant compte des puits.
* Simulations numeriques
  * avec viscosité $\mu_w=1\cdot10^{-10} Pa.s$
  * avec viscosité $\mu_w=3\cdot10^{-10} Pa.s$

---

### Présentation du Modèle

---

### Calcul formule Harmonique

$$F_{K\sigma}=-|\sigma|\kappa_K \frac{(P_{\sigma}-P_K)}{d_{K\sigma}}$$
$$F_{L\sigma}=-|\sigma|\kappa_L \frac{(P_{\sigma}-P_L)}{d_{L\sigma}}$$

On suppose que $F_{K\sigma} =-F_{L\sigma}$ ce qui nous permet d'avoir une expression $P_\sigma=f(P_K,P_L)$

---

### On obtient alors

$$
\begin{align*}
F_{K\sigma} =-F_{L\sigma} &\Leftrightarrow  |\sigma|\kappa_K \frac{(P_{\sigma}-P_K)}{d_{K\sigma}}-|\sigma|\kappa_L \frac{(P_{\sigma}-P_L)}{d_{L\sigma}} &=0 \\
 &\Leftrightarrow  \sigma|\kappa_K \frac{(P_{\sigma}-P_K)}{d_{K\sigma}}-|\sigma|\kappa_L \frac{(P_{\sigma}-P_L)}{d_{L\sigma}} &=0 \\
 &\Leftrightarrow  \sigma|\kappa_K \frac{(P_{\sigma}-P_K)}{d_{K\sigma}}-|\sigma|\kappa_L \frac{(P_{\sigma}-P_L)}{d_{L\sigma}} &=0 \\
 &\Leftrightarrow  \sigma|\kappa_K \frac{(P_{\sigma}-P_K)}{d_{K\sigma}}-|\sigma|\kappa_L \frac{(P_{\sigma}-P_L)}{d_{L\sigma}} &=0
\end{align*}
$$

---

### Expression de $F_{K\sigma}$ sans $P_{\sigma}$

On obtient alors une expression du flux entrant dans la cellule K par la surface $\sigma$ valant:

$$\boxed{F_{K\sigma}=|\sigma| \cdot \frac{\kappa_K \kappa_L}{\kappa_L d_{K\sigma}+\kappa_K d_{L\sigma}}}$$

---

### Simulations numeriques

Dans les prochaines sections, on va simuler numeriquement l'évolution du milieu poreux soumis à des contraites sur 6 mailles dites "perforés" (i.e puits INJECTEUR ou PRODUCTEUR sur la maille).

* Simulation avec viscusité $\mu_o=1\cdot10^{-10} Pa.s$
* Simulation avec viscusité $\mu_o=3\cdot10^{-10} Pa.s$

---

### Resultat avec $\mu_o=1\cdot10^{-10} Pa.s$

<video src="animation_test1_speedup40.mp4" controls width=800></video>

---

### Remarques

On remarque que la récuperation de petrole est plus importante quand la viscosité de l'eau se rapproche de celle de l'huile.

---

![center](https://kroki.io/mermaid/svg/eNpNzr0OgjAcBPCdp7hRB-QNNHz4Masb6VBK0zZg_6S2GhTf3cJgnH93l1OODxrXKgHyGrkdyUowpOl2KrmFlv2ACcUKR4InKON1aDaCbtkYXsFLobPOUWewjgvF3EP5xomec1qQ9c40MbbDJ3q5eFXjLAdy3liFJqg72M_2NS6au1lMK_k_HeK99kGCLz2WfAHnUzg5)
