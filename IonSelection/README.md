# Configuration Selection

## Quasi-Chemical Theory

QCT treats the ion-waters as a molecular species of the system under analysis, then provides a concise format,

![equation](https://latex.codecogs.com/gif.latex?\mu_{X^-}^{(ex)}&space;=&space;-&space;RT&space;\ln&space;K_n^{(0)}&space;&plus;&space;\rho_{\mathrm{H}_2\mathrm{O}}^{~~~~n}&space;&plus;&space;RT&space;\ln&space;p_{X^-}(n)&space;&plus;&space;\left\(\mu_{(\mathrm{H}_2\mathrm{O})_n}^{(ex)}&space;-&space;n&space;\mu_{\mathrm{H}_2\mathrm{O}}^{(ex)}&space;\right&space;\))

for free energies of solution components X. The populations of clusters are established by applying a clustering 
algorithm, according to which proximal ligands of a specific ion are defined as *inner-shell* partners of that ion. 

Application of QCT to X(aq) thus begins with identification of inner shell configurations of the medium relative to the 
ion. Emphasizing the observation above that the challenge in treating anions lies in the variety of H-bond donations 
structures of ![equation](https://latex.codecogs.com/gif.latex?(\mathrm{H}_2\mathrm{O})_n\mathrm{F}^-), a natural 
procedure is to identify water molecules with H atoms within a distance of λ from X as clustered: That is, as 
inner-shell partners with the distinguished X ion. From there, with n water ligands in the cluster and 
directly interacting with the ion.

Determining the "cut-off" distance for interacting H-Bonds is done through radial distributions provided from *ab initio*
molecular dynamics. 

## Isolated Cluster Energetics

To evaluate cluster energetics, and to address the discrepancy in harmonic approximation and free energy experiments, 
we select configurations containing n water molecules within λ of the ion of interest.

These configurations go through a molecule deletion scheme to evaluate the free energy,

![equation](https://latex.codecogs.com/gif.latex?\Delta&space;U&space;=&space;E(\gamma_n\sigma)&space;-&space;E(\gamma_{n-1}\sigma)&space;-&space;E(\gamma&space;\sigma)&space;&plus;&space;E(\sigma))

so that,

![equation](https://latex.codecogs.com/gif.latex?n&space;K_n^{(0)}&space;=&space;\frac{K_1^{(0)}K_{n-1}^{(0)}}{\langle&space;e^{\beta\Delta&space;U}\rangle}_n)

Where γ denote water ligands and σ is the ion of interest.

## Running

```
python3 main.py -f [cp2k .xyz trajectory file] -i [ion of interest] -r [λ distance in angstroms] -A [Age trajectory] 
```

This will create directories for all σγ_n, σγ_n-1, σγ.

For example, if your cluster is F and three waters (n = 3)

The program will create a director called full, for σγ_n, and 3 directories for σγ_n-1, {AB, AC, BC} consisting of 
configurations with 2 waters and solute, and 3 directories for σγ, {C,B,A}, for the single ligand interactions with the 
solute.

