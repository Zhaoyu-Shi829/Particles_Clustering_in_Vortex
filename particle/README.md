## Lagrangian Point Particle Tracking

* **inertial particle parameter space**: <br>
  * The appropriate choice of drag model for individual particle's movement depends on particle properties, such as particle diameter ($d/\eta$), density ratio ($\rho_p/\rho_f$) and volume fraction ($\Phi=N_pV_p/V$).
  * The non-dimensional parameter `Stokes number` is a measure of partical inertia, e.g. $Sk=\tau_p/\tau_f(\tau_p=\rho_pd^2/18\rho_f\nu)$. $\tau_f$ is the characteristic fluid time scale, and it can be nominal (cylinder diameter $D$ and incoming flow velocity $U_0$) or physics-based scale (local flow field $1/|\omega_z|$, thus ***effective Sk***).
    <p align="center">
      <img src="https://github.com/user-attachments/assets/f81576e9-041f-45b5-83ea-79409c03d625" width="400">
    </p>
  * The ***point*** particle approach requires particle size smaller than the smallest scale of flow ($d/\eta<0.1$, kolmogorov scale $\eta$) 
* **Reduced-order Maxey-Riley equation**  <br>
  * We consider heavy particles ($\rho_p/\rho_f>10^3$), thus viscous drag force is the only force with the non-linear Faxen correction for the finite-size effect.
      <p align="center">
        <img src="https://github.com/user-attachments/assets/cffe0306-7384-46ee-bb66-1c607764871e" width="200">
      </p>
  * The drag coefficient is a function of particle Reynolds number $`C_D=f(Re_p)`$, and depends on particle surface roughness, fluid viscosity, particle shape and its rotation, flow compressibility etc. Herein only smooth sphere solid particle considered. 
* **Identification and quantification of Preferential concentration/Cluster**
      <p align="center">
        <img src="https://github.com/user-attachments/assets/c257ba31-292e-4635-ab42-72ab465911c6" width="790">
      </p>
  * Vironoi diagram: measure density inhomogeneity of particle field, based on relative-position of particle pairs; an essemble of polyheron cells.
      <p align="center">
        <img src="https://github.com/user-attachments/assets/004751b0-e069-48f7-992d-2abf6e643d87" width="800">
      </p>
  * Conventionally, the scale of **clusters** and **voids** are determined from PDF of Voronoi cell size/volume, but this does NOT work for vortex-dominated flows.
  * **`A new identification is proposed in this thesis, that is based on the correlation of flow field quantities`** (vorticity $|\omega_{x/z}|$, $\lambda_2$, $Q$-criterion).
  















