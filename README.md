## Numerical Studies of Particle Clustering in Wake Flows 
* 🚦 This project documents [`my doctoral thesis`](https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/2976848). The emphasis of my PhD project is to systmatically investigate particle concentration and dispersion in vortex-dominated flows by leveraging DNS and HPC. This work spans a range of flow regimes, from <ins>laminar</ins> to <ins>turbulent</ins>, and particle properties such as <ins>density</ins>, <ins>size</ins> and <ins>volume fraction</ins>. To explain the intriguing **`clustering patterns`** of particles, I provided an in-depth study of the interaction between the dispersed particle field and coherent votex structures.
* 💡The principal goal is, from a fundamental perspective, to provide a better understanding of particle-flow interaction in many particle-laden processes, such as liquid-fueled combustion efficiency, particle dynamics in drug delivery, powder particle distribution in 3D printing, particle suspension and deposition in processing bed/reactor, source and sink of particulate pollution and many more manufacturing processes. 
* 💻 The DNS/LES accoustic solver [`MGLET`](https://km-turbulenz.de/products-services/) is used to solve particle-fluid flow dynamics
    * 🕸️ Multi-Lever grid mesh generation: `local refinement` in a hierarchical arrangement; <br>
    * :accessibility:  Cut-cell immersed boundary condition:
            <p align="center">
                 <img src="https://github.com/user-attachments/assets/85efa209-f607-4ba0-b07f-233c6ed91704" width="400"> <br>
                 <b>Multi-level grid box consisting of $`N^3`$ uniform cubic cells</b>
            </p>
    * 💧 Particles' movement is tracked by solving reduced-order Maxey-Riley equations;
            <p align="center">
                <img src="https://github.com/user-attachments/assets/6ad9db1a-8d29-40a6-8bb4-7db282992e0f" width="500"> <br>
                <b>Particle distribution in the Transition-in-Wake</b>
            </p>
* My doctoral work includes several very interesting findings, such as the effect of particle-wall collision model resulting in the bow-shock pattern quantified by Voronoı̈ diagram, particle dynamics influenced by local vortices, scale-dependency of clutering in the transitional and turbulence wake, and correlation-based criterion for concentration identification. Check all my publications:
    * [`Bow shock clustering in particle-laden wetted cylinder flow`](https://www.sciencedirect.com/science/article/pii/S0301932220300562)
    * [`Clusters and coherent voids in particle-laden wake flow`](https://www.sciencedirect.com/science/article/pii/S0301932221001269)
    * [`Scale-dependent particle clustering in transitional wake flow`](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/scaledependent-particle-clustering-in-transitional-wake-flow/59A4E829E6ED4E013B208370B7D607C7)
    * [`Different topologies of natural vortex dislocations in Mode A wake`](https://pubs.aip.org/aip/pof/article/34/2/021702/2846495)



  

