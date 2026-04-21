# BLUEPRINT — Kagome HOTI Simulation

Computational Blueprint
Breathing Kagome Higher-Order Topological Insulator Plates
Simulation-Only Phase A: Detailed Step-by-Step Plan with Verification Checklists
Principal Investigator: Raj Rahul Mandal
Faculty Advisor: Dr. Yang (Caldwell Summer Research Fellowship)
Target Venue: Physical Review Applied
Document Date: April 2026
1. Executive Summary and How to Read This Document
This document is the complete technical blueprint for the computational phase of a research project that will produce a first-authored paper on a printable audio-frequency breathing kagome higher-order topological insulator plate. It is written for one reader, the principal investigator, and it assumes you are starting with no existing code, no simulation results, and no prior FEniCSx experience. Every assumption you must test, every number you must verify, and every decision you must make before writing code is stated explicitly.
The document is organized in strict execution order. Section 2 states the research question, hypothesis, and scientific value. Section 3 lists every pre-code decision that is already locked. Section 4 explains the theoretical foundation in enough depth that you can verify your own derivations. Sections 5 through 10 are the six computational phases, each containing atomic steps with individual time estimates, verification criteria, and specific failure modes. Section 11 is the disorder protocol. Section 12 is the figure map for the paper. Section 13 is the complete milestone calendar at five hours per week. Section 14 is the risk register. Section 15 is the list of items to escalate to Dr. Yang versus decide independently.
Read this document once from start to finish before writing any code. Then use it as a reference during execution, returning to each section when you reach the corresponding phase. The document is deliberately dense because density here saves weeks of wasted work later.
1.1 What This Project Is, in One Paragraph
You will build the first open-source, fully three-dimensional, viscoelastically-damped, orthotropically-anisotropic computational framework for a printable audio-frequency breathing kagome higher-order topological insulator plate. You will use that framework to produce quantitative predictions for five specific physical behaviors: the bulk band gap location and width, the three zero-dimensional corner states pinned to 60-degree acute terminations, the Z3 Berry phase bulk invariant, the disorder tolerance threshold for corner state survival, and the non-monotonic gap width behavior under Zener viscoelastic damping. The primary deliverable is a peer-reviewed paper in Physical Review Applied. The secondary deliverable is a GitHub repository with documented FEniCSx code, Gmsh scripts, and post-processing pipelines that any future researcher can extend.
1.2 What Changed From the Original Proposal
The original proposal described a honeycomb valley Hall plate with experimental validation on a 3D printer. Two literature findings invalidated the original framing. First, Rosiek and collaborators demonstrated in 2023 that valley Hall topological waveguides do not actually provide protection against backscattering from realistic structural disorder, collapsing the central protection claim of valley-Hall designs when fabricated with short-range imperfections such as those produced by commodity FDM printing. Second, the breathing kagome higher-order topological insulator relies on crystalline C3 symmetry and a true quantized Z3 Berry phase bulk invariant, which survives microscopic fabrication disorder up to approximately 10 percent geometric perturbation. The project has therefore been reframed around the kagome HOTI and restricted to a pure computational contribution, with experimental validation deferred to future work.
2. Research Framing and Scientific Value
2.1 The Research Question
What is the quantitative behavior of a breathing kagome higher-order topological insulator plate when simulated with the realistic material and fabrication effects that would govern a desktop-FDM-printed realization, and what is the resulting defect tolerance threshold?
This question is narrow on purpose. It is answerable with a bounded set of simulations. It produces numbers rather than qualitative descriptions. It invites direct comparison with existing kagome HOTI literature in aluminum and laser-cut PMMA, while filling the specific gap where no one has yet quantified the design for viscoelastic polymer plates at audio frequencies.
2.2 The Hypothesis, Stated as Five Testable Predictions
The primary hypothesis is that a breathing kagome lattice implemented as a three-dimensional viscoelastic orthotropic thin plate using PLA-representative material parameters will exhibit the following five behaviors, each producing a specific publishable figure.
A topologically non-trivial bulk band gap centered between 3.5 and 4.5 kHz, with width between 500 Hz and 1.5 kHz, whose location and width are predictable from the bending stiffness and lattice constant via a simple analytical estimate from thin-plate theory.
Three zero-dimensional corner states pinned to the 60-degree acute terminations of an equilateral triangular finite sample, each with eigenfrequency within 10 percent of the bulk gap center and spatial inverse participation ratio exceeding 10 times that of bulk modes.
A Z3 Berry phase invariant that takes quantized values of zero in the trivial phase (breathing ratio less than one) and two pi over three in the topological phase (breathing ratio greater than one), with the transition precisely at breathing ratio equal to one.
Corner state survival probability exceeding 90 percent for geometric disorder standard deviations up to approximately 10 percent of the nominal beam width, followed by linear degradation at higher disorder strengths, quantified with at least 30 random realizations per disorder level.
Non-monotonic gap width behavior as a function of loss tangent under the Zener standard linear solid damping model, specifically a regime where the gap widens before eventually narrowing, which is invisible to the simpler Kelvin-Voigt damping model.
Each of these predictions corresponds to one figure in the paper. If all five are confirmed, the design is validated as a computational blueprint. If any one fails in an unexpected way, that itself is a publishable finding. There is no outcome in which the project produces no paper.
2.3 Why This Research Is Valuable
The scientific community has produced a rich literature on kagome HOTIs in aluminum, laser-cut PMMA, and silicon chips at megahertz frequencies. The specific case of a polymer plate at audio frequencies, printable on commodity desktop FDM hardware and including the full realistic physics of viscoelastic damping, orthotropic fabrication anisotropy, and three-dimensional vector elasticity, has not been treated as an integrated simulation framework anywhere in the published record. Closing this gap produces four distinct contributions.
First, a specific numerical prediction for the disorder tolerance threshold. Everyone cites the general robustness of HOTIs. Nobody has quantified where the breakdown happens for the specific case of an FDM-printable plate. This paper will report a percentage with an uncertainty estimate.
Second, a subtle damping physics result. The Zener versus Kelvin-Voigt comparison, specifically the branch overtaking phenomenon where increased damping can widen a gap in certain regimes, has been noted in viscoelastic phononics but has never been applied to HOTIs. This paper will be the first to report whether HOTI corner states benefit or suffer from the hereditary damping kernel.
Third, a reproducible open-source framework. Future researchers can swap PLA for PETG, swap kagome for another HOTI lattice, swap the triangle for a different geometry, and rerun the pipeline. Papers that release working code of this quality get cited for years.
Fourth, positioning for experimental follow-up. If a later group decides to print and measure the plate, the simulation tells them exactly what to expect. The computational blueprint stands as a standalone contribution regardless of whether experiment follows.
3. Locked Design Parameters
Every parameter below is locked before any code is written. These are not optimization variables. They are design decisions that set the physical regime and cannot be changed mid-project without restarting most of the phases. If Dr. Yang overrides one of these, the timeline extends by two to four weeks.
3.1 Geometric Parameters

3.2 Material Parameters (PLA)

3.3 Target Output
The target gap center frequency is 3.5 to 4.5 kHz. The target gap width is 500 Hz to 1.5 kHz, roughly 15 to 30 percent of center. The three corner states are expected at or within 10 percent of gap center. The Z3 Berry phase in the topological phase is exactly two pi over three. The disorder tolerance is expected to be 10 to 15 percent. Each of these is a testable prediction, and each figure in the paper demonstrates one of them.
4. Theoretical Foundation You Must Verify Before Coding
This section is not a replacement for reading the primary literature. It is a derivation trace that tells you what you must be able to reproduce on paper before you trust any simulation output. If you cannot derive the results below from scratch, your simulation outputs will be numbers you cannot interpret.
4.1 The Breathing Kagome Tight-Binding Hamiltonian
The kagome lattice has three sites per primitive unit cell, arranged at the corners of an upward-pointing triangle and its inverted neighbor. In the breathing kagome, bonds within the upward triangles (intra-cell) have hopping amplitude t_a, and bonds connecting to the inverted triangles (inter-cell) have amplitude t_b. When t_a equals t_b, the lattice is standard kagome with a Dirac cone at the K point and a flat band. When t_a does not equal t_b, the flat band remains but the Dirac cone is gapped.
The tight-binding Hamiltonian in momentum space is a three-by-three matrix whose off-diagonal elements carry the hopping amplitudes with Bloch phase factors. The topological phase transition occurs at t_a equal to t_b. For t_b greater than t_a the system is in a topological phase with corner states at acute terminations. For t_a greater than t_b the system is trivial.
Verification task before coding: derive this three-by-three Hamiltonian on paper, diagonalize it at the Gamma and K points for three ratios (0.5, 1.0, 1.5), and plot the resulting bands using NumPy in a twenty-line Python script. Confirm the gap closes at ratio 1.0 and opens symmetrically on each side.
4.2 The Z3 Berry Phase as the Bulk Invariant
The relevant topological invariant for the breathing kagome HOTI is the Z3 Berry phase, computed as a Wilson loop along a closed path in the Brillouin zone. For a C3-symmetric system, the Berry phase is quantized to values of zero or plus or minus two pi over three. The quantization is protected by C3 rotational symmetry alone, without requiring time-reversal or inversion.
The Berry phase is extracted numerically by computing overlaps of Bloch eigenvectors at neighboring k-points along a path, multiplying these link variables around a closed loop, and taking the imaginary part of the logarithm. Gauge fixing is handled automatically because the loop product is gauge invariant.
Verification task: implement the Wilson loop Berry phase for your tight-binding NumPy toy model. Confirm it gives zero in the trivial phase and two pi over three in the topological phase. Verify convergence as you refine the k-mesh from 10 by 10 to 50 by 50. This is the reference against which you will later validate your FEniCSx implementation.
4.3 The Bulk-Boundary Correspondence for HOTIs
Unlike first-order topological insulators where a bulk invariant dictates the existence of edge states, HOTIs exhibit a higher-order bulk-boundary correspondence where the bulk invariant dictates the existence of corner states. For the breathing kagome, three corner states exist whenever the bulk is in the topological phase and the sample is terminated with three 60-degree acute angles. The corner states carry fractional charge one-third each, summing to unity, reflecting the Z3 nature of the topology.
Critically, corner states do not exist at 120-degree obtuse terminations. A hexagonal sample of the topological kagome has six 120-degree corners and zero corner states. An equilateral triangle has three 60-degree corners and three corner states. A parallelogram has two of each, yielding two corner states. Geometry matters.
Verification task: read Ezawa 2018 Physical Review Letters on the breathing kagome HOTI. Understand the proof that the Z3 Berry phase implies corner states at acute terminations. You do not need to reproduce the proof, but you must be able to explain in plain language why the corner count depends on geometry.
4.4 Translation to Continuum Elasticity
The tight-binding model captures the topology but not the full quantitative behavior of a continuum elastic plate. In the continuum description, each tight-binding site becomes a node of the beam network, the hopping amplitudes become effective spring constants determined by beam flexural stiffness, and the eigenvalues of the Hamiltonian become squared angular frequencies omega squared.
The effective spring constant scales as width cubed over length cubed for flexural bending of a beam, so the breathing ratio of spring constants scales as the cube of the beam width ratio. A geometric breathing ratio of 1.5 in beam width therefore corresponds to a dynamical breathing ratio of 3.375 in effective spring constant. This is well into the topological phase.
Verification task: derive the effective spring constant formula for a doubly-clamped Euler-Bernoulli beam in flexure. Confirm the cubic scaling. Use this to predict the dynamical breathing ratio from the geometric breathing ratio, and verify that your chosen nominal geometry puts you unambiguously in the topological phase.
4.5 The Flexural Wave Dispersion Relation
Flexural waves in a thin isotropic plate obey the biharmonic Kirchhoff-Love dispersion relation, omega squared equals D over rho h times k to the fourth, where D is the bending stiffness and k is the wavenumber. The quadratic scaling of frequency with wavenumber means that halving the lattice constant quadruples the gap frequency.
For PLA at 3 mm thickness, plugging in E equal to 3 GPa, rho equal to 1240 kg per cubic meter, and nu equal to 0.36, the dispersion at the K-point of a 4 cm lattice gives a rough frequency estimate of 3.5 to 4.5 kHz. This matches the target and gives you a sanity-check formula to verify FEniCSx output against.
Verification task: compute this dispersion analytically. Produce a single-cell NumPy prediction of the lowest flexural band at the K-point. When FEniCSx eventually produces a band structure, the K-point frequency must match this analytical prediction to within 20 percent.
5. Phase 2: Environment Setup and Tight-Binding Sanity Check
Timeline: Week 1. Budget: 5 hours. Deliverable: working conda environment and a NumPy tight-binding script that reproduces the breathing kagome band structure and Z3 Berry phase. This phase exists to catch environment errors and conceptual errors before you commit to FEniCSx.
5.1 Environment Setup
You are starting from a Linux or macOS machine (Raspberry Pi 5 is usable for this phase, though the main simulation phases will require the campus CUDA GPU workstation).
Atomic steps
□  Install Miniforge or Mambaforge. Miniconda works but conda-forge-pinned environments resolve faster with mamba.
□  Create a conda environment named kagome-hoti pinned to Python 3.11.
□  Install fenics-dolfinx, gmsh, mpi4py, petsc4py, and slepc4py from conda-forge.
□  Install numpy, scipy, matplotlib, h5py, pyvista, and meshio via pip within the conda environment.
□  Verify installation by running the FEniCSx hello world example from the official tutorial.
□  Verify Gmsh installation by scripting a simple 2D square mesh and viewing it.
□  Set up a private GitHub repository named kagome-hoti-simulation.
□  Create the directory structure: src/, scripts/, meshes/, results/, figures/, notebooks/, tests/, docs/.
5.2 Tight-Binding Sanity Check in NumPy
Before writing FEniCSx code, implement the pure tight-binding breathing kagome in forty lines of NumPy. This is your ground truth reference. Every subsequent simulation must reproduce these qualitative features.
Atomic steps
□  Write a function that returns the 3x3 breathing kagome Bloch Hamiltonian given k_x, k_y, t_a, t_b.
□  Diagonalize the Hamiltonian along Gamma-M-K-Gamma for ratios 0.5, 1.0, 1.5.
□  Plot the three bands with matplotlib. Verify Dirac cone at K for ratio 1.0.
□  Verify symmetric gap opening for ratios 0.5 and 1.5 (same gap width, opposite topology).
□  Implement the Z3 Berry phase via discretized Wilson loop on a 20x20 k-mesh.
□  Verify Berry phase equals 0 for ratio 0.5 (trivial).
□  Verify Berry phase equals 2*pi/3 for ratio 1.5 (topological).
□  Test convergence as you refine the k-mesh from 20x20 to 50x50.
□  Commit this working script to the repo as tests/tb_reference.py.
5.3 Verification Criteria (Phase Gate)
You do not proceed to Phase 3 until all of the following are true.
Conda environment activates cleanly on both your development machine and the campus GPU workstation.
FEniCSx hello world runs to completion without errors.
Gmsh produces and displays a valid 2D mesh.
The tight-binding band structure shows a closed Dirac cone at ratio 1.0.
The tight-binding band structure shows symmetric gaps at ratios 0.5 and 1.5.
The Z3 Berry phase computation gives 0 and 2*pi/3 to within 1 percent numerical error on a 50x50 k-mesh.
The repository is initialized with the directory structure and first commit.
5.4 Failure Modes and Fixes
If the conda environment fails to resolve, switch to mamba, which handles FEniCSx dependencies more reliably. If FEniCSx hello world fails with an MPI error, your mpi4py and petsc4py are out of sync; rebuild the environment from a clean state. If your Berry phase does not quantize to 2*pi/3, the most likely cause is a sign error in the Wilson loop product or an insufficient k-mesh density. If your band structure shows the Dirac cone closing at a ratio other than 1.0, you have a bug in the Hamiltonian matrix elements, not a physics surprise.
6. Phase 3: 3D Unit Cell Bloch Solver in FEniCSx
Timeline: Weeks 2 through 6. Budget: 20 hours over four weeks. Deliverable: working FEniCSx Bloch periodic solver that reproduces the tight-binding band structure qualitatively in full 3D continuum elasticity, sweeps breathing ratio across the topological transition, and confirms gap closure and reopening.
6.1 Parametric Gmsh Mesh of the 3D Unit Cell
The unit cell is a rhombic prism with two lattice vectors of length a at 60 degrees and a thickness h. It contains the beam network forming the breathing kagome, with two distinct beam widths for intra-cell and inter-cell connections. The mesh must resolve the flexural wavelength at the highest frequency you care about, which is roughly 2x the upper band edge.
Atomic steps
□  Write a Gmsh Python API script that takes (a, h, w_intra, w_inter, breathing_ratio) as parameters and produces the 3D unit cell.
□  Include the node-to-node beam network with both intra-cell and inter-cell struts.
□  Set mesh size to 2 mm nominal, refined to 1 mm at beam junctions where stress concentrates.
□  Tag lattice vector boundaries (pairs +a1/-a1, +a2/-a2) for periodic boundary conditions.
□  Export as .msh version 2.2 and verify with Gmsh GUI that the mesh looks correct.
□  Convert .msh to .xdmf using meshio for FEniCSx loading.
□  Count mesh elements. Target is 30,000 to 80,000 tetrahedral elements per unit cell.
6.2 FEniCSx Bloch Periodic Boundary Conditions
Bloch periodicity imposes u at x plus lattice vector equal to exp(i k dot lattice vector) times u at x. In FEniCSx 0.7 this is implemented via MultiPointConstraint or via a custom form assembly. The k-dependence of the phase factor means you assemble a new stiffness matrix for every k-point, which is the dominant computational cost.
Atomic steps
□  Define the VectorFunctionSpace with second-order Lagrange elements for displacement.
□  Write the linear elasticity weak form with isotropic constitutive tensor.
□  Implement Bloch periodic boundary conditions via MultiPointConstraint, linking paired boundary DOFs with k-dependent phase.
□  Assemble K and M matrices at k = Gamma point.
□  Solve the generalized eigenvalue problem with SLEPc, requesting the lowest 30 eigenvalues.
□  Verify the lowest 3 modes are zero-frequency rigid body modes (translations).
□  Sweep k along Gamma-M-K-Gamma with 60 points total, solving at each k.
□  Plot the band structure and compare with the tight-binding reference.
6.3 Verify the Topological Transition
Atomic steps
□  Run the band structure at breathing ratio 1.0 (standard kagome). Verify Dirac cone at K closes.
□  Run at ratio 0.5 (trivial). Verify gap opens.
□  Run at ratio 1.5 (topological). Verify gap opens with similar width as trivial.
□  Compute gap center frequency at ratio 1.5. Compare with analytical thin-plate prediction. Must match within 20 percent.
□  Color-code band structure by polarization (flexural vs extensional vs shear horizontal).
□  Identify the polarization of the bands surrounding the gap. Must be predominantly flexural.
6.4 Verification Criteria (Phase Gate)
Mesh has 30,000 to 80,000 tetrahedra per unit cell.
Lowest three modes at Gamma are zero frequency to within 0.01 Hz (rigid body).
Dirac cone closes cleanly at breathing ratio 1.0.
Gap opens at ratio 0.5 and ratio 1.5 with widths within 10 percent of each other.
K-point frequency matches analytical thin-plate prediction to within 20 percent.
Gap bands are predominantly flexural, not extensional or shear.
Total wall-clock time per full band structure is less than 30 minutes on the GPU workstation.
6.5 Failure Modes and Fixes
The single most common failure in this phase is a buggy Bloch periodic boundary condition implementation. Symptom: bands look qualitatively wrong, or the Dirac cone closes at a non-physical location. Fix: verify PBCs with a test problem where the expected answer is known, such as a homogeneous cube with known free vibration modes. Debug before trusting kagome output.
The second failure mode is an undersized mesh. Symptom: gap width is wrong by 30 percent or more, eigenvalues are erratic between neighboring k-points. Fix: refine mesh, especially at beam junctions. Doubling mesh density and observing less than 5 percent change in key eigenvalues is convergence.
The third is eigensolver convergence failure. Symptom: SLEPc returns fewer than 30 converged eigenvalues or produces spurious modes. Fix: switch from shift-invert with sigma=0 to Krylov-Schur with a target frequency near the expected gap, and loosen the tolerance to 1e-8.
7. Phase 4: Viscoelastic Damping Sweep
Timeline: Weeks 6 through 9. Budget: 15 hours over three weeks. Deliverable: gap width vs loss tangent plot for both Kelvin-Voigt and Zener damping models, revealing the branch overtaking phenomenon that only Zener captures.
7.1 Implementing Complex Modulus
Viscoelastic damping is introduced by replacing the real Young's modulus E with a complex modulus E* = E'(1 + i eta(omega)), where eta is the loss tangent. For Kelvin-Voigt, eta is linear in frequency, eta = eta_0 * omega. For Zener (standard linear solid), the loss tangent has a Debye peak at frequency 1/tau, where tau is the relaxation time.
The eigenvalue problem becomes complex-valued. Real eigenfrequencies become complex: the real part is the oscillation frequency, the imaginary part is the decay rate. The gap is now defined as the frequency interval where no real eigenfrequencies exist and where imaginary parts are small enough that any modes near that range are heavily damped.
7.2 Atomic Steps
□  Modify the stiffness matrix assembly to use E* = E' * (1 + i * eta) for Kelvin-Voigt.
□  Implement complex-valued SLEPc eigensolver configuration.
□  Solve at loss tangent 0.0 as reference. Verify match with lossless Phase 3 result.
□  Sweep loss tangent from 0.0 to 0.15 in steps of 0.01 at Kelvin-Voigt.
□  Extract real eigenfrequencies, compute effective gap width at each loss tangent.
□  Plot gap width vs loss tangent for Kelvin-Voigt. Expect monotonic decrease.
□  Implement the Zener SLS constitutive model with parameters tau and E_infinity/E_0.
□  Fit tau from published PLA DMA data using Williams-Landel-Ferry extrapolation to audio frequencies.
□  Sweep loss tangent from 0.0 to 0.15 at Zener SLS.
□  Plot Zener gap width vs loss tangent on the same axes as Kelvin-Voigt.
□  Identify any branch-overtaking regime where Zener gap width exceeds Kelvin-Voigt.
7.3 Verification Criteria
Kelvin-Voigt gap width decreases monotonically with loss tangent.
Zener gap width is non-monotonic in some regime, ideally showing a local maximum or plateau.
At nominal PLA loss tangent 0.02, gap width is reduced by less than 30 percent from the lossless case.
Critical loss tangent at which gap closes entirely is identified and reported with uncertainty.
7.4 Failure Modes and Fixes
If SLEPc returns purely real eigenvalues at nonzero loss tangent, the complex modulus is not being applied correctly. Verify the stiffness matrix is actually complex by checking a few entries. If no branch overtaking appears in the Zener sweep, the relaxation time tau may be outside the audio frequency range; adjust tau so 1/tau falls within 3-5 kHz. If the Williams-Landel-Ferry extrapolation is unreliable, report the damping sweep for a range of plausible tau values as a sensitivity study.
8. Phase 5: Z3 Berry Phase from FEniCSx Eigenvectors
Timeline: Weeks 9 through 11. Budget: 10 hours over two weeks. Deliverable: quantized Z3 Berry phase vs breathing ratio plot, verifying the topological transition at ratio 1.0.
8.1 Wilson Loop Implementation
You will port the NumPy tight-binding Wilson loop from Phase 2 to operate on FEniCSx continuum eigenvectors. The core algorithm is unchanged: compute overlaps between Bloch eigenvectors at neighboring k-points, multiply around a closed loop, take the imaginary part of the logarithm.
The subtlety is that FEniCSx eigenvectors are discretized displacement fields on a finite element mesh, not abstract vectors in a fixed basis. Overlaps require integrating over the unit cell with the mass matrix as the inner product weight: <u1 | M | u2>. Normalization and phase conventions must be enforced consistently across k-points.
8.2 Atomic Steps
□  Set up a 30x30 k-mesh covering the full Brillouin zone.
□  Solve the eigenvalue problem at each k-point, extracting the lowest three bands.
□  Normalize each eigenvector against the mass matrix: v_normalized = v / sqrt(v^dagger M v).
□  Compute link variables U_k->k+dk = <u(k) | M | u(k+dk)> / |<u(k) | M | u(k+dk)>|.
□  Compute the Wilson loop product around the boundary of the reduced Brillouin zone.
□  Take the imaginary part of the logarithm to extract the Berry phase.
□  Verify against the tight-binding reference from Phase 2.
□  Sweep breathing ratio from 0.5 to 2.0 in steps of 0.1.
□  Plot Berry phase vs breathing ratio. Expect a step from 0 to 2*pi/3 at ratio 1.0.
□  Refine k-mesh from 30x30 to 60x60 at three values to confirm quantization.
8.3 Verification Criteria
Berry phase is 0 plus or minus 0.05 radian for all ratios below 1.0.
Berry phase is 2*pi/3 plus or minus 0.05 radian for all ratios above 1.0.
Transition at ratio 1.0 is sharp (jump within ratio 1.0 plus or minus 0.02).
Refining the k-mesh improves the quantization plateau rather than shifting it.
8.4 Failure Modes and Fixes
The largest failure mode is an inconsistent phase convention across k-points, which produces a random-looking Berry phase rather than a quantized value. Fix: enforce a convention where the largest-magnitude component of each eigenvector is made real positive before use. If the transition is fuzzy rather than sharp, the k-mesh is too coarse near the band crossing; refine to 60x60. If the Berry phase gives a sign opposite to expected, re-verify your lattice vector sign convention; this is cosmetic and not a physics error.
9. Phase 6: Finite Triangular Sample and Corner States
Timeline: Weeks 11 through 15. Budget: 20 hours over four weeks. Deliverable: three localized corner states in a 12-cell equilateral triangle sample, visualized in real space, with eigenfrequencies sitting inside the bulk gap.
9.1 Building the Finite Triangle Mesh
The 12-by-12 equilateral triangle has 78 unit cells total, which at 3 sites per cell is 234 network nodes. The full 3D mesh has roughly 1.5 to 2 million degrees of freedom. This is the largest simulation in the project and will stress desktop RAM. Use the campus workstation.
Atomic steps
□  Modify Gmsh script to produce the 12-cell equilateral triangle cut from the breathing kagome lattice.
□  Ensure the three corners terminate at 60-degree acute vertices. Verify visually.
□  Apply free (stress-free) boundary conditions on all external edges.
□  Use coarser mesh density (3-4 mm) than the unit cell to keep DOF count manageable.
□  Count total DOFs. Target is 1 to 2 million.
□  Verify mesh on a paper printout showing all three corners are clearly acute.
9.2 Solving for Eigenmodes
Atomic steps
□  Assemble K and M matrices for the finite sample (no Bloch periodicity).
□  Use SLEPc shift-invert with target frequency at the expected gap center (~4 kHz).
□  Request the lowest 100 eigenvalues around the target.
□  Identify modes with eigenfrequency inside the bulk gap identified in Phase 3.
□  Filter by spatial localization: compute inverse participation ratio for each mode.
□  Identify three modes with IPR at least 10x greater than bulk modes.
□  Plot the three localized modes in 3D. Each must peak at one of the three corners.
9.3 Verifying Corner State Properties
Atomic steps
□  Extract eigenfrequencies of the three corner states. All three should be within 5 percent of each other (degenerate in infinite-sample limit, slightly split by finite-size effects).
□  Verify each mode's peak displacement is at a different corner, matching the three 60-degree vertices.
□  Compute decay length by fitting an exponential to displacement magnitude vs distance from the corner.
□  Decay length should be 1 to 3 unit cells.
□  Compare with the trivial case (breathing ratio 0.5): verify no corner states exist there.
9.4 Verification Criteria
Exactly three modes sit inside the bulk gap for breathing ratio 1.5.
Each mode is localized at a different corner of the triangle.
Each mode has IPR at least 10 times the bulk average.
Decay length is 1 to 3 unit cells.
No equivalent modes exist in the trivial phase (breathing ratio 0.5).
9.5 Failure Modes and Fixes
Shift-invert in SLEPc can be fragile for 2-million-DOF systems. If it fails, switch to the polynomial filtered eigensolver or the contour integral eigensolver. If too few modes are found in the gap, increase the number of requested eigenvalues. If modes in the gap are delocalized, check the mesh at the corners; a sloppy mesh termination can create spurious edge modes. Always compare to the trivial phase as a control.
10. Phase 7: Disorder Sweep and Defect Tolerance
Timeline: Weeks 15 through 20. Budget: 25 hours over five weeks. Deliverable: corner state survival probability as a function of disorder strength, with 30 random realizations per level, yielding a publishable threshold figure.
10.1 Disorder Protocol
Three forms of disorder are physically meaningful for FDM-printed plates: node position perturbations (corresponding to positional errors in the deposited material), beam width perturbations (corresponding to extrusion width variance), and bond removals (corresponding to print failures). Start with beam width perturbations because they are the cleanest to implement and the most representative of real FDM variance.
Atomic steps
□  Parameterize the Gmsh script to accept per-beam width values drawn from a Gaussian distribution.
□  Set up disorder levels at 0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25 percent of mean beam width.
□  For each disorder level, generate 30 independent random realizations (different random seeds).
□  This gives 9 levels times 30 realizations = 270 distinct meshes.
□  Solve for eigenmodes in each; extract corner state count, frequencies, IPRs.
10.2 Survival Statistics
Atomic steps
□  Define corner state survival: mode exists in the gap, localized at one of the three corners, with IPR at least 5x bulk average.
□  Count surviving corner states per realization (0, 1, 2, or 3).
□  Compute survival probability per level = (total surviving corners) / (3 * 30).
□  Plot survival probability vs disorder level with error bars (binomial uncertainty).
□  Fit a logistic or piecewise linear curve to identify the threshold where survival drops below 50 percent.
□  Plot corner state frequency mean and standard deviation vs disorder level.
□  Frequency standard deviation should grow roughly linearly with disorder.
10.3 Computational Efficiency
270 full 3D finite-sample solves at roughly 10 minutes each is 45 hours of wall-clock compute. This must run on the campus GPU workstation, overnight and over weekends, in batches of 30 with a shell script. Parallelize across MPI ranks if the workstation allows.
Atomic steps
□  Write a shell script that iterates through 270 parameter files, calls the FEniCSx solver, and saves eigenfrequencies and IPRs to a CSV.
□  Run locally on 3 realizations at a small disorder level as a smoke test. Confirm CSV output is correct.
□  Launch the full sweep on the GPU workstation. Estimate total runtime.
□  Monitor for crashes and re-run failed cases.
10.4 Verification Criteria
Survival probability at 0 disorder is 100 percent (confirming no bugs).
Survival probability at 10 percent disorder exceeds 85 percent.
Survival probability drops to below 50 percent at some threshold between 10 and 20 percent.
Frequency spread grows roughly linearly with disorder.
The threshold and its uncertainty are both reportable numbers.
10.5 Failure Modes and Fixes
If survival is 100 percent even at 25 percent disorder, either your disorder implementation is buggy (verify by visually inspecting a disordered mesh) or your threshold criterion is too permissive. If survival drops rapidly even at 2.5 percent, either your baseline corner states are fragile (check the nominal case) or your threshold criterion is too strict. Expect to iterate on the survival criterion before the final figure is clean.
11. Paper Figure Map
The paper has six figures, each tied to one phase. Section 11 is the figure budget and the data each figure must contain.

12. Milestone Calendar at 5 Hours per Week
The full computational project runs approximately 25 weeks, or about 6 months, from project start. At 5 hours per week this is 125 hours total, matching the original estimate of 120 to 130 hours for the simulation-only scope. Start date April 2026, target submission late 2026.

13. Risk Register
Four categories of risk that can derail this project. Each is listed with probability, impact, and mitigation.
13.1 Technical Risks
FEniCSx Bloch periodic boundary conditions implementation bugs. Probability high, impact high. Mitigation: verify PBCs on a homogeneous test cube before kagome, compare systematically against the NumPy tight-binding reference, allocate a full week of buffer in Phase 3.
SLEPc eigensolver convergence failures on the 2M-DOF finite sample. Probability medium, impact high. Mitigation: start with coarser mesh, use shift-invert with known-good target frequency, fall back to contour integral method if shift-invert fails.
Disorder sweep wall-clock time exceeds expectations. Probability medium, impact medium. Mitigation: if 270 realizations each at 10 minutes is too long, reduce sample size to 10x10 cells (cuts DOFs by half) or reduce realizations per level to 20.
13.2 Physics Risks
Corner states do not appear cleanly at the chosen parameters. Probability low, impact very high. Mitigation: this is what Phase 6 tests. If it fails, the problem is usually that the breathing ratio is too close to 1.0 and finite-size effects hybridize corner and bulk modes. Increase ratio to 1.8 or sample to 15 cells per side.
Gap is overwhelmed by damping at nominal PLA loss tangent. Probability low, impact high. Mitigation: the literature suggests the gap survives at 2 percent loss tangent. If it does not, the project pivots to report the critical damping as the key result, which is still publishable.
Zener branch overtaking does not materialize. Probability medium, impact low. Mitigation: this is a bonus finding, not a core result. If it does not appear, simply report gap-vs-damping for both models without the headline of branch overtaking.
13.3 Schedule Risks
Five hours per week not sustained due to semester or personal pressure. Probability high, impact high. Mitigation: Phase 3 and Phase 6 have continuity sensitivity. Avoid scheduling these phases over winter break or finals. If weeks are missed, pause at a phase boundary rather than mid-phase.
Dr. Yang unavailable at critical review points. Probability medium, impact medium. Mitigation: set a standing bi-weekly meeting now. Make sure the Phase 3 and Phase 6 gate reviews fall on weeks when Dr. Yang is known to be in town.
13.4 Scoop Risks
Another group publishes similar work during the project. Probability low, impact high. Mitigation: monitor arxiv weekly for 'kagome HOTI plate' and 'viscoelastic higher-order topological'. If scoop occurs, the branch overtaking and the specific disorder threshold are still novel contributions. Reframe rather than abandon.
14. What to Escalate to Dr. Yang
This section separates decisions that belong to Dr. Yang from decisions you make independently. Confusion on this wastes advisor time and delays the project.
14.1 Decisions That Require Dr. Yang's Approval
Final choice of target journal (PRApplied vs Journal of Sound and Vibration).
Any change to the fundamental project scope (adding experiment, switching lattice, adding collaborators).
Authorship order (default: you first, Dr. Yang corresponding, but confirm).
Any decision to release code under a restrictive license rather than MIT or BSD.
Final approval of the submitted paper.
14.2 Decisions You Make Independently
All code structure, naming, and implementation details.
Numerical tolerances, mesh sizes, solver choices.
Parameter sweeps within the locked ranges in Section 3.
Figure aesthetics (colors, fonts, layout).
When to reach out for help with specific bugs.
14.3 Decisions That Need a Check-In but Not Approval
Interpretation of the damping sweep results (is the branch overtaking real?).
Reporting conventions in the paper (how to describe uncertainty, how to cite prior work).
Choice of additional disorder protocols if beam-width disorder alone is not publishable.
15. How to Use This Document
This blueprint is a living document. Treat it as follows.
Before starting any phase, reread Section 3 to remind yourself of the locked parameters, then read the section for that phase in full. Work through the atomic steps in order. Check each box as you complete it. Do not skip a verification criterion even if it feels redundant. The verification criteria exist because projects like this fail silently when unverified assumptions accumulate.
When you encounter an unexpected result, first check whether it matches the 'failure modes and fixes' subsection for the current phase. Roughly 80 percent of problems you will encounter are already anticipated in those subsections. If the issue is not listed, add it to the document with a dated note. The document improves as you execute.
Before each meeting with Dr. Yang, scan Section 14 to confirm which items you need her approval on. Bring only those items. Questions that fall under independent decisions should be solved yourself before the meeting.
At the end of each phase, update Section 12 (milestone calendar) with your actual hours and adjust the remaining schedule. If you are ahead, do not advance a phase; use the time for verification. If you are behind, assess whether a sub-phase can be simplified rather than skipped.
15.1 Daily Workflow in Claude Code
Your primary tool is Claude Code with Opus 4.7. Use the following pattern every session.
Open Claude Code with the kagome-hoti-simulation repo as the working directory.
State which phase and which atomic step you are working on.
Ask Claude to either write, debug, or explain code; be specific.
Review generated code before running. Do not autorun in production.
After each atomic step is verified, commit to git with a descriptive message.
At the end of the session, update your lab notebook (notebooks/log.md) with what you did, what worked, what failed, and what the next step is.
15.2 Communication Discipline
Send Dr. Yang a short status email every two weeks, regardless of whether you are ahead or behind. The email contains three lines: what you finished, what you are working on, what you need. Nothing else. This maintains visibility without consuming her time.
15.3 When to Stop and Re-Plan
Stop execution and revisit this document if any of the following is true. One: you have fallen three or more weeks behind the calendar in Section 12. Two: a verification criterion in any phase fails and the failure mode in the phase's troubleshooting subsection does not describe your problem. Three: the scope of the paper shifts, for instance if an experimental collaboration opportunity arises. Four: Dr. Yang requests a fundamental change.
This document is your complete blueprint. Everything from environment setup through submission is captured here. The next action is to open Claude Code, create the kagome-hoti-simulation repository, and execute Phase 2 Step 5.1.1: install miniforge. The day you do that, the project starts.
