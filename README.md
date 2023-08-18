# NormalModes.jl

This is a simple package that provides utility function to compute normal modes from a Hessian matrix, and then work with them.

The other use of this package is to give me a space to ramble about normal modes, and I will shamelessly use this README for this purpose.

# A proper introduction of normal modes

Because I find the subject confusing, and the introductions describing the problem as well, I think it is worth laying down the basics here.

We are interested in the ground state of an equilibrium Hamiltonian $\mathcal{H} = \mathcal{T} + \mathcal{V}$, whith $\mathcal{T}$ the kinetic energy operator and $\mathcal{V}$ the potential energy operator. We note ${\rm \bf x}$ the vector of concatenated position, that is ${\rm \bf x} = {\rm vcat}({\rm \bf r}_1, {\rm \bf r}_2, {\rm \bf r}_3, \dots)$[^julia], with ${\rm \bf r}_A$ the position of atoms $A$.

Then we do a Taylor expansion of the potential, leading to

$$
\mathcal{V}({\rm \bf x}) = V_0 + \nabla \mathcal{V}({\rm \bf x}) + {\rm \bf x}^T {\rm \bf H} {\rm \bf x} + \mathcal{O}({\rm \bf x}^3),
$$

where ${\rm \bf H}$ is the Hessian matrix of the potential[^potential], and $V_0$ a constant shift that we can set to zero. Also since we are assuming to be at equilibrium, the gradient $\nabla \mathcal{V}({\rm \bf x})$ is zero as well. This leads to the approximation

$$
\mathcal{V}({\rm \bf x}) = {\rm \bf x}^T {\rm \bf H} {\rm \bf x}
$$

and the Hamiltonian

$$
\mathcal{H} = \sum_i \frac{-\hbar^2}{2 m_i} \frac{\partial^2}{\partial x_i^2} + {\rm \bf x}^T {\rm \bf H} {\rm \bf x},
$$

where $m_i$ is the mass of dimension $i$[^mass]. Now this is the Hamiltionian of a set of coupled harmonic oscillator, and we will change basis to rewrite it as a system of *uncoupled* harmonic oscillators.

We introduce the diagonal mass weighting matrix ${\rm \bf M}$ with ${\rm \bf M}_{ii} = m_i^{-\frac{1}{2}}$. We do an eigendecomposition of the weighted Hermitian $\rm \bf M H M$ as

$$
{\rm \bf M H M} = {\rm \bf U \Omega U}^T,
$$

where $\rm \bf U$ is the matrix of vector and $\rm \bf \Omega$ the diagonal matrix of eigenvalues, that we write ${\rm \bf \Omega}_{ii} = \omega_i^2$[^posdef].

We have what we need to define the new coordinates $\rm \bf z$ that we need, namely through

$$
{\rm \bf x} = {\rm \bf M U z}.
$$

Inserting it in the Hamiltonian we get the decoupled system of harmonic oscillators (see the annex for the derivation)

$$
\mathcal{H} = \sum_j \frac{-\hbar^2}{2} \frac{\partial^2}{\partial z_j^2} + \omega_j^2 z_j^2.
$$

We can check wikipedia to get whatever we need from it, like for example the Wigner distribution for sampling, and it is how it was implemented in the package.

# The package

Now what are the normal modes in all that? That's where it gets confusing and the reason I wrote this section to begin with. Candidates are

1. The eigen vectors ${\rm \bf u}_j$ that are the columns of $\rm \bf U$. They are however defined only in the mass weighted space.
2. The columns of $\rm \bf MU$. They are in real space, but they are are no longer neither normed nor orthogonal.
3. The columns of $\rm \bf MU$ normalized. Still not orthogonal to each other but at least they have norm 1. However it loses information about $\rm \bf M$ which is needed to convert between $\rm \bf x$ and $\rm \bf z$ or to perform Wigner sampling. So those are ofter return together with so-called *mode masses*, the norm of the columns of $\rm \bf MU$ before being normed[^lunacy].

As far as I know the last one is what, despite my strong feelings, is known as *normal modes*. However to avoid a slow descent into madness, the package provides a specific object type `NormalDecomposition(hessian)`, that stores the useful information and encapsulate the inner confusing parts. Available are
- `normal_modes` : The (normed) normal modes for example to plot or animate them (possible integration with Makie may become available at some point).
- `sample` : Perform Wigner sampling, and get the displacements from the mean for positions and momenta.
- `frequencies` : Frequencies of the modes (in GHz).
- `wave_number` : Wave numbers of the modes (in inverse cm).
- `reduced_mass` : Reduced masses of the modes (in AMU), following [this question](https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses) and (I believe) similar to what Gamess does.

# Appendix

## Derivation of the uncoupled equation

We want to substitute ${\rm \bf x} = {\rm \bf M U z}$ into

$$
\mathcal{H} = \sum_i \frac{-\hbar^2}{2 m_i} \frac{\partial^2}{\partial x_i^2} + {\rm \bf x}^T {\rm \bf H} {\rm \bf x}.
$$

By the definition of the diagonal matrix ${\rm \bf M}$ and the orthogonal matrix ${\rm \bf U}$ we have

$$
{\rm \bf x}^T {\rm \bf H x} =
    {\rm \bf z}^T {\rm \bf U}^T {\rm \bf M} {\rm \bf HMUz} = {\rm \bf z}^T {\rm \bf \Omega z},
$$

where ${\rm \bf \Omega}$ is diagonal.

For the other term, we use Leibnitz chain rule, starting with the first derivative only

$$
\frac{\partial}{\partial x_i} =
    \sum_j \frac{\partial z_j}{\partial x_i} \frac{\partial}{\partial z_j}.
$$

Next we use ${\rm \bf z} = {\rm \bf U}^T {\rm \bf M}^{-1} {\rm \bf x}$ to compute

$$
\frac{\partial z_j}{\partial x_i} =
    \frac{\partial}{\partial x_i} {\rm \bf e}_j^T {\rm \bf U}^T {\rm \bf M}^{-1} {\rm \bf x} =
    {\rm \bf e}_j^T {\rm \bf U}^T {\rm \bf M}^{-1} {\rm \bf e}_i,
$$

where ${\rm \bf e}_i$ are the unit vectors of the canonical basis. Since no new
dependency in ${\rm \bf x}$ or ${\rm \bf z}$ appear, we can apply it
twice. We transpose one of them however, to make the vector with index $i$
appear together[^transpose trick]

$$
\sum_i \frac{-\hbar^2}{2 m_i} \frac{\partial^2}{\partial x_i^2} =
    -\frac{\hbar^2}{2} \sum_{i, j, k} \frac{1}{m_i}
        \frac{\partial z_k}{\partial x_i}
        \left[\frac{\partial z_j}{\partial x_i} \right]^T
    = -\frac{\hbar^2}{2} \sum_{j, k}
        {\rm \bf e}_k^T {\rm \bf U}^T {\rm \bf M}^{-1}
        \left[ \sum_i \frac{1}{m_i} {\rm \bf e}_i {\rm \bf e}_i^T \right]
        {\rm \bf M}^{-1} {\rm \bf U} {\rm \bf e}_j
        \frac{\partial}{\partial z_j} \frac{\partial}{\partial z_k}.
$$

Now it would be really beautiful if

$$
\sum_i \frac{1}{m_i} {\rm \bf e}_i {\rm \bf e}_i^T = {\rm \bf M}^2,
$$

as it would make everything collapse. Thankfully with our choice of ${\rm \bf M}$ it is exactly what happens.

Putting it in and using the fact that ${\rm \bf e}_j^T {\rm \bf e}_k$ is a Kroenecker delta, we get, as expected,

$$
\sum_i \frac{-\hbar^2}{2 m_i} \frac{\partial^2}{\partial x_i^2}
    = -\frac{\hbar^2}{2} \sum_j \frac{\partial^2}{\partial z_j^2}.
$$

## Footnotes

[^julia]: I am indeed using the julia function `vcat` to describe this mathematical operation, I have never found a better way to express it.

[^potential]: This makes sense because in the position representation the potential is not an operator but just a scalar.

[^mass]: Which is the mass of the atom from which this component comes from. Essentially in the vector of mass $m_i$ each mass is repeated 3 times, once for each spatial dimension of each atom.

[^posdef]: Eigenvalue are in theory positive or zero, because the Hamiltonian is semi-positive definite and symmetric. Numerically the smallest eigenvalues can be negative instead of exactly zero.

[^annex]: The correct transformation of the derivative requires a bit of work to properly simplify, but trust me you can do it.

[^lunacy]: Which is pure lunacy if you ask me, the code I was using was norming the component in one function, and then immediately unnorming them in the next one. 

[^transpose trick] We transpose a scalar, which is a trick I like very much and
which is surprisingly useful.
