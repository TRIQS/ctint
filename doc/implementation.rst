***********************************
Derivation and implementation notes
***********************************

The CTINT algorithm
===================

The following section describes the CTINT algorithm and its derivation.
We first present it for the case of static density-density interactions,
and present the extension to general dynamical interactions thereafter.
For notations and conventions used throughout these notes, the reader is
referred to the attached document (which we will refer to as “vertex
conventions”).

.. _static_dens:

Action with static density-density interaction
----------------------------------------------

Let us begin by considering the example of an impurity system with
static density-density interaction. The action can then be written as

.. math::

   \mathcal{S}=\mathcal{S}_0 + S_{\mathrm{int}}

with

.. math::

   \mathcal{S}_0 = -\bar{c}_{\bar{\alpha}} \mathcal{G}_{\bar{\alpha}\beta}^{-1} c_\beta

and

.. math::

   \mathcal{S}_{\mathrm{int}} = \sum_{ab} U_{ab} \int_\tau n_a(\tau) n_b(\tau).

Interaction expansion
---------------------

In CT-INT we proceed by expanding the partition function

.. math::
   :label: partition_function_def

   Z\equiv\int\mathcal{D}[\bar{c},c]e^{\bar{c}_{\bar{\alpha}}\mathcal{G}_{\bar{\alpha}\beta}^{-1}c_\beta-\mathcal{S}_{\mathrm{int}}}

in powers of :math:`\mathcal{S}_{\mathrm{int}}`. One obtains

.. math::
   :label: Z_full

   Z=Z_0\sum_{k=0}^\infty\frac{(-1)^k}{k!}\sum_{\mathbf{a},\mathbf{b}}\left(\prod_{i=1}^k U_{a_i b_i}\right)\int_{\boldsymbol{\tau}}\langle n_{a_1}(\tau_1)\dots n_{a_k}(\tau_k)n_{b_1}(\tau_1)\dots n_{b_k}(\tau_k)\rangle_0,

where here, and in the following, bold symbols denote vectors. Note that
in practice one typically orders the integral times, which allows to
cancel the :math:`1/k!` prefactor. Using the Wick theorem we can write
the non-interacting expectation value in the equation above compactly
using a determinant

.. math::

   Z=Z_0\sum_{k=0}^\infty\frac{(-1)^k}{k!}\sum_{\mathbf{a},\mathbf{b}}\left(\prod_{i=1}^k U_{a_i b_i}\right)\int_{\boldsymbol{\tau}}\det\hat{\mathcal{G}}_k,

Here we have introduced the :math:`2k\times2k` matrix

.. math::

   [\hat{\mathcal{G}}_k]_{ij}=\mathcal{G}_{x_i\bar{y}_j}.

with vectors

.. math::

   \boldsymbol{x}[\mathbf{a},\mathbf{b},\boldsymbol{\tau}] & \equiv[(c_1,\tau_1),\ldots,(c_{2k},\tau_{2k})]\equiv[(a_1,\tau_1),\ldots,(a_k,\tau_k),(b_1,\tau_1),\ldots,(b_k,\tau_k)],\\
   \bar{\boldsymbol{y}}[\mathbf{a},\mathbf{b},\boldsymbol{\tau}] & \equiv[(\bar{c}_1,\tau_1),\ldots,(\bar{c}_{2k},\tau_{2k})]\equiv[(a_1,\tau_1),\ldots,(a_k,\tau_k),(b_1,\tau_1),\ldots,(b_k,\tau_k)].

In summary we can write

.. math:: Z = Z_0\Xi
  :label: Z_as_Z0Xi

with

.. math::
   :label: Xi_def

   \Xi\equiv\sum_{k=0}^\infty\sum_{\mathbf{a},\mathbf{b}}\int_{\boldsymbol{\tau}}A_k[\mathbf{a},\mathbf{b}]\times\det\hat{\mathcal{G}}_k[\mathbf{a},\mathbf{b},\boldsymbol{\tau}]

and

.. math::

   A_k[\mathbf{a},\mathbf{b}]=\frac{(-1)^k}{k!}\left(\prod_{i=1}^k U_{a_i b_i}\right)

Note that :math:`A_k` is dependent on the interaction :math:`U`, while
the determinant is dependent on the bare propagator :math:`\mathcal{G}`.

The :math:`\alpha` function and the auxiliary spin
--------------------------------------------------

In order to handle the trivial sign problem originating from the
minus-signs associated with the expansion of each static interaction
vertex we shift the quartic part of the action
:math:`\mathcal{S}_{\mathrm{int}}\rightarrow\mathcal{\tilde{S}}_{\mathrm{int}}`
with

.. math::

   \tilde{\mathcal{S}}_{\mathrm{int}}
   &= \sum_l U_l \int_\tau
   \bar{c}_{\bar{a}_l}(\tau) c_{b_l}(\tau) \bar{c}_{\bar{c}_l}(\tau) c_{d_l}(\tau) -
   n_{\bar{a}_l b_l}(\tau) \alpha_{11}^l - n_{\bar{c}_l d_l}(\tau) \alpha_{00}^l +
   n_{\bar{a}_l d_l}(\tau) \alpha_{10}^l + n_{\bar{c}_l b_l}(\tau) \alpha_{01}^l \\
   &= \sum_l U_l \int_\tau
   \left( n_{\bar{a}_l b_l}(\tau) - \alpha_{00}^l \right)
   \left( n_{\bar{c}_l d_l}(\tau) - \alpha_{11}^l \right) +
   n_{\bar{a}_l d_l}(\tau) \alpha_{10}^l + n_{\bar{c}_l b_l}(\tau) \alpha_{01}^l +
   \mathrm{const} \\
   &= \sum_l U_l \int_\tau
   - \left( n_{\bar{a}_l d_l}(\tau) - \alpha_{01}^l \right)
     \left( n_{\bar{c}_l b_l}(\tau) - \alpha_{10}^l \right)
   - n_{\bar{a}_l b_l}(\tau) \alpha_{11}^l
   - n_{\bar{c}_l d_l}(\tau) \alpha_{00}^l
   + \mathrm{const}

This quadratic shift has to be accounted for in both the kinetic and the
interaction term of the action such that the overall action remains
unchanged, i.e.

.. math::

   \mathcal{S}=\tilde{\mathcal{S}}_0+\tilde{\mathcal{S}}_{\mathrm{int}},

with

.. math::

   \tilde{\mathcal{S}}_0=-\bar{c}_{\bar{\alpha}}\tilde{\mathcal{G}}_{\bar{\alpha}\beta}^{-1}c_\beta.

Here we introduced the non-interacting propagator by the quadratic cross
terms of the product above

.. math::

   \tilde{\mathcal{G}}_{\bar{e}f}^{-1}(i\omega)\equiv\mathcal{G}_{\bar{e}f}^{-1}(i\omega)-\sum_l U_l\left(\delta_{\bar{e},\bar{a}_l}\delta_{f,b_l}\alpha_{11}^l+\delta_{\bar{e},\bar{c}_l}\delta_{f,d_l}\alpha_{00}^l-\delta_{\bar{e},\bar{a}_l}\delta_{f,d_l}\alpha_{10}^l-\delta_{\bar{e}\bar{c}_l}\delta_{fb_l}\alpha_{01}^l\right)

The path-integral of :math:`Z` has to be adjusted accordingly

.. math::

   Z=Z_0\sum_{k=0}^\infty\frac{(-1)^k}{k!}\sum_{\mathbf{l}}\left(\prod_{i=1}^k U_{l_i}\right)\int_{\boldsymbol{\tau}}\det\hat{\mathcal{G}}_{\mathbf{l}}[\boldsymbol{\tau}],

with matrix

.. math::
   :label: det_alpha

   \hat{\mathcal{G}}_{\mathbf{l}}[\boldsymbol{\tau}]=\begin{pmatrix}\tilde{G}_{\mathbf{b_l\bar{a}_l}}[\boldsymbol{\tau}]-\operatorname{diag}(\boldsymbol{\alpha_{00}^{\mathbf{l}}}) &  & \tilde{G}_{\mathbf{b_l\bar{c}_l}}[\boldsymbol{\tau}]-\operatorname{diag}(\boldsymbol{\alpha}_{\mathbf{10}}^{\mathbf{l}})\\
   \tilde{G}_{\mathbf{d_l\bar{a}_l}}[\boldsymbol{\tau}]-\operatorname{diag}(\boldsymbol{\alpha}_{\mathbf{01}}^{\mathbf{l}}) &  & \tilde{G}_{\mathbf{d_l\bar{c}_l}}[\boldsymbol{\tau}]-\operatorname{diag}(\boldsymbol{\alpha_{11}^l})
   \end{pmatrix}

To see this, one has to rewrite an expectation value of the kind
:math:`\langle\tilde{n}_{a_{l_1}}(\tau_1)\dots\tilde{n}_{b_{l_k}}(\tau_k)\rangle`
into a determinant. This procedure is again based on the Wick theorem,
and leads to the matrix in Eq. :eq:`det_alpha`, as
can be easily proven, starting from the expression
:math:`\langle n_{a_{l_k}}(\tau_1)\dots n_{b_{l_k}}(\tau_k)\rangle`,
by successive replacements :math:`n\rightarrow\tilde{n}`. We then get

.. math::

   \Xi\equiv\sum_{k=0}^\infty\sum_{\mathbf{l}}\int_{\boldsymbol{\tau}}A_{\mathbf{l}}\times\det\hat{\mathcal{G}}_{\mathbf{l}}[\boldsymbol{\tau}]

and

.. math::

   A_{\mathbf{l}}=\frac{(-1)^k}{k!}\left(\prod_{i=1}^k U_{l_i}\right)

We solve the following set of equations self-consistently

.. math::

   \tilde{\mathcal{G}}_{b_l\bar{a}_l}[\hat{\alpha}^{\mathbf{l}}] & =\alpha_{00}^l\\
   \tilde{\mathcal{G}}_{b_l\bar{c}_l}[\hat{\alpha}^{\mathbf{l}}] & =\alpha_{10}^l\\
   \tilde{\mathcal{G}}_{d_l\bar{a}_l}[\hat{\alpha}^{\mathbf{l}}] & =\alpha_{01}^l\\
   \tilde{\mathcal{G}}_{d_l\bar{d}_l}[\hat{\alpha}^{\mathbf{l}}] & =\alpha_{11}^l

in order to determine alpha based off of the self-consistent Hartree
Fock solution. We finally introduce a small term-wise assymetry to cope
with the trivial part of the sign problem

.. math::

   \hat{\alpha}_l\rightarrow\hat{\alpha}_l+\begin{pmatrix}\delta & \delta\\
   \delta & -\delta
   \end{pmatrix}

General dynamic interactions
----------------------------

Let us now consider an interaction of the most general form

.. math::

   \mathcal{S}_{\mathrm{int}} & =\mathop{\sum\int}_{\bar{\alpha}\beta\bar{\gamma}\delta}U_{\bar{\alpha}\beta\bar{\gamma}\delta}\,\bar{c}_\alpha c_\beta \bar{c}_\gamma c_\delta\\
    & =\sum_{\bar{a}b\bar{c}d}\int_{\boldsymbol{\tau}}U_{\bar{a}b\bar{c}d}(\tau_a,\tau_b,\tau_c,\tau_d)\times\bar{c}_a(\tau_a)c_b(\tau_b)\bar{c}_c(\tau_c)c_d(\tau_d).

In this case the expansion of the partition function takes the form

.. math::

   Z=Z_0\sum_{k=0}^\infty\frac{(-1)^k}{k!}\mathop{\sum\int}_{\bar{\boldsymbol{\alpha}}\boldsymbol{\beta}\bar{\boldsymbol{\gamma}}\boldsymbol{\delta}}\left(\prod_{i=1}^k U_{\bar{\alpha}_i\beta_i\bar{\gamma}_i\delta_i}\right)\langle\bar{c}_{\alpha_1}c_{\beta_1}\bar{c}_{\gamma_1}c_{\delta_1}\ldots\bar{c}_{\alpha_k}c_{\beta_k}\bar{c}_{\gamma_k}c_{\delta_k}\rangle_0,

which again, using the Wick theorem, takes the compact form

.. math::

   Z=Z_0\sum_{k=0}^\infty\frac{(-1)^k}{k!}\mathop{\sum\int}_{\bar{\boldsymbol{\alpha}}\boldsymbol{\beta}\bar{\boldsymbol{\gamma}}\boldsymbol{\delta}}\left(\prod_{i=1}^k U_{\bar{\alpha}_i\beta_i\bar{\gamma}_i\delta_i}\right)\det\hat{\mathcal{G}}_k,

with a matrix

.. math::

   [\hat{\mathcal{G}}_k]_{ij}=\mathcal{G}_{x_i\bar{y}_j},

that now has indices

.. math::

   \boldsymbol{x}[\boldsymbol{\beta},\boldsymbol{\delta}]\equiv[{\beta_1}{\delta_1}\ldots{\beta_k}{\delta_k}],\qquad\bar{\boldsymbol{y}}[\bar{\boldsymbol{\alpha}},\bar{\boldsymbol{\gamma}}]\equiv[\bar{\alpha}_1\bar{\gamma}_1\ldots\bar{\alpha}_k\bar{\gamma}_k].

In summary we find that

.. math::
   :label: qmc_sum

   \Xi=\frac{Z}{Z_0}=\sum_{k=0}^\infty\mathop{\sum\int}_{\bar{\boldsymbol{\alpha}}\boldsymbol{\beta}\bar{\boldsymbol{\gamma}}\boldsymbol{\delta}}A_k[\boldsymbol{\alpha},\boldsymbol{\beta},\boldsymbol{\gamma},\boldsymbol{\delta}]\times\det\hat{\mathcal{G}}_k[\boldsymbol{\alpha},\boldsymbol{\beta},\boldsymbol{\gamma},\boldsymbol{\delta},\boldsymbol{s}]

with

.. math::

   A_k[\boldsymbol{\alpha},\boldsymbol{\beta},\boldsymbol{\gamma},\boldsymbol{\delta}]=\frac{(-1)^k}{k!}\left(\prod_{i=1}^k U_{\bar{\alpha}_i\beta_i\bar{\gamma}_i\delta_i}\right).

Eq. :eq:`qmc_sum` can be understood as a sum over
weights :math:`w_{\mathcal{C}}` for all configurations :math:`\mathcal{C}`

.. math::

   \sum_{\mathcal{C}}w_{\mathcal{C}}

Dynamic density-density and spin-exchange interactions
------------------------------------------------------

Let us now consider retarded interactions of both the density-density
(:math:`\mathcal{D}`) and the spin-exchange
(:math:`\mathcal{J}^\perp`) type. In this case the interacting part of
the action reads

.. math::

   \mathcal{S}_{\mathrm{int}} &= \frac{1}{N_s}\iint_{\tau,\tau'}\sum_{abs}\mathcal{D}_{ab}(\tau-\tau')\tilde{n}_{as}(\tau)\tilde{n}_{bs}(\tau')\\
   & +\frac{1}{2}\iint_{\tau,\tau'}\sum_{uv}\mathcal{J}_{uv}^\perp(\tau-\tau')\left[s_u^+(\tau)s_v^-(\tau')+s_u^-(\tau)s_v^+(\tau')\right]

with

.. math::

   s_u^+(\tau) & =\bar{c}_{u,\uparrow}(\tau)c_{u,\downarrow}(\tau),\\
   s_u^-(\tau) & =\bar{c}_{u,\downarrow}(\tau)c_{u,\uparrow}(\tau).

As done already in Sec. :ref:`static_dens`, we have to adjust
the non-interacting propagator due to the use of the shifted densities
:math:`\tilde{n}`

.. math::
   :label: Gtilde_def

   \tilde{\mathcal{G}}_{\bar{a}b}^{-1}(i\omega)\equiv\mathcal{G}_{\bar{a}b}^{-1}(i\omega)-\delta_{\bar{a},b}\sum_c\frac{1}{N_s}\left[\mathcal{D}_{bc}(i\Omega=0)+\mathcal{D}_{cb}(i\Omega=0)\right]\sum_s\alpha_{cs}.

Monte-Carlo algorithm
---------------------

The sum in Eq. :eq:`Xi_def` is evaluated using the
Metropolis-Hastings algorithm. It can be understood as a sum over
configurations
:math:`\mathcal{C}\equiv(k,\mathbf{a},\mathbf{b},\boldsymbol{\tau})`
(and possibly :math:`{\bf s}`) with different probabilities
:math:`|w_{\mathcal{C}}|`, i.e.

.. math::

   \Xi=\sum_{\mathcal{C}}w_{\mathcal{C}}=\sum_{\mathcal{C}}|w_{\mathcal{C}}|\operatorname{sign}(\mathcal{C}).

The weights :math:`|w_{\mathcal{C}}|` are used to construct a Markov chain
of configurations satisfying detailed balance and ergodicity. Under
these conditions, the expectation value of an observable :math:`\hat{O}`
can be computed as:

.. math::
   :label: MC_average_def

   \langle\hat{O}\rangle=\frac{\sum_{\mathcal{C}}|w_{\mathcal{C}}|\operatorname{sign}(\mathcal{C})\hat{O}_{\mathcal{C}}}{\sum_{\mathcal{C}}|w_{\mathcal{C}}|\operatorname{sign}(\mathcal{C})}=\frac{\langle\operatorname{sign}\hat{O}\rangle_{\mathrm{MC}}}{\langle\operatorname{sign}\rangle_{\mathrm{MC}}}

where we have defined

.. math::
   :label: MC_avg_def

   \langle\hat{A}\rangle_{\mathrm{MC}}\equiv\sum_{\mathcal{C}}|w_{\mathcal{C}}|\hat{A}_{\mathcal{C}}.

Here, :math:`\hat{A}_{\mathcal{C}}` is the value of observable
:math:`\hat{A}` when evaluated in the MC configuration
:math:`\mathcal{C}`.

Correlation functions and measurements
======================================

In this section we discuss the measurement of different correlation
functions in CTInt. Note that in the following we will not explicitly
denote that the propagator has been shifted using a
:math:`\tilde{\mathcal{G}}`, but we will rather use the symbol
:math:`\mathcal{G}`.

Single-particle Green’s function :math:`G`
------------------------------------------

From Eq. :eq:`partition_function_def` we
see that
:math:`\frac{\partial Z}{\partial\mathcal{G}_{\bar{\beta}\alpha}^{-1}}=Z\langle\bar{c}_{\bar{\beta}}c_\alpha\rangle`,
i.e.

.. math::
   :label: G_as_derivative_of_Z

   G_{\alpha\bar{\beta}}=\frac{1}{Z}\frac{\partial Z}{\partial\mathcal{G}_{\bar{\beta}\alpha}^{-1}}

Using
Eqs. :eq:`G_as_derivative_of_Z`, :eq:`Z_as_Z0Xi`
and :eq:`Xi_def` this can be further evaluated to

.. math::
   :label: G_interm

   G_{\alpha\bar{\beta}}=\frac{1}{Z}\left[\frac{\partial Z_0}{\partial\mathcal{G}_{\bar{\beta}\alpha}^{-1}}\Xi+Z_0\frac{\partial\Xi}{\partial\mathcal{G}_{\bar{\beta}\alpha}^{-1}}\right]=\frac{1}{Z}\left[\mathcal{G}_{\alpha\bar{\beta}}Z_0\Xi-Z_0\Xi\mathcal{G}_{\alpha\bar{\delta}}\frac{1}{\Xi}\frac{\partial\Xi}{\partial\mathcal{G}_{\gamma\bar{\delta}}}\mathcal{G}_{\gamma\bar{\beta}}\right].

Hence

.. math::
   :label: G_vs_M

   G_{\alpha\bar{\beta}}=\mathcal{G}_{\alpha\bar{\beta}}+\mathcal{G}_{\alpha\bar{\delta}}M_{\bar{\delta}\gamma}\mathcal{G}_{\gamma\bar{\beta}},

with

.. math::
   :label: M_def

   M_{\bar{\alpha}\beta}\equiv-\frac{1}{\Xi}\frac{\partial\Xi}{\partial\mathcal{G}_{\beta\bar{\alpha}}}.

Instead of measuring :math:`G` directly, one typically measures the
quantity :math:`M_{\bar{\alpha}\beta}` instead, as is done in the CTINT
code of TRIQS. Let us see how it can be written by means of a
Monte-Carlo average. Eq. :eq:`Xi_def` gives

.. math::

   M_{\bar{\alpha}\beta}=-\frac{1}{\Xi}\sum_{k=0}^\infty\sum_{\mathbf{a},\mathbf{b},\mathbf{s}}\int_{\boldsymbol{\tau}}A_k\frac{\partial\det\hat{\mathcal{G}}_k}{\partial\mathcal{G}_{\beta\bar{\alpha}}}.

Using the chain rule

.. math::

   \frac{\partial\det\hat{\mathcal{G}}_k}{\partial\mathcal{G}_{\beta\bar{\alpha}}}=\sum_{ij}\frac{\partial\det\hat{\mathcal{G}}_k}{\partial[\mathcal{\hat{G}}_k]_{ij}}\frac{\partial[\hat{\mathcal{G}}_k]_{ij}}{\partial\mathcal{G}_{\beta\bar{\alpha}}}=\det\hat{\mathcal{G}}_k\sum_{ij}[\hat{\mathcal{G}}_k^{-1}]_{ji}\delta_{x_i\beta}\delta_{\bar{y}_i\bar{\alpha}}

we obtain

.. math::
   :label: M_MC_final

   M_{\bar{\alpha}\beta} & =-\frac{1}{\Xi}\sum_{k=0}^\infty\sum_{\mathbf{a},\mathbf{b},\mathbf{s}}\int_{\boldsymbol{\tau}}A_k\frac{\partial\det\hat{\mathcal{G}}_k}{\partial\mathcal{G}_{\beta\bar{\alpha}}}\nonumber \\
    & =-\frac{1}{\Xi}\sum_{\mathcal{C}}w_{\mathcal{C}}\left\{ \sum_{ij}\delta_{\bar{y}_i\bar{\alpha}}\delta_{x_i\beta}[\hat{\mathcal{G}}_k^{-1}]_{ji}\right\} \nonumber \\
    & =\frac{\Big\langle-\sum_{ij}\delta_{\bar{y}_i\bar{\alpha}}\delta_{x_i\beta}[\hat{\mathcal{G}}_k^{-1}]_{ji}\,\operatorname{sign}(\mathcal{C})\Big\rangle_{\mathrm{MC}}}{\Big\langle\operatorname{sign}(\mathcal{C})\Big\rangle_{\mathrm{MC}}}.

Imaginary time measurement
^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| Just like the Green function, :math:`M_{\bar{\alpha}\beta}` depends
  only on the time-difference
  :math:`\tau=\tau_{\bar{\alpha}}-\tau_\beta`, which, in practice, is
  chosen from a discrete set of points :math:`\{\tau_l\}`. They are
  determined by binning the interval :math:`[0,\beta)` equidistantly,
  such that grid points are seperated by :math:`\Delta\tau`. We thus
  define

  .. math::

    M_{\bar{a}b}^{\Delta\tau}(\tau_l)=-\frac{1}{\beta\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\left\langle \sum_{ij}\delta_{\bar{c}_j\bar{a}}\delta_{c_i b}\left[\delta_{\Delta\tau}(\tau_l-\tau_j+\tau_i)-\delta_{\Delta\tau}(\tau_l-\beta-\tau_j+\tau_i)\right][\hat{\mathcal{G}}_k^{-1}]_{ji}\operatorname{sign}(\mathcal{C})\right\rangle _{\mathrm{MC}}

  where :math:`\delta_{\Delta\tau}` is the broadened delta function

  .. math::

     \delta_{\Delta\tau}(\tau)\equiv\begin{cases}
     \frac{1}{\Delta\tau} & 0<\tau<\Delta\tau\\
     0 & \tau\leq0\quad\mathrm{or}\quad\tau\geq\Delta\tau,
     \end{cases}

  and a factor :math:`1/\beta` should be introduced when making use of
  time-translational invariance
  :math:`(\tau_{\bar{\alpha}},\tau_\beta)\rightarrow\tau=\tau_{\bar{\alpha}}-\tau_\beta`.
  The second :math:`\delta` function with the minus sign takes care of
  shifting negative time-differences
  :math:`\tau_{\bar{\alpha}}-\tau_\beta` to the interval
  :math:`[0,\beta)` using the antiperiodicity of :math:`M`. As
  :math:`M_{\bar{u}v}(\tau_l)=\lim_{\Delta\tau\rightarrow0}{M}_{\bar{u}v}^{\Delta\tau}(\tau_l)`
  we can approximate
  :math:`M_{\bar{u}v}(\tau_l)\approx M_{\bar{u}v}^{\Delta\tau}(\tau_l)`
  for a sufficiently small :math:`\Delta\tau`.

Frequency measurement
^^^^^^^^^^^^^^^^^^^^^

|  
| Alternatively, we can measure :math:`M` directly in frequency space,
  i.e.

  .. math::

     M_{\bar{a}b}(i\omega_n)=-\frac{1}{\beta\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\left\langle \sum_{ij}\delta_{\bar{c}_j\bar{a}}\delta_{c_i b}e^{i\omega_n(\tau_j-\tau_i)}[\hat{\mathcal{G}}_k^{-1}]_{ji}\operatorname{sign}(\mathcal{C})\right\rangle _{\mathrm{MC}}

  Here, we can avoid the binning, while a factor :math:`1/\beta` has to
  be introduced as before. To evaluate this measure, we can make use of
  the NFFT library. To make this more obvious we rewrite

  .. math::

     e^{i\omega_n(\tau_j-\tau_i)}=e^{i\pi(2n+1)\left(\frac{\tau_j-\tau_i}{\beta}-\frac{1}{2}+\frac{1}{2}\right)}=e^{2\pi inx}e^{i\pi\left(n+x+\frac{1}{2}\right)}

  with
  :math:`x\equiv\left(\frac{\tau_j-\tau_i}{\beta}-\frac{1}{2}\right)\in[-0.5,0.5)`
  and :math:`n\in\mathbb{Z}`, as is required by the interface of the
  library.

Two-particle Green’s function :math:`\chi^4`
--------------------------------------------

The same idea can be pursued when considering two-particle Green
functions. From
Eq. :eq:`partition_function_def` we
can directly see that
:math:`\frac{\partial Z}{\partial\mathcal{G}_{\bar{\alpha}\beta}^{-1}\partial\mathcal{G}_{\bar{\gamma}\delta}^{-1}}=Z\langle\bar{c}_\alpha c_\beta \bar{c}_\gamma c_\delta\rangle`
i.e

.. math::
   :label: chi4_as_derivative_of_Z

   \chi_{\bar{\alpha}\beta\bar{\gamma}\delta}^4=\langle\bar{c}_\alpha c_\beta \bar{c}_\gamma c_\delta\rangle=\frac{1}{Z}\frac{\partial Z}{\partial\mathcal{G}_{\bar{\alpha}\beta}^{-1}\partial\mathcal{G}_{\bar{\gamma}\delta}^{-1}}

Starting from Eqs. :eq:`chi4_as_derivative_of_Z` and :eq:`G_as_derivative_of_Z`
we first note:

.. math::
   :label: chi4_interm

   \chi_{\bar{\alpha}\beta\bar{\eta}\kappa}^4=\frac{1}{Z}\frac{\partial(ZG_{\beta\bar{\alpha}})}{\partial\mathcal{G}_{\bar{\eta}\kappa}^{-1}}=-\frac{1}{Z}\mathcal{G}_{\mu\bar{\eta}}\frac{\partial(ZG_{\beta\bar{\alpha}})}{\partial\mathcal{G}_{\mu\bar{\nu}}}\mathcal{G}_{\kappa\bar{\nu}}

Let us now evaluate
:math:`\frac{\partial(ZG_{\beta\bar{\alpha}})}{\partial\mathcal{G}_{\mu\bar{\nu}}}`
using :eq:`G_interm`. We have

.. math::
   :label: chi4_interm_2

   \frac{\partial(ZG_{\beta\bar{\alpha}})}{\partial\mathcal{G}_{\mu\bar{\nu}}} & =\frac{\partial}{\partial\mathcal{G}_{\mu\bar{\nu}}}\left[\mathcal{G}_{\beta\bar{\alpha}}Z_0\Xi-Z_0\mathcal{G}_{\gamma\bar{\alpha}}\frac{\partial\Xi}{\partial\mathcal{G}_{\gamma\bar{\delta}}}\mathcal{G}_{\beta\bar{\delta}}\right]\nonumber \\
    & =\delta_{\beta\mu}\delta_{\bar{\nu}\bar{\alpha}}Z_0\Xi+\mathcal{G}_{\beta\bar{\alpha}}\frac{\partial Z_0}{\partial\mathcal{G}_{\mu\bar{\nu}}}\Xi+\mathcal{G}_{\beta\bar{\alpha}}Z_0\frac{\partial\Xi}{\partial\mathcal{G}_{\mu\bar{\nu}}}\nonumber \\
    & \;\;-\frac{\partial Z_0}{\partial\mathcal{G}_{\mu\bar{\nu}}}\mathcal{G}_{\gamma\bar{\alpha}}\frac{\partial\Xi}{\partial\mathcal{G}_{\gamma\bar{\delta}}}\mathcal{G}_{\beta\bar{\delta}}-Z_0\delta_{\mu\gamma}\delta_{\bar{\alpha}\bar{\nu}}\frac{\partial\Xi}{\partial\mathcal{G}_{\gamma\bar{\delta}}}\mathcal{G}_{\beta\bar{\delta}}-Z_0\mathcal{G}_{\gamma\bar{\alpha}}\frac{\partial\Xi}{\partial\mathcal{G}_{\mu\bar{\nu}}\partial\mathcal{G}_{\gamma\bar{\delta}}}\mathcal{G}_{\beta\bar{\delta}}-Z_0\mathcal{G}_{\gamma\bar{\alpha}}\frac{\partial\Xi}{\partial\mathcal{G}_{\gamma\bar{\delta}}}\delta_{\beta\mu}\delta_{\bar{\delta}\bar{\nu}}\nonumber \\
    & =\delta_{\beta\mu}\delta_{\bar{\nu}\bar{\alpha}}Z_0\Xi-\mathcal{G}_{\beta\bar{\alpha}}\mathcal{G}_{\bar{\nu}\mu}^{-1}Z_0\Xi-\mathcal{G}_{\beta\bar{\alpha}}Z_0\Xi M_{\bar{\nu}\mu}\nonumber \\
    & \;\;-\mathcal{G}_{\bar{\nu}\mu}^{-1}Z_0\mathcal{G}_{\gamma\bar{\alpha}}\Xi M_{\bar{\delta}\gamma}\mathcal{G}_{\beta\bar{\delta}}+\Xi Z_0\delta_{\mu\gamma}\delta_{\bar{\alpha}\bar{\nu}}M_{\bar{\delta}\gamma}\mathcal{G}_{\beta\bar{\delta}}-Z_0\mathcal{G}_{\gamma\bar{\alpha}}\Xi M_{\bar{\nu}\mu\bar{\delta}\gamma}^4\mathcal{G}_{\beta\bar{\delta}}+Z_0\mathcal{G}_{\gamma\bar{\alpha}}\Xi M_{\bar{\delta}\gamma}\delta_{\beta\mu}\delta_{\bar{\delta}\bar{\nu}}\nonumber \\
    & =Z\Big\{\delta_{\beta\mu}\delta_{\bar{\nu}\bar{\alpha}}-\mathcal{G}_{\beta\bar{\alpha}}\mathcal{G}_{\bar{\nu}\mu}^{-1}-\mathcal{G}_{\beta\bar{\alpha}}M_{\bar{\nu}\mu}-\mathcal{G}_{\bar{\nu}\mu}^{-1}\mathcal{G}_{\gamma\bar{\alpha}}M_{\bar{\delta}\gamma}\mathcal{G}_{\beta\bar{\delta}}\nonumber \\
    & \;\;\;+\delta_{\bar{\alpha}\bar{\nu}}M_{\bar{\delta}\mu}\mathcal{G}_{\beta\bar{\delta}}-\mathcal{G}_{\gamma\bar{\alpha}}M_{\mu\bar{\nu}\gamma\bar{\delta}}^4\mathcal{G}_{\beta\bar{\delta}}+\mathcal{G}_{\gamma\bar{\alpha}}M_{\bar{\nu}\gamma}\delta_{\beta\mu}\Big\}

where we have defined

.. math::
   :label: M4_def

   M_{\mu\bar{\nu}\gamma\bar{\delta}}^4\equiv\frac{1}{\Xi}\frac{\partial\Xi}{\partial\mathcal{G}_{\mu\bar{\nu}}\partial\mathcal{G}_{\gamma\bar{\delta}}}.

Plugging :eq:`chi4_interm_2` into
:eq:`chi4_interm`, we obtain

.. math::

   \chi_{\bar{\alpha}\beta\bar{\eta}\kappa}^4 & =-\mathcal{G}_{\beta\bar{\eta}}\mathcal{G}_{\kappa\bar{\alpha}}+\mathcal{G}_{\kappa\bar{\eta}}\mathcal{G}_{\beta\bar{\alpha}}+\mathcal{G}_{\beta\bar{\alpha}}\mathcal{G}_{\kappa\bar{\nu}}M_{\bar{\nu}\mu}\mathcal{G}_{\mu\bar{\eta}}+\mathcal{G}_{\beta\bar{\delta}}M_{\bar{\delta}\gamma}\mathcal{G}_{\gamma\bar{\alpha}}\mathcal{G}_{\kappa\bar{\eta}}\\
    & \;\;-\mathcal{G}_{\beta\bar{\delta}}M_{\bar{\delta}\mu}\mathcal{G}_{\mu\bar{\eta}}\mathcal{G}_{\kappa\bar{\alpha}}+\mathcal{G}_{\mu\bar{\alpha}}\mathcal{G}_{\beta\bar{\nu}}M_{\mu\bar{\nu}\gamma\bar{\delta}}^4\mathcal{G}_{\gamma\bar{\eta}}\mathcal{G}_{\kappa\bar{\delta}}-\mathcal{G}_{\beta\bar{\eta}}\mathcal{G}_{\kappa\bar{\nu}}M_{\bar{\nu}\gamma}\mathcal{G}_{\gamma\bar{\alpha}}

i.e, after relabeling and reordering,

.. math::
   :label: chi4_final

   \chi_{\bar{\alpha}\beta\bar{\gamma}\delta}^4 & =\mathcal{G}_{\mu\bar{\alpha}}\mathcal{G}_{\beta\bar{\nu}}M_{\mu\bar{\nu}\kappa\bar{\eta}}^4\mathcal{G}_{\kappa\bar{\gamma}}\mathcal{G}_{\delta\bar{\eta}}\nonumber \\
    & \;\;+\mathcal{G}_{\beta\bar{\alpha}}\mathcal{G}_{\delta\bar{\gamma}}+[\mathcal{G}M\mathcal{G}]_{\beta\bar{\alpha}}\mathcal{G}_{\delta\bar{\gamma}}+\mathcal{G}_{\beta\bar{\alpha}}[\mathcal{G}M\mathcal{G}]_{\delta\bar{\gamma}}\nonumber \\
    & \;\;-\mathcal{G}_{\beta\bar{\gamma}}\mathcal{G}_{\delta\bar{\alpha}}-[\mathcal{G}M\mathcal{G}]_{\beta\bar{\gamma}}\mathcal{G}_{\delta\bar{\alpha}}-\mathcal{G}_{\beta\bar{\gamma}}[\mathcal{G}M\mathcal{G}]_{\delta\bar{\alpha}}.

Thus, in order to compute :math:`\chi^4`, one needs to compute
:math:`M` and :math:`M^4`.

Connected four-point function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| By inspection of Eq.:eq:`chi4_final` we see that
  the connected part of the two-particle Green function now reads

  .. math::
     :label: chi4tilde_conn_final

     \chi_{\bar{\alpha}\beta\bar{\gamma}\delta}^{4,\mathrm{conn}}=\mathcal{G}_{\mu\bar{\alpha}}\mathcal{G}_{\beta\bar{\nu}}M_{\mu\bar{\nu}\kappa\bar{\eta}}^4\mathcal{G}_{\kappa\bar{\gamma}}\mathcal{G}_{\delta\bar{\eta}}-\left[\mathcal{G}M\mathcal{G}\right]_{\alpha\bar{\beta}}\left[\mathcal{G}M\mathcal{G}\right]_{\gamma\bar{\delta}}+\left[\mathcal{G}M\mathcal{G}\right]_{\alpha\bar{\delta}}\left[\mathcal{G}M\mathcal{G}\right]_{\gamma\bar{\beta}}.

  Instead, we can write

  .. math::

    \chi_{\bar{\alpha}\beta\bar{\gamma}\delta}^{4,\mathrm{conn}}=\mathcal{G}_{\mu\bar{\alpha}}\mathcal{G}_{\beta\bar{\nu}}M_{\mu\bar{\nu}\kappa\bar{\eta}}^{4,\mathrm{conn}}\mathcal{G}_{\kappa\bar{\gamma}}\mathcal{G}_{\delta\bar{\eta}}

  with

  .. math::

    M_{\mu\bar{\nu}\kappa\bar{\eta}}^{4,\mathrm{conn}}\equiv M_{\mu\bar{\nu}\kappa\bar{\eta}}^4-M_{\bar{\nu}\mu}M_{\bar{\eta}\kappa}+M_{\bar{\nu}\kappa}M_{\bar{\eta}\mu}.

Fully reducible vertex
^^^^^^^^^^^^^^^^^^^^^^

|  
| Using the definition of the fully reducible vertex function :math:`F`
  (see vertex conventions), we get

  .. math::

    F_{\alpha\bar{\beta}\gamma\bar{\delta}}=\bigl\{ G_{\bar{\beta}\beta}^{-1}\mathcal{G}_{\beta\bar{\nu}}\bigr\}\left\{ G_{\bar{\delta}\delta}^{-1}\mathcal{G}_{\delta\bar{\eta}}\right\} M_{\mu\bar{\nu}\kappa\bar{\eta}}^{4,\mathrm{conn}}\left\{ \mathcal{G}_{\mu\bar{\alpha}}G_{\bar{\alpha}\alpha}^{-1}\right\} \left\{ \mathcal{G}_{\kappa\bar{\gamma}}G_{\bar{\gamma}\gamma}^{-1}\right\} .

  With Eq. :eq:`G_vs_M` we can calculate :math:`F`
  directly from :math:`M^4` and :math:`M`

  .. math::

    F_{\alpha\bar{\beta}\gamma\bar{\delta}}=\left\{ 1+M\mathcal{G}\right\} _{\bar{\beta}\bar{\nu}}^{-1}\left\{ 1+M\mathcal{G}\right\} _{\bar{\delta}\bar{\eta}}^{-1}M_{\mu\bar{\nu}\kappa\bar{\eta}}^{4,\mathrm{conn}}\left\{ 1+M\mathcal{G}\right\} _{\mu\alpha}^{-1}\left\{ 1+M\mathcal{G}\right\} _{\kappa\gamma}^{-1}.

  From this expression, one can see that :math:`F` and
  :math:`M^{4,\mathrm{conn}}` have the same asymptotics.

Imaginary time measurement of :math:`M^4`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| In order to write :math:`M^4` by means of a Monte-Carlo average we
  have to evaluate

  .. math::

    M_{\alpha\bar{\beta}\gamma\bar{\delta}}^4=\frac{1}{\Xi}\frac{\partial\Xi}{\partial\mathcal{G}_{\alpha\bar{\beta}}\partial\mathcal{G}_{\gamma\bar{\delta}}}=-\frac{1}{\Xi}\sum_{k=0}^\infty\sum_{\mathbf{a},\mathbf{b},\mathbf{s}}\int_{\boldsymbol{\tau}}A_k\,\frac{\partial\det\hat{\mathcal{G}}_k}{\partial\mathcal{G}_{\alpha\bar{\beta}}\partial\mathcal{G}_{\gamma\bar{\delta}}}.

  This derivative evaluates to (neglecting the :math:`k` index for now)

  .. math::

     \frac{\partial\det\hat{\mathcal{G}}}{\partial\mathcal{G}_{\alpha\bar{\beta}}\partial\mathcal{G}_{\gamma\bar{\delta}}} & =\frac{\partial}{\partial\mathcal{G}_{\alpha\bar{\beta}}}\det\hat{\mathcal{G}}\,\sum_{kl}[\hat{\mathcal{G}}^{-1}]_{lk}\delta_{x_k\gamma}\delta_{\bar{y}_l\bar{\delta}}\\
      & =\det\hat{\mathcal{G}}\,\sum_{ijkl}\left\{ [\hat{\mathcal{G}}^{-1}]_{ji}[\hat{\mathcal{G}}^{-1}]_{lk}-[\hat{\mathcal{G}}^{-1}]_{li}[\hat{\mathcal{G}}^{-1}]_{jk}\right\} \delta_{x_i\alpha}\delta_{\bar{y}_i\bar{\beta}}\delta_{x_k\gamma}\delta_{\bar{y}_l\bar{\delta}},

  where we had to use that

  .. math::

     \phantom{=}\frac{\partial}{\partial\mathcal{G}_{\alpha\bar{\beta}}}\sum_{kl}[\hat{\mathcal{G}}^{-1}]_{lk}\delta_{x_k\gamma}\delta_{\bar{y}_l\bar{\delta}}=-\sum_{ijkl}[\hat{\mathcal{G}}^{-1}]_{li}\underbrace{\frac{\partial\hat{\mathcal{G}}_{ij}}{\partial\mathcal{G}_{\alpha\bar{\beta}}}}_{\delta_{x_i\alpha}\delta_{\bar{y}_i\bar{\beta}}}[\hat{\mathcal{G}}^{-1}]_{jk}\delta_{x_k\gamma}\delta_{\bar{y}_l\bar{\delta}}.

  The Monte-Carlo average for :math:`M^4` thus takes the form

  .. math::
     :label: Measure_M4

     M_{\alpha\bar{\beta}\gamma\bar{\delta}}^4=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\sum_{ijkl}\left\{ [\hat{\mathcal{G}}^{-1}]_{ji}[\hat{\mathcal{G}}^{-1}]_{lk}-[\hat{\mathcal{G}}^{-1}]_{li}[\hat{\mathcal{G}}^{-1}]_{jk}\right\} \delta_{x_i\alpha}\delta_{\bar{y}_j\bar{\beta}}\delta_{x_k\gamma}\delta_{\bar{y}_l\bar{\delta}}\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}.

  Using the intermediate scattering matrix

  .. math::

    \bar{M}_{\bar{\beta}\alpha}\equiv\sum_{ij}[\hat{\mathcal{G}}^{-1}]_{ji}\delta_{x_i\alpha}\delta_{\bar{y}_j\bar{\beta}}

  it takes the compact form

  .. math::

     M_{\alpha\bar{\beta}\gamma\bar{\delta}}^4=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Big\langle\left(\bar{M}_{\bar{\beta}\alpha}\bar{M}_{\bar{\delta}\gamma}-\bar{M}_{\bar{\delta}\alpha}\bar{M}_{\bar{\beta}\gamma}\right)\operatorname{sign}(\mathcal{C})\Big\rangle_{\mathrm{MC}}.

  This factorization is directly used in our implementation.

Frequency measurement of :math:`M^4`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| Just like :math:`M`, we can measure :math:`M^4` also directly in
  Matsubara frequencies [1]_

  .. math::

     \begin{split}M_{a\bar{b}c\bar{d}}^4(i\omega_a,i\omega_{\bar{b}},i\omega_c) & =\frac{1}{\beta\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\left(\bar{M}_{\bar{\beta}\alpha}\bar{M}_{\bar{\delta}\gamma}-\bar{M}_{\bar{\delta}\alpha}\bar{M}_{\bar{\beta}\gamma}\right)\\
      & \hspace{2cm}\times e^{-i\omega_a(\tau_a-\tau_d)}\,e^{i\omega_{\bar{b}}(\tau_b-\tau_d)}\,e^{-i\omega_c(\tau_c-\tau_d)}\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}\\
      & =\frac{1}{\beta\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\bar{M}_{\bar{b}a}(\omega_b,\omega_a)\bar{M}_{\bar{d}c}(\omega_a+\omega_c-\omega_b,\omega_c)\\
      & \hspace{2cm}-\bar{M}_{\bar{d}a}(\underbrace{\omega_a+\omega_c-\omega_b}_{\omega_d},\omega_a)\bar{M}_{\bar{b}c}(\omega_b,\omega_c)\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}.
     \end{split}

  Here we have defined the intermediate scattering matrix in the
  Matsubara domain

  .. math::

    \bar{M}_{\bar{b}a}(\omega_b,\omega_a)\equiv\sum_{ij}[\hat{\mathcal{G}}^{-1}]_{ji}\delta_{c_i a}\delta_{\bar{c}_j\bar{b}}\times\,e^{i\omega_b\tau_j}e^{-i\omega_a\tau_i}.

  In the implementation we precompute this intermediate scattering
  matrix, which allows us to reduce the effort from
  :math:`\mathcal{O}(k^4\log(k))` down to
  :math:`\mathcal{O}(k^2\log(k))`.

Equal-time correlation functions :math:`\chi^3`
-----------------------------------------------

The equal-time correlation functions :math:`\chi^{3,r}` can be directly
calculated from the intermediate quantity :math:`M^{3,r}` defined in the
following. In the first step we define three extensions of the
intermediate scattering matrix

.. math::

   \overline{\mathcal{G}M}_{b\alpha}\equiv\overline{\mathcal{G}M}_{ba}(\tau_a)\equiv\sum_{ij}\mathcal{G}_{b\bar{c}_j}(-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\delta_{x_i\alpha}=-\sum_{ij}\mathcal{G}_{b\bar{c}_j}(\beta-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\delta_{x_i\alpha},

.. math::

   \overline{M\mathcal{G}}_{\bar{\beta}\bar{a}}\equiv\overline{M\mathcal{G}}_{\bar{b}\bar{a}}(\tau_b)\equiv\sum_{ij}\delta_{y_j\beta}[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i),

.. math::

   \overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}\equiv\sum_{ij}\mathcal{G}_{b\bar{c}_j}(-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i)=-\sum_{ij}\mathcal{G}_{b\bar{c}_j}(\beta-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i).

Using these, :math:`M^{3,r}` is measured as

.. math::

   M_{abcd}^{3,pp}(\tau_a,\tau_c)=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\overline{\mathcal{G}M}_{ba}(\tau_a)\overline{\mathcal{G}M}_{dc}(\tau_c)-\overline{\mathcal{G}M}_{da}(\tau_a)\overline{\mathcal{G}M}_{bc}(\tau_c)\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}},

.. math::

   M_{a\bar{b}\bar{c}d}^{3,ph}(\tau_a,\tau_{\bar{b}})=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\bar{M}_{\bar{\beta}\alpha}\overline{\mathcal{G}M\mathcal{G}}_{d\bar{c}}-\overline{\mathcal{G}M}_{d\alpha}\overline{M\mathcal{G}}_{\bar{\beta}\bar{c}}\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}.

Frequency measurement of :math:`M^{3,r}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| Again, we can measure :math:`M^{3,r}` also directly in Matsubara
  frequencies. For this we need the scattering matrices

  .. math::

    \overline{\mathcal{G}M}_{ba}(\omega_a)\equiv\sum_{ij}\mathcal{G}_{b\bar{c}_j}(-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}e^{-i\omega_a\tau_i}\delta_{c_i a}=\sum_{ij}\mathcal{G}_{b\bar{c}_j}(\beta-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}e^{i\omega_a(\beta-\tau_i)}\delta_{c_i a},

  .. math::

    \overline{M\mathcal{G}}_{\bar{b}\bar{a}}(\omega_{\bar{b}})\equiv\sum_{ij}e^{i\omega_{\bar{b}}\tau_j}\delta_{\bar{c}_j\bar{b}}[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i),

  .. math::

    \overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}\equiv\sum_{ij}\mathcal{G}_{b\bar{c}_j}(-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i)=-\sum_{ij}\mathcal{G}_{b\bar{c}_j}(\beta-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i),

Then we can measure

.. math::

   M_{abcd}^{3,pp}(\omega_a,\omega_c)=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\overline{\mathcal{G}M}_{ba}(\omega_a)\overline{\mathcal{G}M}_{dc}(\omega_c)-\overline{\mathcal{G}M}_{da}(\omega_a)\overline{\mathcal{G}M}_{bc}(\omega_c)\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}},

.. math::

   M_{a\bar{b}\bar{c}d}^{3,ph}(\omega_a,\omega_{\bar{b}})=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\bar{M}_{\bar{b}a}(\omega_{\bar{b}},\omega_a)\overline{\mathcal{G}M\mathcal{G}}_{c\bar{d}}-\overline{\mathcal{G}M}_{ca}(\omega_a)\overline{M\mathcal{G}}_{\bar{b}\bar{d}}(\omega_{\bar{b}})\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}.

Constructing :math:`\chi^{3,r}` from :math:`M^{3,r}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| With

  .. math::

     M_{abcd}^{3,pp,conn}(\omega_a,\omega_c)=M_{abcd}^{3,pp}(\omega_a,\omega_c)-\overline{\mathcal{G}M}_{ba}(\omega_a)\overline{\mathcal{G}M}_{dc}(\omega_c)+\overline{\mathcal{G}M}_{da}(\omega_a)\overline{\mathcal{G}M}_{bc}(\omega_c),

  .. math::

     M_{a\bar{b}\bar{c}d}^{3,ph,conn}(\omega_a,\omega_{\bar{b}})=M_{a\bar{b}\bar{c}d}^{3,ph}(\omega_a,\omega_{\bar{b}})-\beta\delta_{\omega_a,\omega_{\bar{b}}}\bar{M}_{\bar{b}a}(\omega_a)\overline{\mathcal{G}M\mathcal{G}}_{c\bar{d}}+\overline{\mathcal{G}M}_{ca}(\omega_a)\overline{M\mathcal{G}}_{\bar{b}\bar{d}}(\omega_{\bar{b}}),

  we have that

  .. math::

     \chi_{\bar{a}b\bar{c}d}^{3,pp,conn}(\omega_a,\omega_c)=\mathcal{G}_{\bar{a}i}(\omega_a)\mathcal{G}_{\bar{c}j}(\omega_c)M_{ibjd}^{3,pp,conn}(\omega_a,\omega_c),

  .. math::

     \chi_{\bar{a}b\bar{c}d}^{3,ph,conn}(\omega_{\bar{a}},\omega_b)=\mathcal{G}_{i\bar{a}}(\omega_{\bar{a}})\mathcal{G}_{b\bar{j}}(\omega_b)M_{i\bar{j}\bar{c}d}^{3,ph,conn}(\omega_{\bar{a}},\omega_b).

  Note that :math:`M_{abcd}^{3,pp}(\omega_a,\omega_c)` and
  :math:`M_{a\bar{b}\bar{c}d}^{3,ph}(\omega_a,\omega_{\bar{b}})` can
  be either measured directly in frequencies, or calculated via a
  Fourier transforms

.. math::

   M_{abcd}^{3,pp}(\omega_a,\omega_c)=\int_{\tau_a,\tau_c}e^{-i\tau_a\omega_a}e^{-i\tau_c\omega_c}M_{abcd}^{3,pp}(\tau_a,\tau_c)

.. math::

   M_{a\bar{b}\bar{c}d}^{3,ph}(\omega_a,\omega_{\bar{b}})=\int_{\tau_a,\tau_{\bar{b}}}e^{-i\tau_a\omega_a}e^{i\tau_{\bar{b}}\omega_{\bar{b}}}M_{a\bar{b}\bar{c}d}^{3,ph}(\tau_a,\tau_{\bar{b}})

Mixed fermion-boson notation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of using the frequencies of the fermionic operators, it is often
convenient to work in a mixed fermion-boson notation. It is defined as

.. math::

   \tilde{M}_{abcd}^{3,pp}(\omega,\Omega)=M_{abcd}^{3,pp}(\omega,\Omega-\omega)
    & =\int_{\tau_a,\tau_c}e^{i(\tau_c-\tau_a)\omega}e^{-i\tau_c\Omega}M_{abcd}^{3,pp}(\tau_a,\tau_c)\\
    & =\int_{-\beta}^0d\tau'\int_0^\beta d\tau_a e^{-i(\tau'+\tau_a)\omega}e^{i\tau'\Omega}M_{abcd}^{3,pp}(\tau_a,-\tau')\\
    & =-\int_0^\beta d\tau'\int_0^\beta d\tau_a e^{-i(\tau'+\tau_a)\omega}e^{i\tau'\Omega}M_{abcd}^{3,pp}(\tau_a,\beta-\tau')\\
    & =\int_0^\beta d\tau'\int_{-\tau'}^{-\beta-\tau'}d\tau e^{i\tau\omega}e^{i\tau'\Omega}M_{abcd}^{3,pp}(-\tau-\tau',\beta-\tau')\\
    & =\int_0^\beta d\tau'\int_{-\tau'}^{\beta-\tau'}d\tau e^{i\tau\omega}e^{i\tau'\Omega}M_{abcd}^{3,pp}(\beta-\tau-\tau',\beta-\tau')\\

.. math::

   \tilde{M}_{a\bar{b}\bar{c}d}^{3,ph}(\omega,\Omega)=M_{a\bar{b}\bar{c}d}^{3,ph}(\omega,\Omega+\omega)=\int_{\tau_a,\tau_{\bar{b}}}e^{i(\tau_{\bar{b}}-\tau_a)\omega}e^{i\tau_{\bar{b}}\Omega}M_{a\bar{b}\bar{c}d}^{3,ph}(\tau_a,\tau_{\bar{b}}),

and accordingly for :math:`\chi^{3,r}`.

Equal-time correlation functions :math:`\chi^2`
-----------------------------------------------

We now consider the case that two pairs of times are equal. In this case
we measure in imaginary time, and a Fourier transform is performed in a
post-processing step. We require again an intermediate object

.. math::

   \overline{\mathcal{G}M\mathcal{G}}_{\beta\bar{\alpha}}\equiv\overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}(\tau_b,\tau_a)\equiv\sum_{ij}\mathcal{G}_{b\bar{c}_j}(\tau_b-\tau_j)[\hat{\mathcal{G}}^{-1}]_{ji}\mathcal{G}_{c_i\bar{a}}(\tau_i-\tau_{\bar{a}}).

The measurements then read

.. math::

   M_{\bar{a}b\bar{c}d}^{2,pp}(\tau)=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}(0^+,\tau^+)\overline{\mathcal{G}M\mathcal{G}}_{d\bar{c}}(0,\tau)-\overline{\mathcal{G}M\mathcal{G}}_{d\bar{a}}(0,\tau^+)\overline{\mathcal{G}M\mathcal{G}}_{b\bar{c}}(0^+,\tau)\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}},

.. math::

   M_{\bar{a}b\bar{c}d}^{2,ph}(\tau)=\frac{1}{\langle\operatorname{sign}(\mathcal{C})\rangle_{\mathrm{MC}}}\Bigg\langle\Big[\overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}(\tau,\tau^+)\overline{\mathcal{G}M\mathcal{G}}_{d\bar{c}}(0,0^+)-\overline{\mathcal{G}M\mathcal{G}}_{d\bar{a}}(0,\tau^+)\overline{\mathcal{G}M\mathcal{G}}_{b\bar{c}}(\tau,0^+)\Big]\operatorname{sign}(\mathcal{C})\Bigg\rangle_{\mathrm{MC}}.

Constructing :math:`\chi^{2,r}` from :math:`M^{2,r}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|  
| We have that

  .. math::

     \chi_{\bar{a}b\bar{c}d}^{2,pp,conn}(\tau)=M_{\bar{a}b\bar{c}d}^{2,pp}(\tau)-\overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}(\beta-\tau)\overline{\mathcal{G}M\mathcal{G}}_{d\bar{c}}(\beta-\tau)+\overline{\mathcal{G}M\mathcal{G}}_{d\bar{a}}(\beta-\tau)\overline{\mathcal{G}M\mathcal{G}}_{b\bar{c}}(\beta-\tau),

  .. math::

     \chi_{\bar{a}b\bar{c}d}^{2,ph,conn}(\tau)=M_{\bar{a}b\bar{c}d}^{2,ph}(\tau)-\overline{\mathcal{G}M\mathcal{G}}_{b\bar{a}}(\beta^-)\overline{\mathcal{G}M\mathcal{G}}_{d\bar{c}}(\beta^-)-\overline{\mathcal{G}M\mathcal{G}}_{d\bar{a}}(\beta-\tau)\overline{\mathcal{G}M\mathcal{G}}_{b\bar{c}}(\tau).

.. [1]
   Careful: The Fourier conventions for :math:`M`-objects are indeed
   reversed, meaning that barred indices are transformed to frequencies
   with :math:`\int_\tau e^{iw\tau}`.
