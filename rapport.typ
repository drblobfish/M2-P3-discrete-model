#import "@preview/rubber-article:0.5.2": *
#import "@preview/lovelace:0.3.1": *

#show: article.with(
  cols: none,
  eq-chapterwise: true,
  eq-numbering: "(1.1)",
  lang: "en",
  page-margins: 1in,
  page-paper: "us-letter",
)
#show table.cell.where(y: 0): strong
#set table(stroke : (x,y) => if (y==0) {(bottom : 0.7pt + black)})

// Frontmatter
#maketitle(
  title: [Discrete Models in Biology Homework

Modelling invasive Glioma Cells with a Cellular Potts Model],
  authors: ("Jules Herrmann",),
  date: datetime.today().display("[day]. [month repr:long] [year]"),
)

= Introduction

Glioblastoma is a type of cancer originating from glial cells in the brain that is very invasive.
In their #cite(<fayzullin_time-lapse_2016>, form :"year") article, #cite(<fayzullin_time-lapse_2016>, form : "author"), used time-lapse microscopy technics to study the behavior of invasive glioma cells originating from a tumor grafted on rodent brain slice.
They classified the gliomal cell in two class of phenotype. The non invasive cells stayed in the tumorsphere, while invasive cells migrated away from the bulk of the tumor.
They observed very specific migratory movement, with some type of invasive cell displaying straight radial trajectories away from the core, even after multiple day of invasion, and up to a millimeter of distance between the cells and the core.
Upon complete removal of the tumor core, the invasive cells stopped their centrifugal trajectories and displayed random movement.
However, after incomplete removal of the core, some invasive cell reversed their direction and migrated back to the tumor and participated to the reconstruction of the core.

As hypothesized in the article, these two experiment suggest some sort of signalling between the tumorsphere and the migrating invasive cells. In this project, I first present a potential mechanism that could explain the behavior observed in <fayzullin_time-lapse_2016>.

To test if this mechanism is sufficient to explain empirical data, an _in
silico_ model must be built. I choose to start from a Cellular Potts model for
three reasons. First, they are quite general, have a natural way to represent
cell/cell adhesion surface, and are known to couple well with diffusion-reaction
PDE models. Second, they are 
pedagogically interersting as I've never used them, contrary to agent-based models.
Finally, there are record in the scientific literature of previous use of this
type of model to study the behavior of glioblastoma.

In their article from
#cite(<rubenstein_role_2008>,form:"year"),#cite(<rubenstein_role_2008>,form:"author"),
tested the impact of the extra-cellular matrix on the tumor using a cellular
potts model with a hamiltonian consisting of a term for adhesion potentials and
an elastic potential driving the cell toward a target perimeter. Their cells
were also able to divide upon reach of their target perimeter. Overall, this
model focused on the evolution of the tumor core and didn't account for the
invasive cell phenotype.

In #cite(<szabo_cellular_2013>,form:"year"), #cite(<szabo_cellular_2013>,form:"author") use cellular potts model to study invasive tumor. This time they consider adhesion potentials and an elastic potential driving the area of the cells toward a target area. They extended the model to account for chemotaxis. The diffusion of chemical signals is computed by solving a diffusion equation, and, at each step of the Metropolis-Hasting algorithm, the difference in potential energy is biased proportionaly to the difference of signal concentration in the source and destination points.

I then present in detail the model that was implemented, which was simplified because of time constraints.

= A possible mechanism

Here I describe in biological term, a possible mechanism susceptible to explain the behavior observed by #cite(<fayzullin_time-lapse_2016>,form:"author").

Glioma cells can express two different phenotype, invasive (I) and non-invasive (N), and communicate using two chemical signals that are diffused in the medium: the stress (S) signal, and the attractive (A) signals. In low stress environment, N cells divide and emit signal A. In high stress environment, N cells stop emitting signal A. In presence of signal A, cells I climb up the concentration gradient of A. In absence of this signal, cells I move down the concentration gradient of signal S. I cell change their phenotype to N in low stress and high density environment, N cells change their phenotype to I in high stress and low density setting. This hypothesis is summarized in 

#figure(
  image("fig/CTMC.pdf"),
  caption:[Diagram of the proposed biological mechanism],
)

This model could explain the behavior observed empirically. Upon grafting of a tumor core composed of tightly packed N cells, S signal would rapidly build up and form a radial concentration gradient around the tumor. Cells at the periphery of the tumor, being less connected, would start expressing phenotype I and would start migrating radially from the tumor.
Removal of the tumorsphere would remove all significant potential source of S or A signal, which would result on random movement of remaining I cells.
Removal of a part of the tumorsphere would decrease the equilibrium concentration of S in the tumor. Remaining N cells would start emitting signal A, which will cause nearby I cells to switch direction, go back to the tumor and switch to phenotype N once arrived because of the low S signal and increased cell density.

Note that, if we assume that cells are not able to produce both a positive and
a negative gradient for the same chemical signal, we can consider that a single
chemical signal can carry 1 bit of information. Three different
configurations were shown to have different long range effect (1 : intact core;
2 : completely removed core; 3 : partially removed core), which corresponds to $log_2(3)$ bits of information,
which would require at least two chemical signals, under our previous assumption.
This shows that, while the proposed model is complex, it may satisfy the principle of parcimony.

= Model

The model takes place on a $W times H$ lattice, and we denote by $V$ the
Von-Neumann Neighborhood on this lattice, i.e. $(i,j) in V$ is $i$ and $j$
differ by 1 on only one of their coordinates.

== Continuous time Markov chain

The transition between the $I$ and $N$ phenotype is modelled by a continuous
time Markov chain, modified to have parameters varying with the environment.

$ PP(tau(c)_(t+ d t) = N | tau(c)_t = I ) = t_(I -> N) l(c) d t $

$ PP(tau(c)_(t+ d t) = I | tau(c)_t = N ) = t_(N -> I) (1-l(c)) d t $

where $l(c)$ is the proportion of the cell's perimeter that is contacting another cell, and
$t_(i -> j)$ is the base transition rate of cell type $i$ to cell type $j$.

== Hamiltonian

The Metropolis Hasting (MH) algorithm samples the system according to the Gibbs distribution
$ g prop e^(-H/T) $. For this reason, we need to define the Hamiltonian of the system $H$. It's the sum of different energetical terms :

$ H = H_"adh" + H_"area" + H_"perim" + H_"chemo" $

The adhesion term $H_"adh"$ models the adhesion of cells and prevent them from disolving intoeach others.

$ H_"adh" = sum_(|i,j|<1) J_(tau(i),tau(j)) bb(1)_(sigma(i) != sigma(j)) $

where $J_(i,j)$ is the potential associated with one unit of surface contact between a cell of type $i$ and a cell of type $j$.

We use an elastic potential $H_"area"$ to drive the cells toward a target area.

$ H_"area" = sum_(c in C) lambda_"area" (a(c) - a_t)^2 $

where $a(c)$ is the area of cell $c$ (i.e. number of lattice point taken by $c$), $a_t$ is a target area, and $lambda_"area"$ is the elastic modulus for the area.

We use an elastic potential $H_"perim"$ to drive the cells toward a target perimeter.

$ H_"perim" = sum_(c in C) lambda_"perim" (p(c) - p_t)^2 $

where $p(c)$ is the perimeter of cell $c$ (i.e. number of point of $c$ that have points not in $c$ in their 4-neighborhood), $p_t$ is a target perimeter, and $lambda_"perim"$ is the elastic modulus for the perimeter.

#v(20pt)

As in @szabo_cellular_2013, chemotaxis is incorporated using the method from
@savill_modelling_1997. During a MH step, the change in potential energy is
biased using the concentration $s$ of chemical signal. If the state of point $i$ is
being copied into point $j$, we define the biased energy change as

$ Delta H' = Delta H + xi_(tau(i)) (s(j) - s(i)) $

where $xi_(tau(i))$ is the strengh of the chemotactic motion on cell type $tau(i)$. $xi > 0$ corresponds to a movement down the concentration gradient, and $xi < 0$ up the concentration gradient.

Another way to incorporate chemotaxis that was experimented with is the following potential, which uses the everage concentration of the chemical signal inside the cell.

$ H_"chemo" = sum_(c in C) mu_(tau(c)) 1/a(c) sum_(i in c) s(i) $

== Diffusion equation

The concentration $s$ of chemical signal S is controlled by a diffusion partial differential equation.

$ partial_t s = D_s nabla s - d_s s + e_s bb(1)_N $

Where $D_s$ is the constant diffusion coefficient, $d_s$ is the rate of decay, and $e_s$ is the rate of emission. $bb(1)_N$ is a function equal to 1 inside $N$ cells and 0 everywhere else.

This equation was solved numerically on the same lattice as the Cellular-Potts model, using an explicit scheme.

== Simplified model

Because of time constraint, only a subset of this model could be implemented in time. The effect of signal $A$ was not considred in the implementation of the model.

= Implementation details

The model was implemented in C++, using the standard library, and raylib to
build a graphical user interface.

The core of the cellular potts model is the following MH algorithm.

#align(center,pad(rest:10pt,
pseudocode-list(booktabs: true,title:[A Metropolis Hasting step])[
  + $H$ is the value of the hamiltonian
  + $(x,y)$ is sampled uniformly on the lattice
  + $sigma = L(x,y)$ is the state of $(x,y)$ on the lattice $L$
  + $sigma'$ is sampled uniformly among the states of the neighbors of $(x,y)$
  + $L(x,y) <- sigma'$
  + $H'$ is the new value of the hamiltonian 
  + $r$ is sampled uniformly in $[0,1]$
  + if $r < exp(- (H' - H)/T)$
    + $H <- H'$
  + else
    + $L(x,y) <- sigma$
]))

The main challenge of the implementation is
to reduce the time required to compute the energy function at each iteration of
the MH algorithm. To do this, the attributes of each cells (area, perimeter,
etc.) are kept in memory. Then, during a transition of the system, only the
cells affected by the change are updated, keeping the time complexity of
lattice update constant.

The diffusion PDE is solved using the following explicit scheme :

$ s'_(x,y) = s_(x,y) + d t (D_s (s_(x-1,y) + s_(x+1,y) + s_(x,y-1) + s_(x,y+1) - 4 s_(x,y)) -d_s s_(x,y) + e_s bb(1)_(tau(x,y) = N)) $

At each time step of this scheme, the new value of $s$ is written in a back buffer, which is then switched with the front buffer, which allows a correct implementation of the numerical scheme with minimal data movement.

Updating the PDE is quite expensive compared to the rest of the implementation, which a complexity in $cal(O)(W times H)$. It is not needed to update the PDE after each transition of thesystem. Therefore the PDE is only updated with a specific frequency.

To initialize the tumor, a circular region of the lattice is set to be cells. Then, a basic k-mean algorithm is applied. This partitions the region according to a voronoi diagram in such a way that the sub-regions have similar size, each subregion is then initialized as a cell of type $N$.

= List of parameters

#figure(
table(columns : 3, align : left,
table.header([Parameter], [variable name], [role]),
$t_(I->N)$, [`t_ni`], [base rate of $I -> N$ transition],
$t_(N->I)$, [`t_in`], [base rate of $N -> I$ transition],
$a_t$,[`target_area`], [target area of the cells],
$lambda_"area"$,[`lambda_area`], [elastic modulus of the area force],
$p_t$,[`target_perim`],[target perimeter of the cells],
$lambda_"perim"$,[`lambda_perim`],[elastic modulus of the perimeter force],
$J_(I,I)$,[`J_ii`],[adhesion potential of one length unit between two $I$ cells],
$J_(N,N)$,[`J_nn`],[adhesion potential of one length unit between two $N$ cells],
$J_(I,N)$,[`J_in`],[adhesion potential of one length unit between a $I$ and a $N$ cell],
$D_s$,[`S_diffusion`],[Diffusion coefficient of signal $S$],
$d_s$,[`S_decay`],[decay rate of signal $S$],
$e_s$,[`S_emit`],[emission rate of signal $S$],
$xi_I$,[`xi`],[strengh of chemotactic motion of $I$ cells],
$mu_I$,[`mu_S_inv`],[strengh of chemotactic motion of $I$ cells]
),
caption : [List of model parameters])


#figure(table(columns : 3,align : left,
table.header([Parameter], [variable name], [role]),
$W$, [`lattice.width`], [width of the lattice],
$H$, [`lattice.height`], [height of the lattice],
$T$, [`T`], [temperature of the sampling procedure],
$d t$, [`S_dt`],[time discretization for the PDE solver],
[],[`S_step_frequ`],[frequency of update of the PDE],
[],[`CTMC_dt`],[time discretization of the continuous Markov chain]
),
caption : [List of numerical parameters])


#bibliography("refs.bib")
