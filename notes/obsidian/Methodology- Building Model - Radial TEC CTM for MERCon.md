---
Title: Methodology- Building Model - Radial TEC CTM for MERCon
date: 2026-04-13 20:59
type: publication
tags:
  - publication
  - MEMS
status:
---
We can copy paste existing one here. I see no need or reason to do modifications tbh.


\section{Formation of the Linear System}

To find the temperature values at each node in the thermal network, we can set up a system of linear equations based on the lumped resistances and heat sources. Each node will have an equation representing the balance of heat flow into and out of that node. Here we write the energy balance available 3 types of nodes: chip layer nodes, TEC layer nodes and central node.

\subsection{Chip Layer}

A general node in this layer has 4 terms: 2 lateral heat transfer, 1 vertical heat transfer to the TEC layer and the heat generation from the heat flux applied. Writing the energy balance equation,

\begin{align*}
  \underbrace{\frac{T_{\text{chip}, i-1} - T_{\text{chip}, i}}{R_{\text{chip}, i-1}}}_{\text{Lateral In}} & + \underbrace{\frac{T_{\text{chip}, i+1} - T_{\text{chip}, i}}{R_{\text{chip}, i}}}_{\text{Lateral Out}} \\ &- \underbrace{\frac{T_{\text{chip}, i} - T_{\text{TEC}, i}}{R_{\text{ve}, i}}}_{\text{Vertical Loss}} + Q_{\text{gen}, i} = 0
\end{align*}

Isolating the unknown temperature terms in the above equation gives,

\begin{equation}
  \small
  \begin{aligned}
     & \frac{1}{R_{\text{chip}, i-1}} T_{\text{chip}, i-1}-\left(
    \frac{1}{R_{\text{chip}, i-1}} + \frac{1}{R_{\text{chip}, i}} + \frac{1}{R_{\text{ve}, i}}
    \right) T_{\text{chip}, i}                                                                           \\
     & + \frac{1}{R_{\text{chip}, i}} T_{\text{chip}, i+1}+ \frac{1}{R_{\text{ve}, i}} T_{\text{TEC}, i}
    =- Q_{\text{gen}, i}.
  \end{aligned}
\end{equation}

This layer doesn't have any boundary conditions due to its closed nature. However, it should be noted that when writing the equation for the 1\textsuperscript{st} node, the term $R_{\text{chip},i-1} = R_{\text{chip},0}$ is equal to $R_{\text{chip}, 0\to 1}$ which is given by \cref{eq:R_chip_0_to_1}.

\subsection{TEC Layer}

The nodes in the TEC layer are governed by thermoelectric effect. The standard equation for the TEC module is as follows.

\begin{equation}
  Q_{c}=S I T_{c}-\textstyle{\frac{1}{2}}I^{2}R-K_t\left(T_{h}-T_{c}\right),
  \label{eq:TEC_standard}
\end{equation}

The terms represents heat pump through thermoelectric effects, joule heating due to electrical resistance, and back conduction respectively. Note that we have derived thermal resistance values instead of conductance values, so we need to use their reciprocals when substituting into the equation.


\subsubsection{Resistance Heat generation Assignment}

The 0.5 term in the joule heating term is to divide the heating from the TE legs to both cold side and hot side equally. However here we have interconnect and outerconnect that generates the heat. But their heat doesn't get added to both sides equally. So instead of a single lumped $R_{e,i}$, we have to  split the electrical resistance of Stage $i$ into three distinct physical components:

\begin{itemize}
  \item $R_{e,ic,i}$ (Inner/Cold Interconnect):The resistance of the copper trace at the cold junction (the node itself).
  \item $R_{e,TE-p,i}$ and $R_{e,TE-n,i}$ (TE Legs): The resistance of the P and N semiconductor pillars.
  \item $R_{e,oc,i}$ (Outer/Hot Interconnect): The resistance of the copper trace at the hot junction (the outer rim of this stage).
\end{itemize}

We assign Joule heating based on where the component is located relative to the node $T_{c,i}$.
\begin{itemize}
  \item At Node $T_{c,i}$: We have the Interconnect of Stage i.

        Contribution: 100\% of $I_i^2 R_{e,ic,i}$ adds heat to this node.

  \item  Between Node $T_{c,i}$ and Outer Boundary: We have the Legs of Stage i.

        Contribution: 50\% of $I_i^2 R_{e,TE-p,i}$ + 50\% of $I_i^2 R_{e,TE-n,i}$ flows back to this node (the cold side).

  \item  Coming from the Previous Stage $(i-1)$: The previous stage rejects heat into  current node. This rejected heat includes the generation from its own legs and its own outer interconnect.

        Contribution: 100\% of $I_{i-1}^2 R_{e,oc,i-1}$ (Outerconnect of prev stage) + 50\% of $I_{i-1}^2 R_{e,TE-p,i-1}$ + 50\% of $I_{i-1}^2 R_{e,TE-n,i-1}$.
\end{itemize}

Note that at the edge boundary, the heat generated from the outerconnect, and the half of p and n TE leg heating of the last stage will not be added to any node, as it is assumed to be dissipated to the boundary condition.

\subsubsection{Writing Energy Balance}

Writing the energy balance equation for a general TEC layer node,

\begin{align}
  \underbrace{\frac{T_{\text{chip},i} - T_{c,i}}{R_{\text{ve},i}}}_{\text{From Chip}}+\underbrace{Q_{h, i-1}}_{\text{From Stage } i-1} - \underbrace{Q_{c, i}}_{\text{Into Stage } i} = 0
  \label{eq:TEC_layer_energy_balance}
\end{align}

\paragraph{Heat Rejected from Stage $i-1$ - $Q_{h, i-1}$} This is the heat exiting the hot side of the inner ring. It includes the Peltier heat, the back conduction, and Specific Joule Terms:

\begin{equation}
  \small
  \begin{aligned}
    Q_{h,i-1} & = S_{i-1} I_{i-1} T_{\text{TEC},i-1} - K_{t,i-1}(T_{\text{TEC},i} - T_{\text{TEC},i-1}) \\ &+ \underbrace{\left[ \frac{1}{2}I_{i-1}^2 (R_{\text{TE-p},i-1}+R_{\text{TE-n},i-1})+  I_{i-1}^2 R_{e,\text{oc},i-1} \right]}_{\text{Joule Heat at Hot Side}}
  \end{aligned}
  \label{eq:Q_h_i-1}
\end{equation}

\paragraph{Heat Absorbed by stage $i$ - $Q_{c, i}$} This is the heat entering the cold side of the inner ring. It includes the Peltier heat, the back conduction, and Specific Joule Terms:

\begin{equation}
  \small
  \begin{aligned}
    Q_{c,i} & = S_{i} I_{i} T_{c,i} - K_{t,i}(T_{\text{TEC},i+1} - T_{\text{TEC},i}) \\ &- \underbrace{\left[ \frac{1}{2}I_{i}^2 (R_{\text{TE-p},i}+R_{\text{TE-n},i}) +  I_{i}^2 R_{e,\text{ic},i} \right]}_{\text{Joule Heat at Cold Side}}
    \label{eq:Q_c_i}
  \end{aligned}
\end{equation}

Note that we subtracted the Joule heating term because, in the context of the current node, the jouel heating term should be added to final equation. Since we are subtracting the $Q_c$ term, the sign of the jouel heating term must be negative to be added as a generation term in the final equation.
Substituting \cref{eq:Q_h_i-1} and \cref{eq:Q_c_i} into \cref{eq:TEC_layer_energy_balance} and isolating the unknown temperature terms gives:

\begin{equation}
  \small
  \begin{aligned}
     & \left(S_{i-1}I_{i-1}+K_{t,i-1}\right)T_{\text{TEC},i-1}-\left(\frac{1}{R_{\text{ve},i}}+K_{t,i-1}+S_{i}I_{i}+K_{t,i} \right)T_{\text{TEC},i} \\ &+K_{t,i}T_{\text{TEC},i+1}
    +\frac{1}{R_{\text{ve},i}}T_{\text{chip},i}                                                                                                     \\ &=-I_{i-1}^2 \left( \frac{(R_{\text{TE-p},i-1}+R_{\text{TE-n},i-1})}{2} + R_{\text{oc},i-1} \right) \\
     & \quad - I_{i}^2 \left( \frac{(R_{\text{TE-p},i}+R_{\text{TE-n},i})}{2} + R_{\text{ic},i} \right)
  \end{aligned}
\end{equation}

Similar to chip layer, when writing the equation for the 1\textsuperscript{st} node, the term $K_{t,i-1} = K_{t,0}$ is equal to $\frac{1}{R_{\text{TEC}, 0\to 1}}$ which is given by \cref{eq:R_TEC_0_to_1}. It should also be noted that the $K_{t}$ values are calculated using the lumped thermal resistance values derived in \cref{eq:TEC_element_Resistance}.

\subsubsection{Boundary Condition}

We assume a constant temperature boundary condition at the outer edge of the TEC layer to ease the modelling process. Since the micro channel characteristics aren't a part of the mathematical model, this boundary condition is justified. Therefore, for the last TEC layer node, the therm $T_{\text{TEC}ni+1}$ becomes a known value equal to the boundary temperature $T_{\text{w}}$. The energy balance equation for the last TEC layer node becomes:

\begin{equation*}
  \small
  \begin{aligned}
     & \left(S_{N-1}I_{N-1}+K_{t,N-1}\right)T_{\text{TEC},N-1}\\ &-\left(\frac{1}{R_{\text{ve},N}}+K_{t,N-1}+S_{N}I_{N}+K_{t,N} \right)T_{\text{TEC},N} +K_{t,N}T_{\text{w}}
    +\frac{1}{R_{\text{ve},N}}T_{\text{chip},N}                                                                                                     \\ &=-I_{N-1}^2 \left( \frac{(R_{\text{TE-p},N-1}+R_{\text{TE-p},N-1})}{2} + R_{\text{oc},N-1} \right) \\
     & \quad - I_{N}^2 \left( \frac{(R_{\text{TE-p},N}+R_{\text{TE-n},N})}{2} + R_{\text{ic},N} \right)
  \end{aligned}
\end{equation*}

Where $N$ is the number of stages in the radial TEC. Since $T_{\text{w}}$ is known, it can be carried out to the right hand side of the equation. Then the above equation becomes,

\begin{equation}
  \small
  \begin{aligned}
     & \left(S_{N-1}I_{N-1}+K_{t,N-1}\right)T_{\text{TEC},N-1}\\ &-\left(\frac{1}{R_{\text{ve},N}}+K_{t,N-1}+S_{N}I_{N}+K_{t,N} \right)T_{\text{TEC},N} +\frac{1}{R_{\text{ve},N}}T_{\text{chip},N}                                                                                                     \\ &=-I_{N-1}^2 \left( \frac{(R_{\text{TE-p},N-1}+R_{\text{TE-p},N-1})}{2} + R_{\text{oc},N-1} \right) \\
     & \quad - I_{N}^2 \left( \frac{(R_{\text{TE-p},N}+R_{\text{TE-n},N})}{2} + R_{\text{ic},N} \right) - K_{t,N} T_{\text{w}}
  \end{aligned}
\end{equation}

\subsection{Central Node}

At the exact center ($r=0$ to $r=r_{\text{cyl}}$), the "Chip Layer" and "TEC Layer" are mechanically and thermally merged by the heat spreader/via bundle. Therefore, they share the same temperature. Denoting this temperature as $T_{0}$, we can write the energy balance equation for this central node as:

\begin{equation*}
     \frac{T_{\text{chip},1} - T_{0}}{R_{\text{chip},0\to 1}} + \frac{T_{\text{TEC},1} - T_{0}}{R_{\text{TEC},0\to 1}} - Q_{\text{gen},0} = 0
\end{equation*}

Isolating the unknown temperature terms in the above equation gives,

\begin{equation}
  \begin{aligned}
     & -\left( \frac{1}{R_{\text{chip},0\to 1}} + \frac{1}{R_{\text{TEC},0\to 1}} \right) T_{0} + \frac{1}{R_{\text{chip},0\to 1}} T_{\text{chip},1}\\ &+ \frac{1}{R_{\text{TEC},0\to 1}} T_{\text{TEC},1} = - Q_{\text{gen},0}
  \end{aligned}
\end{equation}

The above equation will work as the first equation in the linear system.

\section{Matrix Formulation}

The system of linear equations derived from the energy balance equations can be represented in matrix form as:

\begin{equation}
  \mathbf{A} \cdot \mathbf{T} = \mathbf{B}
\end{equation}

Where $\mathbf{T}$ is the vector of unknown temperatures at each node, $\mathbf{A}$ is the coefficient matrix derived from the lumped resistances and thermoelectric parameters, and $\mathbf{B}$ is the vector representing the known terms including heat generation and boundary conditions.

\subsection{Block Notation}

\paragraph{The unknown vector} $\mathbf{T}$ can be organized in a block notation as follows:

\begin{equation}
  \small
  \mathbf{T} = \begin{bmatrix}
    T_{0}                   \\
    \mathbf{T_{\text{chip}}}\\
    \mathbf{T_{\text{TEC}}}
  \end{bmatrix}
\end{equation}

with,

\begin{equation}
  \mathbf{T_{\text{chip}}}=
  \begin{bmatrix}
    T_{\text{chip},1} \\
    T_{\text{chip},2} \\
    \vdots \\
    T_{\text{chip},N}
  \end{bmatrix}
  \quad \text{and} \quad
  \mathbf{T_{\text{TEC}}}=
  \begin{bmatrix}
    T_{\text{TEC},1} \\
    T_{\text{TEC},2} \\
    \vdots \\
    T_{\text{TEC},N}
  \end{bmatrix}
\end{equation}

\paragraph{The coefficient matrix} $\mathbf{A}$ can be organized in a block notation as follows:

\begin{equation}
   \mathbf{A} = 
   \begin{bmatrix} \mathbf{A_{00}} & \mathbf{Link_{0 \to Si}} & \mathbf{Link_{0 \to TEC}} \\ \mathbf{Link_{Si \to 0}} & \mathbf{A_{chip}} & \mathbf{A_{ve}} \\ \mathbf{Link_{TEC \to 0}} & \mathbf{A_{ve}} & \mathbf{A_{TEC}} 
  \end{bmatrix}
\end{equation}

Where: 
\begin{itemize}
\item $\mathbf{A_{00}}$ (Scalar): The sum of conductances leaving Node 0.
\item $\mathbf{Link_{0 \to Si}}$ : Coupling to the first chip ring.
\item $\mathbf{Link_{0 \to TEC}}$ : Coupling to the first TEC ring.
\item $\mathbf{A_{chip}}$ : a Tridiagonal matrix representing lateral conduction in the Chip. 
\item $\mathbf{A_{TEC}}$ : a Tridiagonal matrix representing the active TEC network. 
\item $\mathbf{A_{ve}}$ : a Diagonal matrix representing the vertical conduction between the Chip and TEC layers.
\end{itemize}

\paragraph{The known vector} $\mathbf{B}$ can be organized in a block notation as follows:

\begin{equation}
  \mathbf{B}= 
  \begin{bmatrix} 
    \mathbf{B_0} \\  \mathbf{B_{\text{chip}}} \\ \mathbf{B_{\text{TEC}}} 
  \end{bmatrix}
\end{equation}

Where :

\begin{itemize}
\item $\mathbf{B_0}$ (Scalar): The heat generation at Node 0
\item $\mathbf{B_{\text{chip}}}$ : A vector representing heat generation in the chip layer.
\item $\mathbf{B_{\text{TEC}}}$ : A vector representing the Joule heating in the TEC layer.
\end{itemize}

\subsection{Expanded Notation}

\subsubsection{$\mathbf{A_{00}}$ - Center Scalar}
$$ \mathbf{A_{00}} = \left[ -\left(\frac{1}{R_{\text{chip}, 0\to 1}} + \frac{1}{R_{\text{TEC}, 0\to 1}}\right) \right] $$
\textbf{Description:} This scalar represents the sum of all thermal conductances leaving the center node ($T_0$). It connects to the first Silicon ring and the first TEC node.

\subsubsection{$\mathbf{Link_{0 \to Si}}$ - Center to Silicon}
$$ \mathbf{Link_{0 \to Si}} = \begin{bmatrix} \frac{1}{R_{\text{chip}, 0\to 1}} & 0 & \dots & 0 \end{bmatrix}_{1 \times N} $$
\textbf{Description:} A row vector of size $1 \times N$. Only the first element is non-zero, representing the connection to $T_{\text{chip},1}$.

\subsubsection{$\mathbf{Link_{0 \to TEC}}$ - Center to TEC}
$$ \mathbf{Link_{0 \to TEC}} = \begin{bmatrix} \frac{1}{R_{\text{TEC}, 0\to 1}} & 0 & \dots & 0 \end{bmatrix}_{1 \times N} $$
\textbf{Description:} A row vector of size $1 \times N$. Only the first element is non-zero, representing the connection to $T_{\text{TEC},1}$.

\subsubsection{$\mathbf{Link_{Si \to 0}}$ - Silicon to Center}
$$ \mathbf{Link_{Si \to 0}} = \begin{bmatrix} \frac{1}{R_{\text{chip}, 0\to 1}} \\ 0 \\ \vdots \\ 0 \end{bmatrix}_{N \times 1} $$
\textbf{Description:} A column vector of size $N \times 1$. This is the transpose of Region 2.

\subsubsection{$\mathbf{A_{chip}}$ - Lateral Conduction Matrix}
\textbf{Description:} A symmetric \textbf{Tridiagonal} matrix representing passive heat spreading.
\begin{itemize}
    \item \textbf{Diagonal $(i,i)$:} The negative sum of all conductances leaving the node.
    \item \textit{Boundary:} Note that at $i=1$, the term connects to $T_0$. At $i=N$, the term assumes an adiabatic edge (no $R_{\text{chip},N}$).
    \item \textbf{Off-Diagonals:} Represents the lateral resistance $R_{\text{chip}}$ between rings.
\end{itemize}

\subsubsection{$\mathbf{A_{ve}}$ - Vertical Coupling}
$$ \mathbf{A_{ve}} = \begin{bmatrix} \frac{1}{R_{\text{ve},1}} & 0 & \dots & 0 \\ 0 & \frac{1}{R_{\text{ve},2}} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & \frac{1}{R_{\text{ve},N}} \end{bmatrix}_{N \times N} $$
\textbf{Description:} A \textbf{Diagonal} matrix.
\begin{itemize}
    \item This matrix appears twice in the global system: once connecting Silicon to TEC (Region 6), and once connecting TEC to Silicon (Region 8).
\end{itemize}

\subsubsection{$\mathbf{Link_{TEC \to 0}}$ - TEC to Center}
$$ \mathbf{Link_{TEC \to 0}} = \begin{bmatrix} \frac{1}{R_{\text{TEC}, 0\to 1}} \\ 0 \\ \vdots \\ 0 \end{bmatrix}_{N \times 1} $$
\textbf{Description:} A column vector of size $N \times 1$. This is the transpose of Region 3.

\subsubsection{$\mathbf{A_{TEC}}$ - Active Pumping Matrix}
\textbf{Description:} An \textbf{Asymmetric Tridiagonal} matrix.
\begin{itemize}
    \item \textbf{Diagonal $(i,i)$:} Negative sum of Vertical Loss + Back conduction IN + Peltier OUT + Back conduction OUT.
    \item \textbf{Lower Off-Diagonal $(i, i-1)$:} Heat arriving from the previous stage. This contains the Peltier term $S_{i-1}I_{i-1}$ and the back conduction $K_{t,i-1}$.
    \item \textbf{Upper Off-Diagonal $(i, i+1)$:} Only back conduction $K_{t,i}$ comes from the hotter stage. The Peltier term does not flow "backwards".
\end{itemize}

\subsubsection{Vector $\mathbf{B}$ - Source Vector}

\textbf{Description:}
\begin{itemize}
    \item $\mathbf{B_{0}}$: Heat generation at the center.
    \item $\mathbf{B_{\text{chip}}}$: Heat generation in the Silicon rings (inverted sign).
    \item $\mathbf{B_{\text{TEC}}}$: Contains the Joule heating terms.
    \item Current Node heating: $I_i^2 R_{\text{ic},i}$ (Interconnect) + $0.5 \cdot I_i^2 (R_{\text{TE-p},i}+R_{\text{TE-n},i})$ (Half leg).
    \item Previous Stage heating: $I_{i-1}^2 R_{\text{oc},i-1}$ (Outerconnect) + $0.5 \cdot I_{i-1}^2 (R_{\text{TE-p},i-1}+R_{\text{TE-n},i-1})$ (Half leg).
    \item \textbf{Last Row:} Includes the boundary condition term $-K_{t,N} T_{\text{w}}$.
\end{itemize}



---
# 1. Exploration


# 2. References
