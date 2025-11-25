
To solve this numerically, we treat this as a **Finite Difference** problem (which is mathematically equivalent to a resistor network with current sources). We must solve for the temperatures simultaneously because the Peltier effect depends on temperature ($\alpha I T$), making the heat flow dependent on the variable we are solving for. 

--- 
# 1 Definitions & Discretization 

We divide the system radially into $N$ concentric rings (Nodes $1$ to $N$). 

**Variables:** 
* $\mathbf{T_{Si}}$: Vector of Chip Temperatures ($N \times 1$) 
* $\mathbf{T_{c}}$: Vector of TEC Node Temperatures ($N \times 1$). 
	* *Note:* In this radial model, Node $T_{c,i}$ represents the junction between Stage $i-1$ (Hot) and Stage $i$ (Cold). 
	* Node $1$ is the innermost cold side (Center). Node $N+1$ is the outermost hot side (Water edge). To keep matrices square for optimization, we usually treat the boundary $N+1$ as a known boundary condition or include it as the final node. Let's assume **$N$ active nodes**, with boundary conditions handling the edges. 

**Constants for Stage $i$:** 
* $R_{lat}$: Lateral Thermal Resistance of Silicon ring. 
* $R_{v,i}$: Vertical Resistance (The "Switch" for TSVs). 
	* Zone 1: $R_{v,i} \approx 0$ (low) 
	* Zone 2: $R_{v,i} \approx \infty$ (high) 
* $K_i$: Thermal conductance of TEC stage $i$ ($K = \frac{k A}{L}$).
* $S_i$: Seebeck coefficient term ($\alpha \times N_{couples}$). 
* $R_{e,i}$: Electrical resistance of stage $i$. 
* $I_i$: Current supplied to stage $i$. 

--- 
# 2 General Nodal Equations 

We write the energy balance ($\sum Q_{in} = 0$) for a generic node $i$. 
## 2.1 Layer 1: The Chip Node ($T_{Si,i}$) 

Heat enters from the inner silicon ring, leaves to the outer silicon ring, enters from generation, and leaves vertically to the TEC. 
$$\underbrace{\frac{T_{s, i-1} - T_{s, i}}{R_{lat, i-1}}}_{\text{Lateral In}} + \underbrace{\frac{T_{s, i+1} - T_{s, i}}{R_{lat, i}}}_{\text{Lateral Out}} - \underbrace{\frac{T_{s, i} - T_{c, i}}{R_{v, i}}}_{\text{Vertical Loss}} + Q_{gen, i} = 0$$
$$ \frac{T_{Si,i-1} - T_{Si,i}}{R_{lat,i-1}} + \frac{T_{Si,i+1} - T_{Si,i}}{R_{lat,i}}  - \frac{T_{Si,i} - T_{c,i}}{R_{v,i}} + Q_{gen,i}= 0 $$
**Isolating the Unknowns ($T$):** 
$$ \left( \frac{1}{R_{lat,i-1}} \right)T_{Si,i-1} - \left( \frac{1}{R_{lat,i-1}} + \frac{1}{R_{lat,i}} + \frac{1}{R_{v,i}} \right)T_{Si,i} + \left( \frac{1}{R_{lat,i}} \right)T_{Si,i+1} + \left( \frac{1}{R_{v,i}} \right)T_{c,i} = -Q_{gen,i} $$
## 2.2 Layer 2: The TEC Node ($T_{c,i}$) 

This node represents the thermal mass of Stage $i$. It receives heat from the previous stage ($i-1$), sends heat to the next stage pumping ($i$), and receives heat vertically from the Chip. We consider a node as the interconnect copper junction. 

The standard equation for the TEC module is as follows.
$$Q_{c}=\alpha I T_{c}-\textstyle{\frac{1}{2}}I^{2}R-K\left(T_{h}-T_{c}\right),$$
The terms represents heat pump, joule heating, and back conduction respectively. the 0.5 term in the joule heating term is to divide the heating from the TE legs to both cold side and hot side equally. However here we have interconnect and outerconnect that generates the heat. But their heat doesn't get added to both sides equally. Therefore we cannot add a $R_{total}$ term and call it a day. Instead of a single lumped $R_{e,i}$, we have to  split the electrical resistance of **Stage i** into three distinct physical components:

1. **$R_{ic,i}$ (Inner/Cold Interconnect):** The resistance of the copper trace at the cold junction (the node itself).
2. **$R_{leg,i}$ (TE Legs):** The resistance of the P and N semiconductor pillars.
3. **$R_{oc,i}$ (Outer/Hot Interconnect):** The resistance of the copper trace at the hot junction (the outer rim of this stage).

We assign Joule heating based on where the component is located relative to the node $T_{c,i}$.
- **At Node $T_{c,i}$:** You have the **Interconnect of Stage i**.
    - Contribution: **100\%** of $I_i^2 R_{ic,i}$ adds heat to this node.
        
- **Between Node $T_{c,i}$ and Outer Boundary:** You have the **Legs of Stage i**.
    - Contribution: **50\%** of $I_i^2 R_{leg,i}$ flows back to this node (the cold side).
        
- **Coming from the Previous Stage $(i-1)$:** The previous stage rejects heat into your current node. This rejected heat includes the generation from its own legs and its own outer interconnect.
    - Contribution: **100\%** of $I_{i-1}^2 R_{oc,i-1}$ (Outerconnect of prev stage) + **50\%** of $I_{i-1}^2 R_{leg,i-1}$.


**Balance Equation ($Q_{vert} + Q_{h,i-1} - Q_{c,i} = 0$):** 
$$Q_{in} + Q_{generated} - Q_{out} = 0$$
$$\underbrace{\frac{T_{s,i} - T_{c,i}}{R_{v,i}}}_{\text{From Chip}}+\underbrace{Q_{h, i-1}}_{\text{From Stage } i-1} - \underbrace{Q_{c, i}}_{\text{Into Stage } i} = 0$$

**Heat Rejected from Stage $i-1$ ($Q_{h, i-1}$)**
This is the heat exiting the hot side of the inner ring. It includes the Peltier heat, the back conduction, and Specific Joule Terms:
$$Q_{h,i-1} = S_{i-1} I_{i-1} T_{c,i-1} - K_{i-1}(T_{c,i} - T_{c,i-1}) + \underbrace{\left[ \frac{1}{2}I_{i-1}^2 R_{leg,i-1} + 1.0 \cdot I_{i-1}^2 R_{oc,i-1} \right]}_{\text{Joule Heat at Hot Side}}$$
**Heat Absorbed by Stage $i$ ($Q_{c, i}$)**
This is the "Net Cooling Capacity". It is the Peltier cooling minus the heat generated locally and minus back conduction:
$$Q_{c,i} = S_{i} I_{i} T_{c,i} + K_{i}(T_{c,i+1} - T_{c,i}) - \underbrace{\left[ \frac{1}{2}I_{i}^2 R_{leg,i} + 1.0 \cdot I_{i}^2 R_{ic,i} \right]}_{\text{Joule Heat at Cold Side}}$$
The Combined Equation Becomes,

$$\begin{align} \frac{T_{Si,i} - T_{c,i}}{R_{v,i}} + \left[S_{i-1} I_{i-1} T_{c,i-1} - K_{i-1}(T_{c,i} - T_{c,i-1}) + \left( \frac{1}{2}I_{i-1}^2 R_{leg,i-1} + I_{i-1}^2 R_{oc,i-1} \right)\right] &\\- \left[ S_{i} I_{i} T_{c,i} + K_{i}(T_{c,i+1} - T_{c,i}) - \left( \frac{1}{2}I_{i}^2 R_{leg,i} + I_{i}^2 R_{ic,i} \right) \right] = 0 \end{align}$$
When terms are isolated,
$$\begin{align}\left(S_{i-1}I_{i-1}+K_{i-1}\right)T_{c,i-1}-\left(\frac{1}{R_{v,i}}+K_{i-1}+S_{i}I_{i}-K_i \right)T_{c,i}+K_iT_{c,i+1}&\\+\frac{1}{R_{v,i}}T_{Si,i}=-I_{i-1}^2 \left( \frac{R_{leg,i-1}}{2} + R_{oc,i-1} \right) - I_{i}^2 \left( \frac{R_{leg,i}}{2} + R_{ic,i} \right)\end{align}$$

---
# 3 Boundary Condition Treatment

## 3.1 Node 0

As discussed in [[TEC Thermal Network#3 Node 0 treatment|this section]] , the $T_0$ for both layers are the same., and therefore a separate equation is required to add that coupling and the associated heat generation term. Therefore, following equation will become the entries for the first row,
$$\left(\frac{1}{R_{Si,0\to 1}}+\frac{1}{R_{TEC,0\to 1}}\right)T_0-\frac{1}{R_{Si,0\to 1}}T_{Si,1}-\frac{1}{R_{c,0\to 1}}T_{c,1}=Q_{gen,0}$$
## 3.2 Last Node - Boundary Condition

A detailed description on the reasoning behind these equations are shown in [[Basic Radial TEC Design#4 Boundary Condition|this]] section. Here the equations get separated into the standard linear form.
### 3.2.1 Silicon / Chip Layer
Only the $T_{Si,i+1}$ term vanishes from this layer equation and the resulting equation looks as follows.
$$\frac{T_{Si,N-1} - T_{Si,N}}{R_{lat, N-1}} - \frac{T_{Si,N} - T_{c,N}}{R_{v,N}} + Q_{gen,N} = 0$$
When we separate into temperature terms,
$$ \left( \frac{1}{R_{lat,N-1}} \right)T_{Si,N-1} - \left( \frac{1}{R_{lat,N-1}}  + \frac{1}{R_{v,N}} \right)T_{Si,N} + \left( \frac{1}{R_{v,N}} \right)T_{c,N} = -Q_{gen,i} $$
### 3.2.2 TEC Layer

#### 3.2.2.1 Constant Temperature BC
For this condition $T_c,N+1\to T_w$, and the equation becomes,
$$\begin{align}\left(S_{N-1}I_{N-1}+K_{N-1}\right)T_{c,N-1}-\left(\frac{1}{R_{v,N}}+K_{N-1}+S_{N}I_{N}-K_N \right)T_{c,N}+K_NT_{w}&\\+\frac{1}{R_{v,N}}T_{Si,N}=-I_{N-1}^2 \left( \frac{R_{leg,N-1}}{2} + R_{oc,N-1} \right) - I_{i}^2 \left( \frac{R_{leg,N}}{2} + R_{ic,N} \right)\end{align}$$
The known term can be carried out to the RHS.
$$\begin{align}\left(S_{N-1}I_{N-1}+K_{N-1}\right)T_{c,N-1}-\left(\frac{1}{R_{v,N}}+K_{N-1}+S_{N}I_{N}-K_N \right)T_{c,N}+\frac{1}{R_{v,N}}T_{Si,N}&\\=-K_NT_{w}-I_{N-1}^2 \left( \frac{R_{leg,N-1}}{2} + R_{oc,N-1} \right) - I_{N}^2 \left( \frac{R_{leg,N}}{2} + R_{ic,N} \right)\end{align}$$


--- 
# 4 The Matrix Structure 

We form the linear system $\mathbf{M} \mathbf{x} = \mathbf{B}$. Since we have two layers, we use a **Block Matrix** approach. 
$$ \mathbf{x} = \begin{bmatrix} \mathbf {T_0} \\ \mathbf{T_{Si}} \\ \mathbf{T_{c}} \end{bmatrix} $$
## 4.1 Block Notation 
$$ \begin{bmatrix} \mathbf{M_{00}} & \mathbf{Link_{0 \to Si}} & \mathbf{Link_{0 \to TEC}} \\ \mathbf{Link_{Si \to 0}} & \mathbf{A_{Silicon}} & \mathbf{A_{Vertical}} \\ \mathbf{Link_{TEC \to 0}} & \mathbf{A_{Vertical}} & \mathbf{A_{TEC}} \end{bmatrix} \begin{bmatrix} \mathbf{T_0} \\ \mathbf{T_{Si}} \\ \mathbf{T_{c}} \end{bmatrix} = \begin{bmatrix} \mathbf{B_0} \\ \mathbf{B_{Si}} \\ \mathbf{B_{TEC}} \end{bmatrix} $$
Where: 
* $\mathbf{M_{00}}$ (Scalar): The sum of conductances leaving Node 0.
* $\mathbf{Link_{0 \to Si}}$ : Coupling to the first Silicon ring.
* $\mathbf{Link_{0 \to TEC}}$ : TCoupling to the first TEC ring.
* $\mathbf{A_{Silicon}}$ : a **Tridiagonal** matrix representing lateral conduction in the Chip. 
* $\mathbf{A_{TEC}}$ : a **Tridiagonal** matrix representing the active TEC network. 
* $\mathbf{A_{Vertical}}$ : a **Diagonal** matrix representing the vertical TSV connections. 


---
## 4.2 Detailed Coefficients 

### 4.2.1 Region 1: $\mathbf{M_{00}}$ (Center Scalar) 
$$ \mathbf{M_{00}} = \left[ \frac{1}{R_{Si, 0\to 1}} + \frac{1}{R_{c, 0\to 1}} \right] $$ **Description:** This scalar represents the sum of all thermal conductances leaving the center node ($T_0$). It connects to the first Silicon ring and the first TEC node. 
### 4.2.2 Region 2: $\mathbf{Link_{0 \to Si}}$ (Center to Silicon) 
$$ \mathbf{Link_{0 \to Si}} = \begin{bmatrix} -\frac{1}{R_{Si, 0\to 1}} & 0 & \dots & 0 \end{bmatrix}_{1 \times N} $$
**Description:** A row vector of size $1 \times N$. Only the first element is non-zero, representing the connection to $T_{Si,1}$. * **Loop Logic:** Set index `[0, 0]` to $-\frac{1}{R_{Si, 0\to 1}}$. 
### 4.2.3 Region 3: $\mathbf{Link_{0 \to TEC}}$ (Center to TEC) 
$$ \mathbf{Link_{0 \to TEC}} = \begin{bmatrix} -\frac{1}{R_{c, 0\to 1}} & 0 & \dots & 0 \end{bmatrix}_{1 \times N} $$**Description:** A row vector of size $1 \times N$. Only the first element is non-zero, representing the connection to $T_{c,1}$. * **Loop Logic:** Set index `[0, 0]` to $-\frac{1}{R_{c, 0\to 1}}$.

### 4.2.4 Region 4: $\mathbf{Link_{Si \to 0}}$ (Silicon to Center) 
$$ \mathbf{Link_{Si \to 0}} = \begin{bmatrix} -\frac{1}{R_{Si, 0\to 1}} \\ 0 \\ \vdots \\ 0 \end{bmatrix}_{N \times 1} $$ **Description:** A column vector of size $N \times 1$. This is the transpose of Region 2. * **Loop Logic:** Set index `[0, 0]` to $-\frac{1}{R_{Si, 0\to 1}}$. 
### 4.2.5 Region 5: $\mathbf{A_{Silicon}}$ (Lateral Conduction Matrix) 
$$ \mathbf{A_{Si}} = \begin{bmatrix} -\left( \frac{1}{R_{Si,0\to 1}} + \frac{1}{R_{lat,1}} + \frac{1}{R_{v,1}} \right) & \frac{1}{R_{lat,1}} & 0 & \dots \\ \frac{1}{R_{lat,1}} & -\left( \frac{1}{R_{lat,1}} + \frac{1}{R_{lat,2}} + \frac{1}{R_{v,2}} \right) & \frac{1}{R_{lat,2}} & \dots \\ \vdots & \ddots & \ddots & \vdots \\ 0 & \dots & \frac{1}{R_{lat,N-1}} & -\left( \frac{1}{R_{lat,N-1}} + \frac{1}{R_{v,N}} \right) \end{bmatrix}_{N \times N} $$ **Description:** A symmetric **Tridiagonal** matrix representing passive heat spreading. 
* **Diagonal $(i,i)$:** The negative sum of all conductances leaving the node. 
* *Boundary:* Note that at $i=1$, the term connects to $T_0$. At $i=N$, the term assumes an adiabatic edge (no $R_{lat,N}$). 
* **Off-Diagonals:** Represents the lateral resistance $R_{lat}$ between rings. 
### 4.2.6 Region 6 & 8: $\mathbf{A_{Vertical}}$ (Vertical Coupling) 
$$ \mathbf{A_{Vert}} = \begin{bmatrix} \frac{1}{R_{v,1}} & 0 & \dots & 0 \\ 0 & \frac{1}{R_{v,2}} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & \frac{1}{R_{v,N}} \end{bmatrix}_{N \times N} $$ **Description:** A **Diagonal** matrix. 
* This matrix appears twice in the global system: once connecting Silicon to TEC (Region 6), and once connecting TEC to Silicon (Region 8). 
* **Loop Logic:** For each $i$, set `[i, i]` to $\frac{1}{R_{v,i}}$. This is the "Switch" variable—if a zone has no TSVs, $R_{v,i}$ is set to a very large number (approx $\infty$), making this term effectively 0. 
### 4.2.7 Region 7: $\mathbf{Link_{TEC \to 0}}$ (TEC to Center) 
$$ \mathbf{Link_{TEC \to 0}} = \begin{bmatrix} -\frac{1}{R_{c, 0\to 1}} \\ 0 \\ \vdots \\ 0 \end{bmatrix}_{N \times 1} $$ **Description:** A column vector of size $N \times 1$. This is the transpose of Region 3. * **Loop Logic:** Set index `[0, 0]` to $-\frac{1}{R_{c, 0\to 1}}$. 
### 4.2.8 Region 9: $\mathbf{A_{TEC}}$ (Active Pumping Matrix) 
$$ \tiny \mathbf{A_{TEC}} = \begin{bmatrix} -\left(\frac{1}{R_{v,1}} + \frac{1}{R_{c,0\to 1}} + S_1 I_1 - K_1 \right) & K_1 & 0 & \dots \\ (S_1 I_1 + K_1) & -\left(\frac{1}{R_{v,2}} + K_1 + S_2 I_2 - K_2 \right) & K_2 & \dots \\ \vdots & \ddots & \ddots & \vdots \\ 0 & \dots & (S_{N-1} I_{N-1} + K_{N-1}) & -\left(\frac{1}{R_{v,N}} + K_{N-1} + S_N I_N - K_N \right) \end{bmatrix}_{N \times N} $$**Description:** An **Asymmetric Tridiagonal** matrix. 
* **Diagonal $(i,i)$:** Negative sum of Vertical Loss + Back conduction IN + Peltier OUT + Back conduction OUT. 
* **Lower Off-Diagonal $(i, i-1)$:** Heat arriving from the previous stage. This contains the Peltier term $S_{i-1}I_{i-1}$ and the back conduction $K_{i-1}$. 
* **Upper Off-Diagonal $(i, i+1)$:** Only back conduction $K_i$ comes from the hotter stage. The Peltier term does not flow "backwards". 
### 4.2.9 Vector $\mathbf{B}$ (Source Vector) 
$$ \mathbf{B} = \begin{bmatrix} Q_{gen,0} \\ \hline -Q_{gen,1} \\ \vdots \\ -Q_{gen,N} \\ \hline -I_{1}^2 \left( \frac{R_{leg,1}}{2} + R_{ic,1} \right) \\ \vdots \\ -K_N T_{w} - I_{N-1}^2 \left( \frac{R_{leg,N-1}}{2} + R_{oc,N-1} \right) - I_{N}^2 \left( \frac{R_{leg,N}}{2} + R_{ic,N} \right) \end{bmatrix}_{(2N+1) \times 1} $$**Description:** 
* **$\mathbf{B_{0}}$:** Heat generation at the center. 
* **$\mathbf{B_{Si}}$:** Heat generation in the Silicon rings (inverted sign). 
* **$\mathbf{B_{TEC}}$:** Contains the Joule heating terms. 
* Current Node heating: $I_i^2 R_{ic,i}$ (Interconnect) + $0.5 \cdot I_i^2 R_{leg,i}$ (Half leg). 
* Previous Stage heating: $I_{i-1}^2 R_{oc,i-1}$ (Outerconnect) + $0.5 \cdot I_{i-1}^2 R_{leg,i-1}$ (Half leg). 
* **Last Row:** Includes the boundary condition term $-K_N T_w$.

--- 
## 4.3 Expanded Matrix 

The LHS Matrix ($\mathbf{M}$) Expanded Form:
$$\tiny \mathbf{M} = \left[ \begin{array}{c|ccccc|ccccc} 
% REGION 1: Node 0 Self 
\left(\frac{1}{R_{Si,0\to 1}}+\frac{1}{R_{TEC,0\to 1}}\right) & 
% REGION 2: Node 0 -> Si 
-\frac{1}{R_{Si,0\to 1}} & \dots & 0 & \dots & 0 & 
% REGION 3: Node 0 -> TEC 
-\frac{1}{R_{c,0\to 1}} & \dots & 0 & \dots & 0 \\ \hline 
% % REGION 4: Si -> Node 0 % Si Row 1 
-\frac{1}{R_{Si,0\to 1}} & 
% REGION 5: Si Internal (Tri-diagonal) 
-\left( \frac{1}{R_{Si,0\to 1}} + \frac{1}{R_{lat,1}} + \frac{1}{R_{v,1}} \right) & \frac{1}{R_{lat,1}} & 0 & \dots & 0 & 
% REGION 6: Si -> TEC (Diagonal) 
\frac{1}{R_{v,1}} & 0 & 0 & \dots & 0 \\ 
% Si Row General 
\vdots & \frac{1}{R_{lat,i-1}} & -\left( \frac{1}{R_{lat,i-1}} + \frac{1}{R_{lat,i}} + \frac{1}{R_{v,i}} \right) & \frac{1}{R_{lat,i}} & \dots & 0 & 0 & \frac{1}{R_{v,i}} & 0 & \dots & 0 \\ 
% Si Row N 
0 & 0 & \dots & \frac{1}{R_{lat,N-1}} & -\left( \frac{1}{R_{lat,N-1}} + \frac{1}{R_{v,N}} \right) & 0 & 0 & 0 & \dots & \frac{1}{R_{v,N}} \\ \hline 
% % REGION 7: TEC -> Node 0 % TEC Row 1
-\frac{1}{R_{c,0\to 1}} & 
% REGION 8: TEC -> Si (Diagonal) 
\frac{1}{R_{v,1}} & 0 & 0 & \dots & 0 & 
% REGION 9: TEC Internal (Tri-diagonal)
-\left(\frac{1}{R_{v,1}} + \frac{1}{R_{c,0\to 1}} + S_1 I_1 - K_1 \right) & K_1 & 0 & \dots & 0 \\ 
% TEC Row General 
\vdots & 0 & \frac{1}{R_{v,i}} & 0 & \dots & 0 & (S_{i-1}I_{i-1}+K_{i-1}) & -\left(\frac{1}{R_{v,i}}+K_{i-1}+S_{i}I_{i}-K_i \right) & K_i & \dots & 0 \\ 
% TEC RowN 
0 & 0 & 0 & \dots & 0 & \frac{1}{R_{v,N}} & 0 & 0 & \dots & (S_{N-1}I_{N-1}+K_{N-1}) & -\left(\frac{1}{R_{v,N}}+K_{N-1}+S_{N}I_{N}-K_N \right) \end{array} \right]
$$

Variable Vector:
$$ \mathbf{T} = \begin{bmatrix} \mathbf{T_0} \\ \hline \mathbf{T_{Si}} \\ \hline \mathbf{T_{TEC}} \end{bmatrix} =   \left[ \begin{array}{c} T_0 \\ \hline T_{Si,1} \\ \vdots \\ T_{Si,i} \\ \vdots \\ T_{Si,N} \\ \hline T_{c,1} \\ \vdots \\ T_{c,i} \\ \vdots \\ T_{c,N} \end{array} \right] $$
The RHS Vector ($\mathbf{B}$) 
$$\mathbf{B}= \begin{bmatrix} \mathbf{B_0} \\ \hline \mathbf{B_{Si}} \\ \hline \mathbf{B_{TEC}} \end{bmatrix} =  =  \left[ \begin{array}{c} Q_{gen,0} \\ \hline -Q_{gen,1} \\ \vdots \\ -Q_{gen,i} \\ \vdots \\ -Q_{gen,N} \\ \hline -I_{0}^2 \left( \frac{R_{leg,0}}{2} + R_{oc,0} \right) - I_{1}^2 \left( \frac{R_{leg,1}}{2} + R_{ic,1} \right) \\ \vdots \\ -I_{i-1}^2 \left( \frac{R_{leg,i-1}}{2} + R_{oc,i-1} \right) - I_{i}^2 \left( \frac{R_{leg,i}}{2} + R_{ic,i} \right) \\ \vdots \\ -K_NT_{w}-I_{N-1}^2 \left( \frac{R_{leg,N-1}}{2} + R_{oc,N-1} \right) - I_{N}^2 \left( \frac{R_{leg,N}}{2} + R_{ic,N} \right) \end{array} \right] $$

--- 

# 5 Example 3 Stage System

Here is the **full expanded matrix formulation** for a 3-Stage, 3-Node system ($N=3$). This represents the specific case where $i=1, 2, 3$, incorporating the Center Node ($T_0$) and the Water Boundary Condition ($T_w$) explicitly. The system size is **7x7**: 
* 1 Row for Center Node ($T_0$) 
* 3 Rows for Silicon ($T_{Si,1}, T_{Si,2}, T_{Si,3}$) 
* 3 Rows for TEC ($T_{c,1}, T_{c,2}, T_{c,3}$) 
### 5.1.1 Variable & Source Vectors 
$$ \mathbf{x} = \begin{bmatrix} T_0 \\ T_{Si,1} \\ T_{Si,2} \\ T_{Si,3} \\ T_{c,1} \\ T_{c,2} \\ T_{c,3} \end{bmatrix} \quad \quad \mathbf{B} = \begin{bmatrix} Q_{gen,0} \\ -Q_{gen,1} \\ -Q_{gen,2} \\ -Q_{gen,3} \\ - I_{1}^2 \left( \frac{R_{leg,1}}{2} + R_{ic,1} \right) \\ -I_{1}^2 \left( \frac{R_{leg,1}}{2} + R_{oc,1} \right) - I_{2}^2 \left( \frac{R_{leg,2}}{2} + R_{ic,2} \right) \\ -K_3 T_w -I_{2}^2 \left( \frac{R_{leg,2}}{2} + R_{oc,2} \right) - I_{3}^2 \left( \frac{R_{leg,3}}{2} + R_{ic,3} \right) \end{bmatrix} $$ --- 
### 5.1.2 The Global Conductance Matrix (7x7) 

This matrix is divided into the **9 Regions**. Note that zeros are omitted in empty spots for clarity, but the structure preserves the tridiagonal/diagonal nature. 

$$ \tiny \mathbf{M} = \left[ \begin{array}{c|ccc|ccc} 
% ===================== ROW 1: NODE 0 ===================== 
% REGION 1 (0->0) 
\left(\frac{1}{R_{Si,0\to 1}}+\frac{1}{R_{c,0\to 1}}\right) & 
% REGION 2 (0->Si) 
-\frac{1}{R_{Si,0\to 1}} & 0 & 0 & 
% REGION 3 (0->TEC) 
-\frac{1}{R_{c,0\to 1}} & 0 & 0 \\ \hline 
% % ===================== ROWS 2-4: SILICON LAYER ===================== 
% REGION 4 (Si->0) 
-\frac{1}{R_{Si,0\to 1}} & 
% REGION 5 (Si Internal - Tridiagonal) 
-\left( \frac{1}{R_{Si,0\to 1}} + \frac{1}{R_{lat,1}} + \frac{1}{R_{v,1}} \right) & \frac{1}{R_{lat,1}} & 0 & 
% REGION 6 (Si->TEC - Diagonal) 
\frac{1}{R_{v,1}} & 0 & 0 \\ 
% 
0 & \frac{1}{R_{lat,1}} & -\left( \frac{1}{R_{lat,1}} + \frac{1}{R_{lat,2}} + \frac{1}{R_{v,2}} \right) & \frac{1}{R_{lat,2}} & 0 & \frac{1}{R_{v,2}} & 0 \\ 
%
0 & 0 & \frac{1}{R_{lat,2}} & -\left( \frac{1}{R_{lat,2}} + \frac{1}{R_{v,3}} \right) & 0 & 0 & \frac{1}{R_{v,3}} \\ \hline 
% % ===================== ROWS 5-7: TEC LAYER ===================== 
% REGION 7 (TEC->0) 
-\frac{1}{R_{c,0\to 1}} & 
% REGION 8 (TEC->Si - Diagonal) 
\frac{1}{R_{v,1}} & 0 & 0 & 
% REGION 9 (TEC Internal - Tridiagonal) 
-\left(\frac{1}{R_{v,1}} + \frac{1}{R_{c,0\to 1}} + S_1 I_1 - K_1 \right) & K_1 & 0 \\ 
% 
0 & 0 & \frac{1}{R_{v,2}} & 0 & (S_1 I_1 + K_1) & -\left(\frac{1}{R_{v,2}}+K_1+S_2 I_2 - K_2 \right) & K_2 \\ 
% 
0 & 0 & 0 & \frac{1}{R_{v,3}} & 0 & (S_2 I_2 + K_2) & -\left(\frac{1}{R_{v,3}}+K_2+S_3 I_3 - K_3 \right) \end{array} \right] 
$$ 

### 5.1.3 Region Description for the 3-Node System

1. **$\mathbf{M_{00}}$ (Row 1, Col 1):** 
	Sum of conductances leaving the center node to both Silicon ring 1 and TEC node 1. 
2. **$\mathbf{Link_{0 \to Si}}$ (Row 1, Cols 2-4):** 
	Direct coupling to $T_{Si,1}$. Zeros for $T_{Si,2}, T_{Si,3}$. 
3. **$\mathbf{Link_{0 \to TEC}}$ (Row 1, Cols 5-7):** 
	Direct coupling to $T_{c,1}$. Zeros for $T_{c,2}, T_{c,3}$. 
4. **$\mathbf{Link_{Si \to 0}}$ (Rows 2-4, Col 1):** 
	Reciprocal coupling from $T_{Si,1}$ to Center. 
5. **$\mathbf{A_{Silicon}}$ (Rows 2-4, Cols 2-4):** 
	* **Row 2 ($Si,1$):** Includes connection to Center ($R_{Si,0\to1}$) and Node 2. 
	* **Row 3 ($Si,2$):** General internal node (connected to 1 and 3). 
	* **Row 4 ($Si,3$):** Boundary node. Note the absence of $1/R_{lat,3}$ (lateral out) as this is the edge of the chip model. 
6. **$\mathbf{A_{Vert, Si\to c}}$ (Rows 2-4, Cols 5-7):** 
	Perfectly diagonal matrix containing $1/R_{v,1}, 1/R_{v,2}, 1/R_{v,3}$. This couples the chip directly to the TEC underneath it. 
7. **$\mathbf{Link_{TEC \to 0}}$ (Rows 5-7, Col 1):** 
	Reciprocal coupling from $T_{c,1}$ to Center. 
8. **$\mathbf{A_{Vert, c\to Si}}$ (Rows 5-7, Cols 2-4):** 
	Perfectly diagonal matrix containing $1/R_{v,1}, 1/R_{v,2}, 1/R_{v,3}$. 
9. **$\mathbf{A_{TEC}}$ (Rows 5-7, Cols 5-7):** 
	* **Row 5 ($c,1$):** Coldest stage. Connects to Center and pumps to Node 2. 
	* **Row 6 ($c,2$):** Middle stage. Receives heat from Node 1, pumps to Node 3. 
	* **Row 7 ($c,3$):** Hot stage. Receives heat from Node 2, pumps to $T_w$. Note that the term for the water connection ($K_3$) appears in the diagonal self-term, while the fixed temperature portion ($K_3 T_w$) has moved to vector **B**.
