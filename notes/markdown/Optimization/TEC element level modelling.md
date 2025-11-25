Applying **Steady-State Energy Balance**, where the net heat flux into any thermal node must be zero:
$$\sum Q_{\text{in}} - \sum Q_{\text{out}} + \sum Q_{\text{generated}} = 0$$
For a single TEC element,
$$\underbrace{Q_{\text{chip}}}_{\text{External Load (Chip Flux)}} + \underbrace{K(T_{h}-T_{c})}_{\text{Back Conduction}}  - \underbrace{\alpha I T_{c}}_{\text{Peltier Cooling}}+ \underbrace{Q_{\text{Joule, deposited at } T_c}}_{\text{Internal Generation}} = 0$$

#### The Joule Heating Component
We need to define $Q_{\text{Joule, deposited at } T_c}$ based on your component-level resistance breakdown:
1. TE Leg Joule Heat: Half of the heat generated in the TE legs: $$Q_{\text{Joule, TE}} = \frac{1}{2}I^2 R_{\text{TE}}$$
2. Interconnect Joule Heat: All the heat generated in the interconnect ($R_{e,ic}$) is deposited at the junction: $$Q_{\text{Joule, Interconnect}} = I^2 R_{e,ic}$$
3. Connecting Wire Joule Heat (if applicable): If the connecting wire runs between two isothermal nodes ($T_c$ and $T_i$), half of its heat is deposited at $T_c$:$$Q_{\text{Joule, Wire}} = \frac{1}{2}I^2 R_{e,w}$$

