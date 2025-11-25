#### 1) [Planar-Radial Structured Thermoelectric Cooler for](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Papers\imeche\Planar-Radial Structured Thermoelectric Cooler for.pdf>)

This paper used 5 wedges to form a thin radial TEC and it showed mere $\sim2.4^{\circ}C$ difference at the junction with modest current density.

Not Much Details are Provided here.

---
## Varying number of elements per stage

### (a) **Variable number of elements per ring**

- You increase the **number of thermoelectric couples (Nₑ)** as the radius grows (so the fill density stays roughly uniform).
- Wiring becomes more complex because you have to series-connect more elements per stage or handle parallel groups.

**Pros:**

- Fine-grained control over local heat pumping capacity.
- Can keep same element size throughout.

**Cons:**

- Wiring and contact resistance can get messy.
- Harder to simulate in axisymmetric COMSOL model (since symmetry breaks).

---

### (b) **Constant number of elements, variable angular width (Δθ)**

- You keep the same _number of annular elements_ all the way through.
- Each one grows in area as ( r ) increases, maintaining continuous annular “slices”.

So each element subtends an angle $\Delta\theta = 2\pi / N_\text{elem}$  — constant across stages.

That means for each stage:  
$$
A_i(r) = w_i(r) \times h_i  
$$
with $w_i(r) = r_i \Delta\theta$ .

**Pros:**

- Perfect for axisymmetric simulation (just model one annular section).
- Simple wiring, symmetric heat spreading.
- Scales naturally with radius — physically meaningful continuity.

**Cons:**

- You lose independent control over packing density; outer region elements are larger (may reduce local current density).

---
