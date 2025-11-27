
# 1 **Concept of 3D Stacking**

In 3D stacking, you take multiple active layers (dies) and stack them vertically rather than just laying everything out in one plane. This allows:

- **Shorter interconnects**: Signals between layers can travel vertically (through-silicon vias, TSVs), which can reduce latency and power consumption.
- **Heterogeneous integration**: You can mix different types of chips (logic, memory, analog, sensors) that may be fabricated using different processes.
- **Higher density**: More functional units per footprint area.

# 2 **Fabrication Approaches**

## 2.1 **Through-Silicon Vias (TSVs)**

- TSVs are vertical electrical connections that pass through the silicon substrate.
- Layers are fabricated separately and then bonded together (wafer-to-wafer or die-to-wafer).
- Example: HBM (High Bandwidth Memory) uses TSVs to stack multiple DRAM dies on top of an interposer.

## 2.2 **Wafer Bonding**

- Entire wafers or individual dies are bonded together using techniques such as:
    - **Direct silicon bonding**: Silicon surfaces are bonded without intermediate layers.
    - **Hybrid bonding**: Combines metal-to-metal and oxide-to-oxide bonding.
- After bonding, TSVs or micro-bumps connect the layers electrically.

## 2.3 **Heterogeneous Node Stacking**

- Often, different layers use **different process nodes**:
    - Logic layers (CPU/GPU cores) → advanced nodes (5 nm, 3 nm).
    - Memory layers (DRAM, SRAM, flash) → more mature, larger nodes (20–40 nm).
    - Analog/RF layers → specialized nodes.
        
- This allows each die to be optimized for its function instead of forcing everything to the same node.

# 3 Example

| **Company / Product**                  | **Type of 3D Stack**                                     | **Notes**                                                                         |
| -------------------------------------- | -------------------------------------------------------- | --------------------------------------------------------------------------------- |
| Intel Foveros (Lakefield, Meteor Lake) | 3D die stacking with logic-on-logic                      | CPU cores at 10 nm, IO and graphics at different nodes, connected via micro-bumps |
| AMD 3D V-Cache                         | SRAM stacked on CPU cores                                | Extra cache dies stacked on top of CPU dies, using TSVs                           |
| HBM (SK Hynix, Samsung)                | DRAM stacking                                            | Multiple DRAM layers bonded, connected via TSVs for high bandwidth                |
| Apple M1/M2                            | Not fully 3D, but interposer + package-level integration | Uses silicon interposer for GPU/DRAM connectivity (heterogeneous integration)     |


![[2021-08-22-image-5-j_1100-1.webp]]


