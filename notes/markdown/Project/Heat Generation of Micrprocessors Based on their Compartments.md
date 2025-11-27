Microprocessors are compartmentalized into regions based on their role in the overall system, therefore, the rate of heat generation of these regions vastly differ. This aspect can be exploited to the advantage when building 3D stacked microchips.

[Source](https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.exxactcorp.com%2Fblog%2FHPC%2Fwhat-is-a-cpu-die-shrink-and-what-does-it-mean-for-the-future-&psig=AOvVaw0l7TFawlQvQhxo--zSmkRR&ust=1764241374824000&source=images&cd=vfe&opi=89978449&ved=0CBgQjhxqFwoTCLCDiLjVj5EDFQAAAAAdAAAAABAE)
![[intel-i7-3960x-die-diagram-1.jpg|7 Die Map of a Hexa-Core Coffee Lake Processor | Download Scientific Diagram]]
[Source](https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FDie-Map-of-a-Hexa-Core-Coffee-Lake-Processor_fig6_332543387&psig=AOvVaw0l7TFawlQvQhxo--zSmkRR&ust=1764241374824000&source=images&opi=89978449)
![[Die-Map-of-a-Hexa-Core-Coffee-Lake-Processor.png|7 Die Map of a Hexa-Core Coffee Lake Processor | Download Scientific Diagram]]

# 1 Parts of a Chip

## 1.1 Traditional Micoprocessors

**Key regions in a typical CPU (traditional desktop/server microprocessor):**

- **Cores:** Independent execution units that run program instructions. Each core has its own integer and floating-point pipelines.
    
- **Cache hierarchy:**
    - **L1 Cache:** Very small, very fast (per core).
    - **L2 Cache:** Larger, still per core or shared by small clusters.
    - **L3 Cache:** Large, shared across all cores.
        
- **Memory Controller:** Interfaces with external DRAM.
- **Interconnect / Ring bus / Mesh:** Connects cores, caches, and I/O blocks.
- **Media / SIMD Engines (e.g., AVX, NEON):** Accelerators for vector operations, multimedia, signal

- **NPUs:**
    - Traditional PC/server CPUs do **not** include NPUs.
    - Modern mobile-oriented or hybrid processors (e.g., Apple Silicon, Qualcomm Snapdragon, Intel Core Ultra) **often include NPUs** on the same die.
        
- **DRAM:**
    - Standard CPUs do **not** integrate DRAM. They only integrate the memory controller; DRAM is external.
    - Some SoCs use **stacked DRAM** (HBM) on the package, but not on the CPU die itself.

## 1.2 System on Chip (SOC)

A processor is called a **System-on-Chip (SoC)** when it integrates many components typically found on a full computer motherboard:

- CPU cores
- GPU
- NPU or DSP
- memory controller
- media codecs
- I/O controllers (USB, PCIe, networking, etc.)
- sometimes modem, ISP, audio hardware, sensor hubs

**Regular CPU = mostly cores + cache + memory controller**  
**SoC = many subsystems integrated, often including GPU/NPU and I/O**

**Examples of well-known SoCs:**
- **Apple Silicon** (M1, M2, M3, A-series): CPU + GPU + NPU + ISP + media engine + unified memory on-package
- **Qualcomm Snapdragon** (mobile SoC with CPU, GPU, NPU, ISP, modem).
- **Samsung Exynos**
- **MediaTek Dimensity**
- **NVIDIA Tegra / Jetson**
- **Intel Core Ultra “Meteor Lake / Lunar Lake”** (CPU tile + GPU tile + NPU tile, considered SoC-like).
- **AMD Ryzen “APU” chips** (CPU + GPU + media engines integrated; sometimes called SoCs).

It can be seem that most modern chips such newest generation intel chips are getting more "SOC like" by adding NPUs and other components directly into the chip.

## 1.3 Chiplets / Multi-die packaging

SOCs can be either made from single photomask with same process node (i.e. 3nm) or they can be made from different process nodes and be combined together. Which is called milti-die approach or chiplets. High-end SoCs increasingly use _chiplets_ or _multi-die_ packaging:

**Examples of multi-die SoCs**
- AMD Ryzen & EPYC CPUs (chiplet design)
- Apple M1 Ultra (two dies joined by [[UltraFusion]])
- Nvidia Grace CPU + LPDDR5X on-package
- Intel Meteor Lake (CPU tile + GPU tile + IO tile + SoC tile)

These use **different process nodes** for different chiplets:

- CPU tile on 4 nm
- GPU tile on 5 nm
- IO tile on 6 nm
- etc.

# 2 Configurations of [[3D Stacked Microchips]]

Usually the different layers (dies) are usually separated **by function**, not simply clones of the same layout stacked vertically. That is, each layer is often a different “module”: e.g. logic, cache, memory, I/O, etc

## Why separate functional blocks into different layers

- **Heterogeneous integration**: Different parts of a system (logic, memory, analog, I/O, cache) have different design requirements (transistor density, leakage tolerance, power, etc.). By fabricating them on different dies (often on different process nodes) and then stacking them, you can optimize each die for its function. [HandWiki+2hsienhsinlee.github.io+2](https://handwiki.org/wiki/Engineering%3AThree-dimensional_integrated_circuit?utm_source=chatgpt.com)
    
- **Yield and cost**: If everything were on one huge monolithic die the yield (fraction of usable chips per wafer) would drop as die size increases. By splitting into smaller dies and stacking, you get higher yield. [HandWiki+2hsienhsinlee.github.io+2](https://handwiki.org/wiki/Engineering%3AThree-dimensional_integrated_circuit?utm_source=chatgpt.com)
    
- **Interconnect efficiency**: Vertical stacking reduces the distance signals must travel between functional blocks (e.g. CPU core ↔ cache, CPU ↔ memory), which decreases latency and power compared to long on‑die wires or off‑chip connections. [HandWiki+2IC Online+2](https://handwiki.org/wiki/Engineering%3AThree-dimensional_integrated_circuit?utm_source=chatgpt.com)


[Source]()
![[b3eb8f68b29fca004d06cfda53f200fa.jpg]]

The most common configuration seems to be stacking DRAMs on top of each other and maintaining other modules on the same plane. Ths is currently used in chiplet configurations.

[Source](https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FChiplet-partitioning-concept_fig1_347759902&psig=AOvVaw0FRJGmeItQGI9QsOCl6x4F&ust=1764247443330000&source=images&cd=vfe&opi=89978449&ved=0CBgQjhxqFwoTCODygYPsj5EDFQAAAAAdAAAAABAg)
![[Chiplet-partitioning-concept.png]]




# 3 Heat Fluxes Based on the Regions



| Component                                     | Typical Heat Flux (W/cm²)         | Notes/Examples                                                                                                                                     | Source Citation                                                                           |
| --------------------------------------------- | --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| **CPU/GPU Cores** (e.g., x86, CUDA/Tensor/RT) | 100-300 (hotspots); avg. 150-250  | Intel i7-8700K: 250 W/cm² max core; i9-9900K: >100 W/cm² overclocked. GPU Tensor cores ~200 W/cm² in high-load AI.                                 | [1](https://eps.ieee.org/images/files/HIR_2023/ch20_thermal.pdf)                          |
| **Cache (SRAM L1/L2/L3)**                     | 5-50 (avg. 10-20)                 | On-chip SRAM: 0.75 W/cm² (older); modern uncore cache ~20 W/cm². Thermal-aware designs reduce to <10 W/cm². Hotspots from core proximity +10-15°C. | [2](https://users.ece.northwestern.edu/~memik/papers/micro05.pdf)                         |
| **Media Engines** (GPU video/decode)          | 20-50                             | Fixed-function; lower than compute cores. Implied in GPU whitepapers as ~20% of total power, spread over larger area.                              |                                                                                           |
| **NPUs** (AI accelerators, MAC arrays)        | 50-100 (MAC hotspots); avg. 20-50 | NPU MAC arrays 2-3x uncore flux; overall lower than CPU cores but thermal stress from parallelism.                                                 |                                                                                           |
| **IO Controllers** (PCIe, memory ctrl.)       | 10-50                             | Part of uncore; ~20 W/cm² avg. in multi-core maps. Voltage regulators (IO-related) ~50-100 W/cm² local.                                            | [3](https://www.researchgate.net/publication/346069484_NPU_Thermal_Management)            |
| **SRAM** (standalone/on-chip)                 | 5-30                              | Similar to cache; 15 W/cm² in processor sims, but lower leakage.                                                                                   |                                                                                           |
| **DRAM** (stacked, e.g., HBM/LPDDR)           | 1-45 (avg. 5-20)                  | Low inherent; 0.15 W/cm² base, up to 45 W/cm² in stacks from crosstalk. Bottom layers +5-15°C from cores.                                          | [4](https://infoscience.epfl.ch/bitstreams/423fdd45-e361-492a-b6ae-1e095cc6aff7/download) |

# 4 Extra

- The thing abot [[media engines]] and how modern SOCs directly incoperate coddecs is interesting. Check them out Later.
- [[x265 and H.265]] first is the software encoder for the second one which is the codec. Seems like media engines directly use hardware as the codec.