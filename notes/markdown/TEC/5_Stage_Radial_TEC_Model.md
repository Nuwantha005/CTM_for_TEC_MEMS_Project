# 5-Stage Radial TEC Model - SolidWorks Implementation

## 1. Overview

This document describes the geometry and parameters for a 5-stage radial thermoelectric cooler (TEC) model for MEMS chip cooling applications. The model follows the parameterization defined in the thermal network analysis.

## 2. Parameter Summary

### 2.1 Global Design Parameters

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| Number of Stages | $N_s$ | 5 | - | Total TEC stages |
| Wedge Angle | $\theta$ | 30 | degrees | Angular width of one TEC wedge |
| Inner Cylinder Radius | $R_{cyl}$ | 1000 | μm | Central heat spreader radius |
| Chip Width | $w_{chip}$ | 10000 | μm | Square chip dimension |
| Radial Expansion Factor | $k_r$ | 1.3 | - | Ratio $L_{i+1}/L_i$ |
| Inter-Stage Insulation Width | $w_{is}$ | 40 | μm | Ceramic insulator width |
| Azimuthal Insulation Width | $w_{az}$ | 20 | μm | Gap between P/N legs |
| TE Thickness | $t$ | 50 | μm | Thermoelectric leg thickness |
| Chip Layer Thickness | $t_{chip}$ | 50 | μm | Silicon die thickness |
| SOI Layer Thickness | $t_{SOI}$ | 100 | μm | SOI layer for TSVs |

### 2.2 Interconnect Parameters (Ratios)

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Interconnect Ratio | $\alpha_{ic}$ | 0.15 | $w_{ic}/L_i$ |
| Outerconnect Ratio | $\alpha_{oc}$ | 0.15 | $w_{oc}/L_i$ |
| Interconnect Angle | $\beta_{ic}$ | 5° | Azimuthal angle |
| Outerconnect Angle | $\beta_{oc}$ | 5° | Azimuthal angle |
| Interconnect Thickness | $t_{ic}$ | 20 | μm |
| Outerconnect Thickness | $t_{oc}$ | 20 | μm |

### 2.3 TSV Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| TSV Radius | $R_{TSV}$ | 10 | μm |
| TSV Pitch | $P_{TSV}$ | 20 | μm |

## 3. Calculated Geometry (5 Stages)

### 3.1 Derived Values

$$R_{base} = \frac{w_{chip}}{\sqrt{2}} = \frac{10000}{1.4142} = 7071.07 \text{ μm}$$

$$L_{total\_active} = (R_{base} - R_{cyl}) - (N_s + 1) \cdot w_{is}$$
$$L_{total\_active} = (7071.07 - 1000) - 6 \times 40 = 5831.07 \text{ μm}$$

$$geo\_sum = \frac{1 - k_r^{N_s}}{1 - k_r} = \frac{1 - 1.3^5}{1 - 1.3} = 9.0431$$

$$L_1 = \frac{L_{total\_active}}{geo\_sum} = \frac{5831.07}{9.0431} = 644.77 \text{ μm}$$

### 3.2 Stage-by-Stage Dimensions

| Stage | $L_i$ (μm) | $r_{start}$ (μm) | $r_{end}$ (μm) | $w_{ic}$ (μm) | $w_{oc}$ (μm) |
|-------|------------|------------------|----------------|---------------|---------------|
| 1 | 644.77 | 1040.00 | 1684.77 | 96.72 | 96.72 |
| 2 | 838.20 | 1724.77 | 2562.97 | 125.73 | 125.73 |
| 3 | 1089.66 | 2602.97 | 3692.63 | 163.45 | 163.45 |
| 4 | 1416.56 | 3732.63 | 5149.19 | 212.48 | 212.48 |
| 5 | 1841.52 | 5189.19 | 7030.71 | 276.23 | 276.23 |

### 3.3 Radii Verification

$$\sum_{i=1}^{5} L_i + R_{cyl} + 6 \cdot w_{is} = 5830.71 + 1000 + 240 = 7070.71 \approx R_{base}$$

## 4. SolidWorks Model Structure

### 4.1 Feature Tree

```
Part1
├── CentralCylinder
│   └── Circle R=1000μm, Extrude 50μm
├── Insulator_Center (w_is = 40μm)
│
├── Stage1
│   ├── Stage1_P_Leg (wedge)
│   ├── Stage1_N_Leg (wedge)
│   ├── Stage1_Interconnect (copper)
│   └── Stage1_Outerconnect (copper)
├── Insulator_S1_S2
│
├── Stage2
│   ├── Stage2_P_Leg
│   ├── Stage2_N_Leg
│   ├── Stage2_Interconnect
│   └── Stage2_Outerconnect
├── Insulator_S2_S3
│
├── Stage3
│   ├── Stage3_P_Leg
│   ├── Stage3_N_Leg
│   ├── Stage3_Interconnect
│   └── Stage3_Outerconnect
├── Insulator_S3_S4
│
├── Stage4
│   ├── Stage4_P_Leg
│   ├── Stage4_N_Leg
│   ├── Stage4_Interconnect
│   └── Stage4_Outerconnect
├── Insulator_S4_S5
│
├── Stage5
│   ├── Stage5_P_Leg
│   ├── Stage5_N_Leg
│   ├── Stage5_Interconnect
│   └── Stage5_Outerconnect
└── Insulator_Outer
```

### 4.2 Sketch Geometry Details

#### Wedge Shape for TEC Legs
Each TEC leg is a wedge-shaped extrusion defined by:
- Inner radius: $r_{start,i} + w_{ic,i}$
- Outer radius: $r_{end,i} - w_{oc,i}$
- Angular extent: $\theta/4 - w_{az}/r_{avg}$ (for each leg)
- Extrusion height: $t$ (50 μm)

#### P-Type Leg Position
- Angular start: $+w_{az}/(2 \cdot r_{avg})$
- Angular end: $+\theta/2 - w_{az}/(2 \cdot r_{avg})$

#### N-Type Leg Position  
- Angular start: $-\theta/2 + w_{az}/(2 \cdot r_{avg})$
- Angular end: $-w_{az}/(2 \cdot r_{avg})$

#### Interconnect (Copper)
- Inner radius: $r_{start,i}$
- Outer radius: $r_{start,i} + w_{ic,i}$
- Angular extent: $\pm\beta_{ic}/2$
- Extrusion height: $t_{ic}$ (20 μm)

#### Outerconnect (Copper)
- Inner radius: $r_{end,i} - w_{oc,i}$
- Outer radius: $r_{end,i}$
- Angular extent: $\pm\beta_{oc}/2$
- Extrusion height: $t_{oc}$ (20 μm)

## 5. Manual Creation Steps in SolidWorks

### Step 1: Create Central Cylinder
1. Select **Top Plane** → Insert Sketch
2. Draw circle at origin with radius = 1 mm (1000 μm)
3. Exit sketch → Boss-Extrude 0.05 mm (50 μm)
4. Rename to "CentralCylinder"

### Step 2: Add Equations
Go to **Tools → Equations** and add:
```
"N_stages" = 5
"theta_deg" = 30
"R_cyl" = 1000um
"k_r" = 1.3
"w_is" = 40um
"w_az" = 20um
"t_TE" = 50um
"alpha_ic" = 0.15
"alpha_oc" = 0.15
"L_1" = 644.77um
"L_2" = 838.20um
"L_3" = 1089.66um
"L_4" = 1416.56um
"L_5" = 1841.52um
```

### Step 3: Create Stage 1
1. **Insulator Ring** (R=1000 to R=1040, θ=±15°)
2. **P-Leg** (R=1136.72 to R=1588.05, θ=+1° to +14°)
3. **N-Leg** (R=1136.72 to R=1588.05, θ=-14° to -1°)
4. **Interconnect** (R=1040 to R=1136.72, θ=±2.5°)
5. **Outerconnect** (R=1588.05 to R=1684.77, θ=±2.5°)

### Step 4: Repeat for Stages 2-5
Use the calculated dimensions from Section 3.2.

### Step 5: Create Circular Pattern
1. Select all stage features
2. **Insert → Pattern/Mirror → Circular Pattern**
3. Axis: Z-axis (vertical through origin)
4. Instances: 360°/30° = 12 copies

## 6. Materials Assignment

| Component | Material | κ (W/m·K) | ρ (Ω·m) |
|-----------|----------|-----------|---------|
| TE Legs (P/N) | Bi₂Te₃ | 1.2 | 1.0×10⁻⁵ |
| Interconnects | Copper | 400 | 1.7×10⁻⁸ |
| Insulators | AlN | 170 | >10¹⁰ |
| Azimuthal Gap | SiO₂ | 1.4 | >10¹⁴ |
| Central Cylinder | Silicon | 150 | - |

## 7. VBA Macro

A complete VBA macro for automated model generation is available at:
`scripts/solidworks/RadialTEC_5Stage.swp.bas`

To use:
1. Open SolidWorks
2. **Tools → Macro → New**
3. Copy the VBA code
4. Run `Main()` subroutine

## 8. Design Table Integration

For parametric studies, create a Design Table with these columns:

| Configuration | k_r | N_stages | theta_deg | w_is | w_az |
|--------------|-----|----------|-----------|------|------|
| Baseline | 1.3 | 5 | 30 | 40 | 20 |
| High_Expansion | 1.5 | 5 | 30 | 40 | 20 |
| Compact | 1.2 | 5 | 45 | 30 | 15 |
| 4_Stage | 1.3 | 4 | 30 | 40 | 20 |

## 9. Export for COMSOL

For FEA analysis, export as:
- **STEP** format for geometry import
- **Parasolid** for better compatibility
- Include material assignments in export

## 10. Verification Checklist

- [ ] All 5 stages created with correct radii
- [ ] P and N legs have azimuthal gap (w_az)
- [ ] Inter-stage insulators at correct positions
- [ ] Interconnects span between stages
- [ ] Circular pattern creates full 360° coverage
- [ ] Materials assigned correctly
- [ ] Equations linked to dimensions
