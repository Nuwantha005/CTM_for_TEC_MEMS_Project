#project
### Description
This project aims to design a system to cool modern stacked microchips. More importantly, the goal is to extract heat from the bottom dies in horizontal direction. Project aims to use a combination of thermoelectric coolers [[TEC]] and [[micro channel heat-sink]] to extract the heat from the bottom dies of the stacked chip.

### Existing Literature
Most proposes in existing literature falls under two main branches. One aims to use [[TFTEC]] to extract heat in vertical direction, toward the top of the chip. method suffers from the fact that heat from the bottom dies has to pass through top dies, needing the top TECs to have higher pumping power as well as maintaining a steep temperature gradient across the chip. Also this puts a maximum limit on how many dies can be stacked.

And the other method is to integrate [[micro channel heat-sink]] directly between the dies and extract heat through fluid. Directly integrating a fluid channel into the dies is problematic and will cause maintenance / reliability issues. Furthermore, this heat sink layer will be thick and it will increase the distance between dies, causing the [[TSV]]s between dies to to become taller.

### Proposed Solution
I propose to use a [[Multi-Stage TEC]] arranged as a  [[Radial TEC]] to extract heat from the middle of the chip (or [[Local Hot-spots]] like cores) and bring heat to the perimeter of the chip. This heat can then be taken away using a [[micro channel heat-sink]] that goes through the perimeter of the TEC. 

![[6574.png]]

![[Pasted image 20250928123936.png]]

Therefore, the TEC array consists of $N$ number of wedges, which is on itself is a multi-stage TEC that flows heat linearly. When all $N$ wedges are combined, it forms a closed circle or a polygon, meaning each wedge takes $360\degree/N$ angular dimensions.

## Objectives

The main objective is to extract heat flux equivalent to modern high end chip with power rating around $\sim400W$ effectively and efficiently, without sacrificing much power for the cooling effect.

[[MEMS Project Plan]]

