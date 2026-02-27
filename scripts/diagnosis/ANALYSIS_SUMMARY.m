%% Summary of MATLAB vs COMSOL Discrepancy Analysis
% Date: December 7, 2025
%
% FINDING: The ThermalNetwork connects TEC stages in SERIES when they
%          should be in PARALLEL. This causes a ~25x overestimate of
%          thermal resistance.
%
% EVIDENCE:
%   - COMSOL (single wedge, I=0): R_th = 85.8 K/W
%   - MATLAB series model: R_th = 2148 K/W (25x higher)
%   - MATLAB parallel model: R_th = 207 K/W (2.4x higher)
%
% PHYSICAL UNDERSTANDING:
%   In a radial multi-stage TEC:
%   - Heat enters through the chip bottom surface (distributed over entire area)
%   - Chip layer spreads heat laterally
%   - Heat flows VERTICALLY through thin insulator into each TEC stage
%   - Each TEC stage conducts heat RADIALLY to its cold side
%   - All cold sides connect to the water/heat sink in PARALLEL
%
%   ASCII diagram:
%
%                        WATER (cold sink)
%                    ______|______|______
%                   |      |      |      |
%                   |  TEC1  TEC2  TEC3  |
%                   |  (r1)  (r2)  (r3)  |
%                   |______|______|______|
%                          |
%                       INSULATOR
%                          |
%                    ______|______
%                   |             |
%                   |    CHIP     | ← Heat flux input
%                   |_____________|
%
%   Heat flows: Chip → Insulator → TEC stages (in parallel) → Water
%
% CURRENT BUG IN ThermalNetwork.m:
%   Lines ~260-265:
%     if i < N
%         M(idx_c, idx_c_next) = K_stages(i);  % WRONG: series connection
%     else
%         B(idx_c) -= K_stages(i) * T_water;   % Only last stage to water
%     end
%
% PROPOSED FIX:
%   Each stage cold side should connect directly to water:
%     B(idx_c) = B(idx_c) - K_stages(i) * T_water;  % ALL stages to water
%
%   The inter-stage connections (via K_stages) should be REMOVED or
%   reinterpreted as lateral connections (if any).
%
% ADDITIONAL PARAMETER MISMATCH FOUND:
%   - MATLAB was using t_TEC = 500 µm
%   - COMSOL uses LL_t_TEC = 1000 µm
%   This has been corrected in the diagnosis scripts.
%
% REMAINING 2.4x DISCREPANCY (after parallel fix):
%   Could be due to:
%   1. Material property differences (k_Bi2Te3)
%   2. 3D heat spreading effects not in 1D model
%   3. Boundary condition differences
%   4. Fill factor interpretation

fprintf('See comments in this file for detailed analysis.\n');
fprintf('Main finding: TEC stages connected in SERIES should be PARALLEL.\n');
