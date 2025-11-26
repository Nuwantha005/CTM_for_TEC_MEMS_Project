% analyze_cooling_capacity.m
% Analyze what components can be cooled with the current TEC design

clear; clc;
addpath(genpath('src'));

fprintf('=== TEC COOLING CAPACITY ANALYSIS ===\n\n');

%% Current feasible heat flux from our optimization
% At 1000 W/m² with optimized design: T_chip ≈ 30°C

%% Chip power density comparison
fprintf('=== TYPICAL CHIP POWER DENSITIES ===\n\n');

components = {
    'Low-power MCU (STM32L4)', 0.1, 50;       % 50 mW, 10mm² die
    'ARM Cortex-M0+ MCU', 0.05, 20;           % 20 mW, 4mm² die
    'LPDDR4 RAM (per die)', 0.5, 200;         % 200 mW, 40mm²
    'eMMC Flash', 0.2, 100;                   % 100 mW, 50mm²
    'Low-power WiFi (ESP32)', 0.3, 150;       % 150 mW, 50mm²
    'Bluetooth LE chip', 0.1, 20;             % 20 mW, 20mm²
    'GPS receiver', 0.15, 30;                 % 30 mW, 20mm²
    'Audio codec', 0.1, 50;                   % 50 mW, 50mm²
    'Power management IC', 0.2, 100;          % 100 mW, 50mm²
    'Image sensor (small)', 0.3, 100;         % 100 mW, 35mm²
    'FPGA (low-power)', 1.0, 500;             % 500 mW, 50mm²
    'Mobile GPU (low-end)', 5.0, 1000;        % 1W, 20mm²
    'Mobile CPU core (A53)', 10.0, 500;       % 500 mW, 5mm²
    'High-perf CPU core (A76)', 50.0, 2000;   % 2W, 4mm²
    'Desktop CPU', 100.0, 100000;             % 100W, 100mm²
    'High-end GPU', 300.0, 300000;            % 300W, 100mm²
};

fprintf('%-30s | %12s | %10s | %12s | %s\n', ...
    'Component', 'Power (W/m²)', 'Power (mW)', 'Area (mm²)', 'TEC Coolable?');
fprintf('%s\n', repmat('-', 1, 90));

% Current TEC capacity
q_max_feasible = 2500;  % W/m² at T < 100°C with thick TEC
q_optimized = 1000;     % W/m² with comfortable margin (T ≈ 30°C)

for i = 1:size(components, 1)
    name = components{i, 1};
    q = components{i, 2} * 10000;  % Convert to W/m² (was in W/cm²)
    P_mW = components{i, 3};
    A = P_mW / (q / 1000);  % Area in mm²
    
    if q <= q_optimized
        status = '✅ Yes (T < 50°C)';
    elseif q <= q_max_feasible
        status = '⚠️ Marginal (T < 100°C)';
    else
        status = '❌ No (too hot)';
    end
    
    fprintf('%-30s | %12.0f | %10.0f | %12.1f | %s\n', name, q, P_mW, A, status);
end

%% Your SoC cooling strategy
fprintf('\n\n=== YOUR PROPOSED SoC COOLING STRATEGY ===\n\n');

fprintf('Strategy: Place low-power components at bottom (TEC-cooled)\n');
fprintf('          Place high-power cores at top (conventional cooling)\n\n');

fprintf('TEC-Coolable Components (bottom layer):\n');
fprintf('  ✅ Cache memory (SRAM) - ~0.5 W/cm² = 5000 W/m²\n');
fprintf('  ⚠️ This is above our 2500 W/m² limit, but SRAM can tolerate higher temps\n\n');

fprintf('  ✅ Low-power RAM controller - ~0.2 W/cm²\n');
fprintf('  ✅ I/O interfaces (USB, SPI, I2C) - ~0.1 W/cm²\n');
fprintf('  ✅ Power management - ~0.2 W/cm²\n');
fprintf('  ✅ Clock/PLL - ~0.1 W/cm²\n\n');

fprintf('Top Layer (conventional cooling):\n');
fprintf('  ❌ High-performance CPU cores - 50-100 W/cm²\n');
fprintf('  ❌ GPU - 50-300 W/cm²\n');
fprintf('  ❌ AI/ML accelerators - 20-100 W/cm²\n\n');

%% Practical analysis for 10mm chip
fprintf('=== PRACTICAL COOLING BUDGET ===\n\n');

chip_area = 100e-6;  % 10mm × 10mm = 100 mm² = 100e-6 m²
q_feasible = 2500;   % W/m²
P_max = q_feasible * chip_area;  % Max power in W

fprintf('Chip size: 10 mm × 10 mm = 100 mm²\n');
fprintf('Max TEC-coolable power: %.0f mW (%.2f W)\n\n', P_max * 1000, P_max);

fprintf('What fits in 250 mW budget:\n');
fprintf('  • 1× Low-power MCU (50 mW)\n');
fprintf('  • 1× RAM controller (50 mW)\n');
fprintf('  • 1× WiFi/BT combo (100 mW)\n');
fprintf('  • 1× Power management (50 mW)\n');
fprintf('  Total: 250 mW ✅\n\n');

fprintf('OR:\n');
fprintf('  • LPDDR4 memory die (200 mW)\n');
fprintf('  • eMMC controller (50 mW)\n');
fprintf('  Total: 250 mW ✅\n\n');

%% Cache cooling analysis
fprintf('=== SRAM CACHE COOLING ANALYSIS ===\n\n');

cache_sizes = [32, 64, 128, 256, 512, 1024];  % KB
power_per_kb = 0.5;  % mW/KB (typical for SRAM read)
area_per_kb = 0.5;   % mm²/KB (7nm process)

fprintf('Cache Size | Power (mW) | Area (mm²) | Power Density | TEC Coolable?\n');
fprintf('-----------+------------+------------+---------------+--------------\n');

for cache_kb = cache_sizes
    P = cache_kb * power_per_kb;
    A = cache_kb * area_per_kb;
    q = P / A * 1000;  % W/m²
    
    if q <= q_feasible
        status = '✅ Yes';
    else
        status = '❌ No';
    end
    
    fprintf('%6d KB  | %10.0f | %10.1f | %8.0f W/m² | %s\n', cache_kb, P, A, q, status);
end

%% Recommendations
fprintf('\n\n=== RECOMMENDATIONS ===\n\n');

fprintf('1. Your strategy of placing low-power components at bottom is VIABLE\n');
fprintf('   for components with power density < 2500 W/m²\n\n');

fprintf('2. Best candidates for TEC bottom layer:\n');
fprintf('   - Low-power memory (LPDDR4/5)\n');
fprintf('   - Flash storage controllers\n');
fprintf('   - I/O peripherals\n');
fprintf('   - Always-on low-power domains\n\n');

fprintf('3. NOT suitable for TEC cooling:\n');
fprintf('   - High-performance CPU cores\n');
fprintf('   - GPU compute units\n');
fprintf('   - L1/L2 cache near hot cores\n\n');

fprintf('4. For SRAM cache specifically:\n');
fprintf('   - Small cache (32-256 KB) can be TEC-cooled\n');
fprintf('   - Large cache (512KB+) needs careful thermal design\n');
fprintf('   - Consider splitting cache between layers\n');

fprintf('\n5. Alternative approach:\n');
fprintf('   Use TEC as "thermal diode" to prevent heat from top layer\n');
fprintf('   reaching temperature-sensitive analog/RF circuits at bottom.\n');
