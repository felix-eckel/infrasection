# infrasection
plot global infrasound sections and perform semblance analysis

## Scripts

### section.py
Plots a global section of waveforms either as raw waveform, root-mean-square amplitudes or STA/LTA characteristic function.

### semblance.py
Plots a semblance analysis for various velocities. Shifts traces by distance and velocity from a given target and calculates a semblance out of the stacked traces.

### semblance_map.py
Similar to semblance.py but instead of searching for velocities from a given origin it searches for the origin itself with a given velocity.

## Usage
define parameters in infrasection.cfg. Run __main__ on the scripts.