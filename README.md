# WaveResolvingCROCO

Some CROCO code and implementations for wave-resolving simulations.

### 01-spectrum_input_file
- Allow you to use a specific spectrum as input via a file text ("freq[i] S[i]" format). Check to have the adequate number of frequency via the provided Python code. Cppkey: WAVE_MAKER_FILE

### 02-onlineparticletracking
- Test of the module ```FLOATS``` in a wave-resolving simulation.
- Deprecated now, and removed from the CROCO source code. New version soon?

### 03-frequencydependantangle
- Extend the wavemaker capacity to have frequency-dependant wave angles.

### 04-crocowithopendrift
- First try with OpenDrift. Routines adds lon/lat/angle and masks to files.
  
### 05-onlinevorticitybudget
- Add first try of online vorticity budget with CROCO.

### 06-offlinevorticitybudget
- 3D vorticity diagnostic as offline tool.

### 07-modifiedonlinevorticitybudget
- 3D vorticity budget implemented in CROCO (similar to DIAGNOSTICS_UV and DIAGNOSTICS_VRT, but with key DIAGNOSTICS_VRT3D)

### To do: 
- [ ] Energy-based frequency decomposition
- [x] Online particle tracking
- [ ] Frequency-dependant angles
- [ ] Spectrum input file: clean
- [ ] Validate infragravity waves?
- [ ] OBC with angle (land masking?)
- [ ] Add averages vorticity budget
- [ ] Add time decomposition vorticity budget
- [x] Fix WAVE_MAKER_EAST (committed on CROCO dev version)
- [ ] Fix stations and float when NO_TRACER
- [ ] progressive resolution
- [x] Online vorticity budget
- [x] Offline vorticity budget?
