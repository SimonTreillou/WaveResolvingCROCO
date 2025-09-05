# WaveResolvingCROCO

Some CROCO code and implementations for wave-resolving simulations.

### 01-spectrum_input_file
- Allow you to use a specific spectrum as input via a file text ("freq[i] S[i]" format). Check to have the adequate number of frequency via the provided Python code.

### 02-onlineparticletracking
- Test of the module ```FLOATS``` in a wave-resolving simulation.
- Deprecated now, and removed from the CROCO source code.

### 03-frequencydependantangle
- Extend the wavemaker capacity to have frequency-dependant wave angles.

### 04-crocowithopendrift
- Add first try with OpenDrift. Routines adds lon/lat/angle and masks to files.
- 
### 05-onlinevorticitybudget
- Add first try of online vorticity budget with CROCO.

### To do: 
- [ ] Energy-based frequency decomposition
- [x] Online particle tracking
- [ ] Frequency-dependant angles
- [ ] Spectrum input file: clean
- [ ] Validate infragravity waves?
- [ ] OBC with angle (land masking?)
- [ ] Fix WAVE_MAKER_EAST
- [ ] Fix stations and float when NO_TRACER
- [ ] progressive resolution
- [x] Online vorticity budget
- [ ] Offline vorticity budget?
