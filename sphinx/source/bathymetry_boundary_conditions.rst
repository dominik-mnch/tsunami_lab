3.3.1
- changes: main(after dambreak()), setup header, setup cpp, supercritical and subcritical get header and cpp setups, bathymetry now works with t_real
- location 3.3.1: 10.01
- Froude number 3.3.1: 0.584356
- location 3.3.2: 10.01
- Froude number 3.3.2: 1.22602
./build/tsunami_lab 1000 25 1 20

3.3 
[nix-shell:~/Project/tsunami_lab]$ ./build/tsunami_lab 1000 25 1 200
Hydraulic jump is where flow switches from Fr > 1 to Fr < 1, so where Froude number crosses 1
Question: should I be looking for the strongest jump? Or at each jump at every timestep?
Jump detected at x = 11.475 at time t = 3.82103
19.9735,11.4625,0.0125,0.0600494,-0.235125,0.124749
19.9735,11.4875,0.0125,0.146527,-0.238781,0.153369
19.9735,11.5125,0.0125,0.208031,-0.2425,0.124748
Die anderen Werte unterscheiden sich nur in der sechsten Nachkommastelle, diese findet sich in solutions 115