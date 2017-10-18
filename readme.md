## PCR efficiency  
Forked from [PCR efficiency Calculator](http://srvgen.upct.es/efficiency.html) by [Izaskun mallona](https://scholar.google.es/citations?user=kMfxCIYAAAAJ&hl=es)

Run testing.py to get a prediction

**Synopsys**

`python testing.py Target_sequence Forward_primer_sequence Reverse_primer_sequence`

**Output**

Number of amplicon sequenced from 1 template at each cycle

**Sidenote**

Geometric progression to find initial number of templates

```
u(n) = q^n * u(0)
u(0) = u(n) / q^n

n : number of PCR cycles
u(n) : number of amplicon at the end of the PCR
u(0) : number of templates at beginning of PCR
q : efficiency
```
