DukeCEvNS

CEvNS recoil response version: 1.0

SNOwGLoBES version: 1.2

#################################################################################################################################################################################

## CEvNS RECOIL RESPONSE INSTRUCTIONS - FOR USE IN SNOWGLOBES

################################################################
### Content:

This subdirectory is home to all the DukeCEvNS code, some of which SNOwGLoBES needs 
to generate a smearing matrix for detectors utilizing Coherent Elastic Neutrino Nucleus 
Scattering (CEvNS) to observe supernova neutrinos. After SNOwGLoBES is made, this subdirectory
should contain an executable cevns_recoil_response.


################################################################
### How the CEvNS code works:

cevns_recoil_response.cc requires 8 arguments:

1: Detector material (Ar, Ge, or Xe, natural abundance only)

2: Form factor (Helm, Klein-Nystrand, Horowitz)

3: Incident neutrino energy, low bound (float, MeV)

4: Incident neutrino energy, high bound (float, MeV)

5: Incident neutrino energy, number of bins (int, dimensionless)

6: Detected energy, low bound (float, MeV)

7: Detected energy, high bound (float, MeV) 

8: Detected energy, number of bins (int, dimensionless)

Cevns_recoil_response.cc uses these arguments to create a whitespace-delimited matrix with M = detected energies 
rows  and N = incident neutrino energies columns. For each entry, cevns_recoil_response.cc calculates the 
differential cross section of interaction for the corresponding incident neutrino energy and detected
energy. The calculated cross-sections are saved as a whitespace-delimited MxN matrix in the subdirectory
/recoil_response_matrices/ with a name corresponding to the input arguments.

The user CAN interact directly with cevns_recoil_response.cc - it is operated using the command

        ./cevns_recoil_response [arguments]

but does not NEED to unless they are trying to diagnose a failure.


################################################################
### How to use the CEvNS code:

Ideally, the user should never need to interact directly with the CEvNS code. Everything is handled 
by the script create_smearing_matrix.py one directory up. When create_smearing_matrix.py (requires Python
3.6+ environment) is directed to create a CEvNS smearing matrix, it gathers the 8 arguments
cevns_recoil_response.cc needs from the user.

Armed with the arguments, create_smearing_matrix.py first searches /dukecevns/recoil_response_matrices/ for 
a file with a name corresponding to the required arguments. If it finds one, it will read it in and 
move on. If it does not find one, it will use the subprocess module and direct cevns_recoil_response.cc to create it. 
Once the required file is found, create_smearing_matrix.py reads it in as a numpy array. 

The array is first normalized by column, turning it into an unsmeared probability matrix - for a 
given (incident neutrino energy) column, each entry represents the probability the corresponding 
detected energy will be observed.

Once normalized, the array is smeared according to the user-selected parameters, formatted as
SNOwGLoBES is expecting, and saved to a .dat file in the subdirectory /snowglobes/smear/new.


################################################################
### Modifying the isotope abundances:

cevns_recoil_response.cc expects, by default, the natural abundances of Xe, Ge, and Ar. It is possible to
change this default by modifying the mixtures.h file.

Open /dukecevns/mixtures.h and locate the element whose abundance you wish to change. Let us take Argon as
an example - on lines 112 and 113, you will see

  isotopes["Ar"].push_back("Ar36");
  molar_fraction["Ar"].push_back(0.003336);

where "Ar36" is the isotope and 0.003336 is its molar fraction. Edit the molar fraction number to reflect
the abundance in your detector. Repeat the process for the other isotopes. 

BE CAREFUL! The molar fractions MUST sum to unity or cevns_recoil_response.cc will return an incorrect answer.


################################################################
### About the Makefile in this subdirectory:

The user should not need this in the course of normal operations. The SNOwGLoBES main Makefile, 
located in /snowglobes/src/, will build all files required by the CEvNS code when SNOwGLoBES is 
compiled. We provide the Makefile here in case the user needs to modify or debug the CEvNS code.

The CEvNS code can be built independently of the rest of SNOwGLoBES using the command

        make cevns_recoil_response

in the /dukecevns/ subdirectory. Also available to the user is the command

        make clean

if cevns_recoil_response.cc is modified.


################################################################
### Future work:

Currently the CEvNS code supports the Helm, Klein-Nystrand, and Hoferichter (numerical) form factors.
Also available as of 2020 is a form factor based on Chiral Effective Field Theory by Hoferichter et al. 
- https://arxiv.org/abs/2007.08529

It is of some interest to compare the Helm and Hoferichter form factors in terms of predicted event 
rate. A user interested in making this comparison should modify cevns_recoil_response.cc and the associated
code to support this new form factor - a place on line 312 of cevns_recoil_response.cc has been designated
for this change.

Given the complexity of the chiral EFT form factor, we suggest contacting Martin Hoferichter and
collaborators and requesting the form factor as a piece of code rather than a mathematical expression.



#################################################################################################################################################################################
#################################################################################################################################################################################

## GENERAL INSTRUCTIONS:
# dukecevns
Simple code for sharing within Duke group for CEvNS rate checks

# Pre-requisition

- https://github.com/nlohmann/json.git
- https://code.ornl.gov/COHERENT/COHERENTProposal2018/tree/master/assumptions (private)

One can run the included `bootstrap.sh` to grap the required packages automatically:

```sh
./bootstrap.sh
```

One can also do it manually:

First, install json:
```sh
cd /path/to/dukecevns
git clone https://github.com/nlohmann/json.git
```

Secondly, grab the following files from ORNL GitLab and put them in `/path/to/dukecevns`:

- `get_flavor_weight.cc`
- `sns_out_BERT_convolved.root` and `sns_out_BERT.root`
- everything in `eff/`, `gs/`, `qf/` and `jsonfiles/`

# Compilation

```sh
cd /path/to/dukecevns
make sns_rates
```

This is done automatically if one uses `bootstrap.sh`.

# Usage

```sh
cd /path/to/dukecevns
./sns_rates base-name-of-a-file-in-the-folder-named-jsonfiles
```

Example `json` files can be found in ORNL GitLab repository `COHERENTProposal2018/assumptions/jsonfiles`.

The output files are located in a new folder named `out`.

## Experimental setup

```json
"timewindow": {
  "start": 1400.0,
  "end": 7400.0
},
```

Kate recommends using the convolved flux histogram. The convolved flux histogram is shifted in time; nominal window to use is 1400 to 7400, although that may get optimized if one knows something about the background.

## Detector response

```json
"detectorresponse": {
  "efftype": "qc", <- qc: charge (q) collected, it can also be # of PEs, ADCs, etc.
  "efficiencyfile": "none", <- file in eff/ folder
  "stepthresh": 2.0, <- lower energy threshold
  "upperthresh": 6.0, <- upper limit of energy
  "qftype": "poly",
  "qfname": "nai",
  "qcperkeVee": 30.0, <- ionization (or light) yield
  "qcbinning": 1, <- bin width
  "gsname": "none", <- electron-equivalent energy resolution (file in gs/)
  "qcsmearing": "poisson", <- "qcgsname", "gsname" are not needed if this is set
  "qcgsname": "none" <- energy resolution in terms of Qc
},
```

### gs file format

`polysqrt` means that it's a polynomial in terms of the square root of energy. The second line defines the range, and the third line includes the coefficients of the formula defined here: https://coherent.phy.duke.edu/wiki/Assumptions#NaI_Smearing (private).

If `efftype = "qc"`, `qcgsname` or `qcsmearing` should be used. If `efftype="eee"`, `gsname` should be used. They cannot be mixed.

### efficiency file format

Column one: energy, # of PEs, or Qc, etc.
Column two: efficiency

### qf file format

Row one: energy range
Row two: quenching factor

If `qftype==poly`, the second row contains polynomial coefficients: C0, C1, C2, etc. so that QF = C0 + C1*E + C2*E^2.

# Output

## sns_diff_rates_quenched-alliso-*.out
Column one: Energy [MeVee]
Column two: Number of events per MeVee



