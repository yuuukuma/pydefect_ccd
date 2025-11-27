pydefect_ccd
===========

pydefect_ccd is a tool to generate input files to calculate configuration 
coordinate diagrams for point defects in semiconductors and insulators.
The supported package is only VASP so far.

The detailed theory is written in the following paper:



Requirements
------------
- Python 3.12 or higher
- nonrad 
- pymatgen
- [pydefect](https://github.com/kumagai-group/pydefect)
- see requirements.txt for others

[vise](https://github.com/kumagai-group/vise) is also recommended to generate 
VASP input files and analyze VASP output files.

- License
-----------------------
This code is licensed under the MIT License.


Workflow
-----------------------------------------
The workflow is depicted below.
Here, I show an example using C-on-N defect in GaN.

1. Create a `ccd_init.json` file from two directories containing pydefect files. 
If the excited state has one more (less) charge state, n-type (p-type) is assumed.
```bash
pydefect_ccd make_ccd_init -u ../unitcell.yaml -pbes ../perfect/perfect_band_edge_state.json -fd ../C_N1_-1 -sd ../C_N1_0 -em ../effective_mass.json
```
The defect in `first_dir` needs to show higher formation energy than that in `second_dir`.
For example, in the above command, the formation energy of C_N1_-1

You can always check the json files using the `pydefect_print` command in pydefect.

2. We next construct the directories for CCD calculations.
```bash
pydefect_ccd make_ccd_dirs --ccd_init ccd_init.json 
```

3. After finishing the VASP calculations in each directory, 
run the following commands to generate the pydefect `calc_results.json`, 
`band_edge_orbital_infos.json`, and `band_edge_states.json` files 
in each directory.
```bash
pydefect_vasp cr -d disp*
pydefect ds -d disp*
pydefect deoi -d disp* 
````

4. We can also evaluate the corrections for the CCD calculations.
The details are written in 
[this paper](https://doi.org/10.1103/PhysRevB.107.L220101).
```bash
pydefect_ccd make_ccd_corrections 
```

5. We then create `single_point_info.json` in each directory, which
summarize the calculation result for each single point, with the following command.
```bash
pydefect_ccd make_single_point_results
```

6.
```bash
pydefect_ccd make_potential_curve_result
```


Citing pydefect_ccd
---------------
If pydefect_ccd has been used in your research, please temporary cite 
the following paper.

Yu Kumagai<br>

Please cite the following papers:
- Theory: 
[Alkauskas, Yan, Van de Walle, PRB (2014).](https://doi.org/10.1103/PhysRevB.90.075202)

- Code:
[Turiansky et al., omput. Phys. Commun. (2021).](https://www.sciencedirect.com/science/article/abs/pii/S0010465521001685)

Contact info
--------------
Yu Kumagai<br>
yukumagai@tohoku.ac.jp<br>

Tohoku University (Japan)
