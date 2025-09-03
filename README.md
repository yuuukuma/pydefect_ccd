pydefect_ccd
===========

Requirements
------------
- Python 3.10 or higher
- [pydefect](https://github.com/kumagai-group/pydefect)
- see requirements.txt for others

License
-----------------------
Python code is licensed under the MIT License.


Workflow
-----------------------------------------
The workflow is depicted above.

1. Create a `ccd_init.json` file from two directories containing pydefect files. 
If the excited state has one more (less) charge state, n-type (p-type) is assumed.

2. Make directories to calculate configuration coordinate diagrams for ground and excited states.

3. Update `single_point_info.json` in each directory.  
Before running this command, make sure that the `calc_results.json`, 
`band_edge_orbital_infos.json`,  and `band_edge_states.json` files 
have been created using pydefect.

4. Create `single_point_info.json` in each directory.  
Before running this command, make sure that `calc_results.json` and  
`band_edge_states.json` files have been created using pydefect.

Development notes
-------------------
### Bugs, requests and questions
Please use the [Issue Tracker](https://github.com/kumagai-group/pydefect_2d/issues) to report bugs, request features.

### Code contributions
Although pydefect_2d is free to use, we sincerely appreciate if you help us to improve this code.
The simplest but most valuable contribution is to send the feature requests and bug reports.

Please report any bugs and issues at PyDefect's [GitHub Issues page](https://github.com/kumagai-group/pydefect_2d).
Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style follows [PEP8](http://www.python.org/dev/peps/pep-0008) and [Google's writing style](https://google.github.io/styleguide/pyguide.html).
- Add unittests wherever possible including scripts for command line interfaces.

### Tests
Run the tests using `pytest pydefect_2d`.

Citing pydefect_2d
---------------
If pydefect_2d has been used in your research, please cite the following paper.

[Corrections on Formation Energies and Eigenvalues of Point Defect Calculations in Two-Dimensional Materials]<br>
Yu Kumagai<br>


Contact info
--------------
Yu Kumagai<br>
yukumagai@tohoku.ac.jp<br>

Tohoku University (Japan)
