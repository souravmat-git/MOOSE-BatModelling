# MOOSE-BatModelling
A diffuse interface methodology for modelling interfacial voids in all-solid-state batteries utilizing MOOSE.

In this project, a diffuse interface model is proposed to simulate the concurrent evolution of interfacial voids and the negative electrode during the stripping and plating processes. The model is applied to an all-solid-state sodium (Na) metal battery with homogeneous and polycrystalline solid electrolyte (SE) separators. Three cases are presented. First, the evolution of the Na negative electrode during stripping and plating with perfect electrode/SE separator interface. Second, the evolution of a single interfacial void placed at the electrode/SE separator interface during stripping and plating. Third, the coalescence of multiple interfacial voids along the electrode/SE separator interface during stripping.

The governing equations of the model are implemented using [MOOSE](https://mooseframework.inl.gov/), an open-source finite element framework.
For running the code, the user can build a MOOSE-based [application](https://mooseframework.inl.gov/getting_started/new_users.html) by enabling the phase field and solid mechanics modules. The user can then run the input files within this repository using the following command:

`mpiexec -n (# of CPUs) /dir/to/your/app/yourapp-opt -i /dir/to/your/inputfile/input_file.i`

This repository also contains additional source codes specifically developed for this project. These codes are located within the 'src' folder. It is necessary to compile these codes within your MOOSE-based application before executing the aforementioned command. 
