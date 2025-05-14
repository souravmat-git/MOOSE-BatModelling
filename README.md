# MOOSE-BatModelling
A diffuse interface methodology for modelling interfacial voids in all-solid-state batteries utilizing MOOSE.

In this project, a diffuse interface model is proposed to simulate the concurrent evolution of interfacial voids and the negative electrode during the stripping and plating processes. The model is applied to an all-solid-state sodium (Na) metal battery with homogeneous and polycrystalline solid electrolyte (SE) separators. Three cases are presented. First, the evolution of the Na negative electrode during stripping and plating with perfect electrode/SE separator interface. Second, the evolution of a single interfacial void placed at the electrode/SE separator interface during stripping and plating processes. Third, the coalescence of multiple interfacial voids along the electrode/SE separator interface during stripping.

The governing equations of the model are implemented using [MOOSE](https://mooseframework.inl.gov/), an open-source finite element framework.
For running the code, the user can create a MOOSE-based [application](https://mooseframework.inl.gov/getting_started/new_users.html) and enable the phase field and solid mechanics modules. The user can then run the input files within this repository using the following command:

`mpiexec -n (# of CPUs) $EXECUTABLE_DIR/yourapp-opt -i $INPUTFILE_DIR/input_file.i`
