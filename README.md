# MOOSE-BatModelling
A diffuse interface methodology for modelling interfacial voids in all-solid-state batteries utilizing MOOSE.

In this project, a diffuse interface model is proposed to simulate the concurrent evolution of interfacial voids and the negative electrode during the stripping and plating processes. The model is applied to an all-solid-state sodium (Na) metal battery with homogeneous and polycrystalline solid electrolyte (SE) separators. Three cases are presented. First, the evolution of the Na negative electrode during stripping and plating with perfect electrode/SE separator interface. Second, the evolution of a single interfacial void placed at the electrode/SE separator interface during stripping and plating. Third, the coalescence of multiple interfacial voids along the electrode/SE separator interface during stripping.

The governing equations of the model are implemented using [MOOSE](https://mooseframework.inl.gov/), an open-source finite element framework.
For running the code, the user can build a MOOSE-based [application](https://mooseframework.inl.gov/getting_started/new_users.html) by enabling the phase field and solid mechanics modules. Once the executable (your-app-opt) of your MOOSE-based application is generated, the user can run the input files using the following command:

`mpiexec -n <# of CPUs> /dir/to/your/executable/your-app-opt -i /dir/to/your/inputfile/input_file.i`

This repository also contains additional custom objects (like `Materials` and `Functions`) specifically developed for this project. The header and source files of these custom objects are located within the `include` and `src` directory of this repository, respectively. It is necessary to compile these custom objects within your MOOSE-based application before executing the aforementioned command. To compile these objects, the user needs to register the objects using the following macro

`registerMooseObject(<AppName>, <ClassName>)`
