#Case: Single void at the electrode/SE separator interface.
#SE separator: Polycrystalline with low-conductivity GBs.

################################################
# Input parameters to setup the mesh
interface_width        = 0.5e-6  #m
characteristic_length  = 1e-9    #m
longitudnal_length     = 50e-6   #m
lateral_length         = 80e-6   #m
       
grid_spacing  = ${fparse interface_width/3.0} #nm

#Energy scale
T = 300      #K
R = 8.314    #J/mol K
R_calper_mole = 1.9872 #cal/(mol K)
F = 96487    #C/mol
kb_eV = 8.62e-5    #eV/K
Vm =${units 23.78 cm^3/mol -> m^3/mol}    # m3/mol
Ec = ${fparse R*T}  #J/mol

#######################################
#Self diffusivity of Na in Na metal
D_Na1 = ${units 8.0 cm^2/s -> m^2/s}
Q_Na1 = 14000  #cal/(mol K)

D_Na2 = ${units 0.144 cm^2/s -> m^2/s}
Q_Na2 = 10840 #cal/(mol K)

D_Na3 = ${units 0.009 cm^2/s -> m^2/s}
Q_Na3 = 8720 #cal/(mol K)

DC_metal = ${fparse D_Na1 * exp(-Q_Na1/(R_calper_mole * T))
                  + D_Na2 * exp(-Q_Na2/(R_calper_mole * T)) 
                  + D_Na3 * exp(-Q_Na3/(R_calper_mole * T))} #m2/s

characteristic_time = ${fparse characteristic_length^2/DC_metal}                  

#Input parameters to calculate kappa and mu
sigma = 0.22  #J/m2

##################################################

real_time = ${units 1 h -> s} #s

###################################################################
# Parameters needed for setting the BC
###################################################################

applied_current_density = 5 #A/m2
#applied_flux = ${fparse applied_current_density/F} #mol/(m2 s)

flux_calc = ${fparse DC_metal/(Vm*characteristic_length)} #mol/(m2 s)
ic_calc   = ${fparse (DC_metal*F)/(Vm*characteristic_length)} #A/m2

################################################################
#Non-dimensional paramters to set up the initial simulation
initial_void_metal = ${fparse -longitudnal_length/characteristic_length/2.0 + 3000}
initial_metal_sei  = ${fparse longitudnal_length/characteristic_length/2.0  - 18500}

[Mesh]

  [gmg]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = ${fparse int(longitudnal_length/grid_spacing)}
    ny   = ${fparse int(lateral_length/grid_spacing)}
    xmin = ${fparse -longitudnal_length/2/characteristic_length}
    xmax = ${fparse  longitudnal_length/2/characteristic_length}
    ymin = ${fparse -lateral_length/2/characteristic_length}
    ymax = ${fparse  lateral_length/2/characteristic_length}
    elem_type = QUAD4
  []

  [separator]
   type = ParsedSubdomainMeshGenerator
   input = gmg
   combinatorial_geometry = 'x >= ${replace initial_metal_sei} & x <= ${replace gmg/xmax}
                            & y >= ${replace gmg/ymin} & y <= ${replace gmg/ymax}'
   block_id = 1
  []

  [iface_metal]
    type = SideSetsBetweenSubdomainsGenerator
    primary_block = 0
    paired_block = 1
    new_boundary = metal_sep_int
    input = separator
  []

[]

[Problem]
  type = FEProblem
  solve = true
[]

[Variables]

  [v_metal]
    block = '0'
  []

  [v_sep]
    block = '1'
  []

  #Order parameter representing the void phase
  [etaa0]
    block = '0'
  []

  #Order parameter representing the metal phase
  [etab0]
    block = '0'
  []

  #Order parameter representing the free space phase
  [eta_space]
    block = '0'
  []

  [w]
    block = '0'
  []

[]

[GlobalParams]
  dy_etas = 'etaa0 etab0 eta_space'
  gamma_names = 'gab gab'
  file_name = 'center_pos.txt'
  op_num = 24
  eta_num = ${op_num}
  grain_num = ${op_num}
  var_name_base = etad
  int_width = ${fparse 1 * interface_width/characteristic_length}
  enable_jit = false
  enable_ad_cache = false
[]

[AuxVariables]

  [xLi]
   order = CONSTANT
   family = MONOMIAL
   block = '0'
  []

  [bnds]
  []

  #Order parameter representing 
  #the separator/electrolyte phase
  [MultiAuxVariables]
   data_type = Real
   variable_base = etad
   grain_num = ${replace GlobalParams/op_num}
   block = 1
   outputs = none
  []

  [ix_metal]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  []

  [iy_metal]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  []

  [ix_sep]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  []

  [iy_sep]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  []

[]

[ICs]

 [ic_eta_space]
   type = FunctionIC
   variable = eta_space
   function = func_eta_space
   block = 0
 []

 [ic_etaa0]
  type = FunctionIC
  variable = etaa0
  function = func_etaa0
  block = 0
 []

 [ic_etab0]
   type = FunctionIC
   variable = etab0
   function = func_etab0
   block = 0
 []

 [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
      block = 1
    []
 []

[]

[Functions]

  [func_eta_space]
    type = ParsedFunction
    expression = 'if(x<${replace initial_void_metal},1,0)'
  []

  [func_etab0]
    type = ParsedFunction
    symbol_values = func_etaa0
    symbol_names = pores
    expression = 'if(x>=${replace initial_void_metal},1-pores,0)'
    #value = '1-pores'
  []

  [func_etaa0]
    type = CircleFunctionIC
    inside_value = 1.0
    outside_value = 0.0
    table_name = 'pore_pos.txt'
  []

[]

[AuxKernels]

  [molefraction_auxkernel]
    type = MaterialRealAux
    variable = xLi
    property = xLi_metal
    block = '0'
    execute_on = 'initial timestep_end'
  []

  [ix_metal]
    type = DiffusionFluxAux
    diffusivity = 'epsilon'
    variable = ix_metal
    diffusion_variable = v_metal
    component = x
    block = 0
  []

  [iy_metal]
    type = DiffusionFluxAux
    diffusivity = 'epsilon'
    variable = iy_metal
    diffusion_variable = v_metal
    component = y
    block = 0
  []

  [ix_sep]
    type = DiffusionFluxAux
    diffusivity = 'epsilon_sep'
    variable = ix_sep
    diffusion_variable = v_sep
    component = x
    block = 1
  []

  [iy_sep]
    type = DiffusionFluxAux
    diffusivity = 'epsilon_sep'
    variable = iy_sep
    diffusion_variable = v_sep
    component = y
    block = 1
  []

  [bnds_sep]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end final'
    block = 1
  []

  [bnds_metal]
    type = BndsCalcAux
    variable = bnds
    v = ${GlobalParams/dy_etas}
    execute_on = 'initial timestep_end final'
    block = 0
  []

[]


[Kernels]

  [potential_metal]
    type = MatDiffusion
    variable = v_metal
    diffusivity = epsilon
    args = ${replace GlobalParams/dy_etas}
    block = 0
  []

  [potential_sep]
    type = MatDiffusion
    variable = v_sep
    diffusivity = epsilon_sep
    block = 1
  []

  ######################################################
  # Phase-field problem will be solved only in block 0
  ######################################################

  #Kernels for order parameter eta_a0
  [ACa0_bulk]
   type = ACGrGrMulti
   variable = etaa0
   v =  'etab0 eta_space'
   mob_name = L
   block = 0
  []

  [ACa0_sw]
   type = ACSwitching
   variable = etaa0
   Fj_names  = 'omegaa omegab omega_space'
   hj_names  = 'ha     hb     h_space'
   coupled_variables = 'etab0 eta_space  w'
   mob_name = L 
   block = 0
  []

  [ACa0_int]
   type = ACInterface
   variable = etaa0
   kappa_name = kappa
   mob_name = L
   variable_L =  false
   block = 0
  []

  [ea0_dot]
    type = TimeDerivative
    variable = etaa0
    block = 0
  []

  # Kernels for order parameter eta_b0 representing the metal region
  [ACb0_bulk]
    type = ACGrGrMulti
    variable = etab0
    v =  'etaa0 eta_space'
    mob_name = L
    block = 0
  []

  [ACb0_sw]
    type = ACSwitching
    variable = etab0
    Fj_names  = 'omegaa omegab omega_space'
    hj_names  = 'ha     hb     h_space'
    coupled_variables = 'etaa0 eta_space w '
    mob_name = L  
    block = 0
  []

  [ACb0_int]
    type = ACInterface
    variable = etab0
    kappa_name = kappa
    mob_name = L
    variable_L =  false
    block = 0
  []

  [eb0_dot]
    type = TimeDerivative
    variable = etab0
    block = 0
  []

  #Kernels for order parameter eta_space
  [ACspace_bulk]
   type = ACGrGrMulti
   variable = eta_space
   v =  'etaa0 etab0'
   mob_name = L
   block = 0
  []

  [ACspace_sw]
   type = ACSwitching
   variable = eta_space
   Fj_names  = 'omegaa omegab omega_space'
   hj_names  = 'ha     hb     h_space'
   coupled_variables = 'etab0 etaa0  w'
   mob_name = L  
   block = 0
  []

  [ACspace_int]
   type = ACInterface
   variable = eta_space
   kappa_name = kappa
   mob_name = L
   variable_L =  false
   block = 0
  []

  [eta_space_dot]
    type = TimeDerivative
    variable = eta_space
    block = 0
  []

  #############################
  # Diffusion problem
  #############################

  [w_dot]
    type = SusceptibilityTimeDerivative
    variable = w
    f_name = chi
    coupled_variables = ${replace GlobalParams/dy_etas}
    block = 0
  []

  [solve_diffusion_metal]
    type = MatDiffusion
    variable = w
    diffusivity = L_BB
    args = ${replace GlobalParams/dy_etas}
    block = 0
  []

  [coupled_etaa0dot]
    type = CoupledSwitchingTimeDerivative
    variable = w
    v = etaa0
    Fj_names = 'rhoa rhob rho_space'
    hj_names = 'ha   hb   h_space'
    coupled_variables = ${replace GlobalParams/dy_etas}
    block = 0
  []

  [coupled_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = w
    v = etab0
    Fj_names = 'rhoa rhob rho_space'
    hj_names = 'ha   hb   h_space'
    coupled_variables  = ${replace GlobalParams/dy_etas}
    block = 0
  []

  [coupled_eta_space_dot]
    type = CoupledSwitchingTimeDerivative
    variable = w
    v = eta_space
    Fj_names = 'rhoa rhob rho_space'
    hj_names = 'ha   hb   h_space'
    coupled_variables  = ${replace GlobalParams/dy_etas}
    block = 0
  []

[]

[InterfaceKernels]

  [equality_of_current_density]
    type = InterfaceDiffusion
    variable = v_metal
    neighbor_var = v_sep
    boundary = 'metal_sep_int'
    D = 'epsilon'
    D_neighbor = 'epsilon_sep'
  []

[]


[BCs]

  ##################################
  #Interfacial boundary conditions #
  ##################################

  [equality_of_electric_potential]
    type = MatchedValueBC
    variable = v_metal
    boundary = 'metal_sep_int'
    v = v_sep
  []

  ####################################
  # Boundary conditions
  ###################################

  [v_left]
    type = DirichletBC
    preset = true
    variable = v_metal
    boundary = left
    value = 0
  []

  [constantcurrent_right]
    type = NeumannBC
    variable = v_sep
    boundary = right
    value =  ${fparse -applied_current_density/ic_calc}
  []  

  [flux_at_interface]
    type = MatNeumannBC
    variable = w
    boundary = 'metal_sep_int'
    boundary_material = 'applied_flux_mat'
    value = '-1'
  []  


[]

[Materials]

  #Parameters to calculate the equilibrium vacancy mole fraction
  vibrational_entropy =  -2.0
  formation_energy_of_vacancy = 1.57e-1 #Sullivan (1974)

  #Parabolic factors dimensionless
  kvoid  = 3.1e2    # void phase
  kspace = ${kvoid} # free space
  kmetal = 1.1e1    # metal phase

  #Equilbrium mole fractions of Na
  equilibrium_vacancy_concentration = ${fparse exp(vibrational_entropy) * exp(-formation_energy_of_vacancy/(kb_eV*T))}

  exLi_void  = 1e-8
  exLi_space = ${exLi_void}
  exLi_metal = ${fparse 1-equilibrium_vacancy_concentration}

  #thermodynamic factors in J/mol
  TF_void  = ${fparse kvoid*Ec}
  TF_space = ${fparse kspace*Ec}
  TF_metal = ${fparse kmetal*Ec}

  #Diffusivities of Li and free space phase
  DC_space = ${fparse DC_metal}
  DC_void  = ${fparse 0}

  #Diffusivities of Na ion in separator phase
  #Parameters to calculate the conductivity of Na+ ion in Na beta alumina ()
  #Grain conductivity from Heinz (2021)                  
  prefactor = ${fparse 8.537e3 * 1e2}     # (S K/m) 
  activation_energy_grain = 0.20 # eV
  cond_sep_measured = ${fparse (prefactor/T) * exp(-activation_energy_grain/(kb_eV*T))} #S/m

  #Properties of the GB
  prefactor_GB = ${fparse 4.786e3 * 1e2}     # (S K/m) 
  activation_energy_GB = 0.35  # eV
  cond_gb_measured = ${fparse  1 * (prefactor_GB/T) * exp(-activation_energy_GB/(kb_eV*T))} #S/m

  ############################################################################
  # We then non-dimensionalize this by defining a characteristic Onsager constant
  # This Onsager constant is the equal to Li_DC_Metal/(R*T) (molm2)/Js
  ##########################################################################

  LB_characteristic  = ${fparse DC_metal/Ec}

  ###########################################################################
  # These are the non-dimensional Onsager mobilities within each phase
  # which is obtained by dividing the diffusivity with thermodynamic factor#
  # The thermodyanic factors for parabolic free energies are the curvature kb*Ec
  ##############################################################################
  LB_void  = ${fparse DC_void/TF_void/LB_characteristic}
  LB_space = ${fparse DC_space/TF_space/LB_characteristic}
  LB_metal = ${fparse DC_metal/TF_metal/LB_characteristic}

  ######################################################################
  #Conductivity is non-dimensionalized
  # by (F^2*D)/(R*T*Vm) which is equal to 0.2756

  #For the void phase, the conductivity was assumed to be identical to
  #that of the metal phase to avoid potential drop at the void/metal interface
  #####################################################################

  factor =  ${fparse 1}

  cond_characteristic = ${fparse (factor * F*F*DC_metal)/(Vm*Ec)}

  #To scale the electronic conductivity on the left, we define another 
  #characteristic conductivity that is much higher than the characteristic ionic conductivity (cond_char)
  #Note that the characteristic current thus is different from "ic_char" as defined in the input file
  #It should be ic_char (anode) =  electronic_cond_char*(R*T)/(F*lc)

  electronic_cond_char = ${fparse cond_characteristic} #S/m

  cond_metal = ${fparse 2.1e7/electronic_cond_char}       #S/m
  cond_sep   = ${fparse cond_sep_measured/cond_characteristic}   #S/m
  cond_void  = ${fparse 2.1e-7/electronic_cond_char}        #S/m
  cond_space = ${cond_metal}
  cond_gb    = ${fparse cond_gb_measured/cond_characteristic} 

  #########################################################################
  # Calculate the interfacial parameters: kappa, m and L
  ########################################################################
  kappa_dim  = ${fparse (3/4) * sigma * interface_width} # in J/m
  mu_dim     = ${fparse  6 *  sigma / interface_width}   # in J/m3
  #zeta      = ${fparse (exLi_metal - exLi_void)^2*TF_metal/(Vm*DC_metal)} # in Js/m5
  #L_dim     = ${fparse  (4 * mu_dim) / (3*kappa_dim*zeta)}   #m3/Js
  L_dim      = ${fparse 1e-4} #m3/Js

  #Reduced values
  kappa_calc = ${fparse kappa_dim / characteristic_length^2 / mu_dim} #dimensionless
  L_calc     = ${fparse L_dim  * characteristic_time * mu_dim}
  
  [interfacial_parameters]
    type = GenericConstantMaterial
    block = 0
    prop_names =  ' kappa          mu      L      gab  '
    prop_values = '${kappa_calc}  1.0    ${L_calc} 1.5 '
  []

  ######################################
  # Material properties for void phase #
  ######################################

  [omegaa]
    type = DerivativeParsedMaterial
    block = 0
    coupled_variables = 'w'
    property_name = omegaa
    expression = '-0.5*w^2/${kvoid} -w*${exLi_void}'
    derivative_order = 2
  []

  [rhoa]
    type = ParsedMaterial
    block = 0
    property_name = rhoa
    material_property_names = 'domegaa_dw:=D[omegaa(w),w]'
    expression = '-domegaa_dw'
  []

  [ha]
    type = SwitchingFunctionMultiPhaseMaterial
    block = 0
    h_name = ha
    phase_etas = 'etaa0'
    all_etas = ${GlobalParams/dy_etas}
  []

  #######################################
  # Material properties for metal phase #
  #######################################

  [omegab]
    type = DerivativeParsedMaterial
    block = 0
    coupled_variables = 'w '
    property_name = omegab
    expression = '-0.5*w^2/${kmetal}-w*${exLi_metal}'
    derivative_order = 2
  []

  [rhob]
    type = ParsedMaterial
    block = 0
    property_name = rhob
    material_property_names = 'domegab_dw:=D[omegab(w),w]'
    expression = '-domegab_dw'
  []

  [hb]
    type = SwitchingFunctionMultiPhaseMaterial
    block = 0
    h_name = hb
    phase_etas = 'etab0'
    all_etas = ${GlobalParams/dy_etas}
  []

  ######################################
  # Material properties for free space #
  ######################################

  [omega_space]
    type = DerivativeParsedMaterial
    block = '0'
    coupled_variables = 'w'
    property_name = omega_space
    expression = '-0.5*w^2/${kspace} -w*${exLi_space}'
    derivative_order = 2
  []

  [rho_space]
    type = ParsedMaterial
    block = '0'
    property_name = rho_space
    material_property_names = 'domega_space_dw:=D[omega_space(w),w]'
    expression = '-domega_space_dw'
  []

  [h_space]
    type = SwitchingFunctionMultiPhaseMaterial
    block = '0'
    h_name = h_space
    phase_etas = 'eta_space'
    all_etas = ${GlobalParams/dy_etas}
  []

  ###################################
  #Define the overall mobilities  ##
  ##################################

  [chi]
    type = DerivativeParsedMaterial
    property_name = chi
    material_property_names = 'ha hb h_space'
    expression = 'ha/${kvoid} + hb/${kmetal} + h_space/${kspace}'
    block = '0'
  []

  #Overall Onsager mobilities
  [interpolated_mobilities_in_left_side]
    type = DerivativeParsedMaterial
    property_name = L_BB
    material_property_names = 'ha hb h_space'
    expression = '${LB_void}*ha + ${LB_metal}*hb + ${LB_space}*h_space'
    block = '0'
    #outputs = exodus
  []

  [epsilon]
    type = ParsedMaterial
    block = 0
    property_name = epsilon
    material_property_names = 'ha hb h_space'
    expression = '${cond_void}*ha + ${cond_metal}*hb + ${cond_space}*h_space'
    #outputs = none
  []

  [epsilon_out_left]
    type = ParsedMaterial
    block = '0'
    property_name = epsilon_out
    material_property_names = 'epsilon'
    expression = 'epsilon * ${cond_characteristic}'
    outputs = exodus
  []

  [epsilon_out_right]
    type = ParsedMaterial
    block = '1'
    property_name = epsilon_out
    material_property_names = 'epsilon_sep'
    expression = 'epsilon_sep * ${cond_characteristic}'
    outputs = exodus
  []
  
  [xLi_metal_left]
    type = ParsedMaterial
    block = 0
    material_property_names = ' ha hb h_space rhoa rhob rho_space'
    expression = 'ha * rhoa + hb * rhob + h_space * rho_space'
    property_name = xLi_metal
    outputs = none
  []

   [xLi_metal_right]
    type = ParsedMaterial
    property_name = xLi_metal
    block = '1'
    expression = 0
  []

  ##############################################
  # Material properties for the separator phase #
  ###############################################

  [hd_sep]
    type = ABSwitchingFunctionMaterial
    h_name = hd_sep
    block = 1
  []

  [hgb_sep]
    type = CTBarrierFunctionMaterial
    function_name = hgb_sep
    g_order = SIMPLE
    block = 1
  []

  [separator_conductivity]
    type = ParsedMaterial
    property_name = epsilon_sep
    material_property_names = 'hd_sep hgb_sep'
    expression = '${cond_sep}*(hd_sep - hgb_sep) + ${cond_gb} * hgb_sep'
    block = '1'
  []

  #In this case, we multiply the non-dimenional current
  # with ic_calc to obtain the dimensional value. Then, we divide 
  # it by F constant to get the dimensional flux. We then divide it by the flux_calc
  # to obtain the net flux. 

  [applied_flux_mat]
    type = ParsedMaterial
    coupled_variables = 'ix_sep etaa0'
    property_name = applied_flux_mat
    expression = 'ix_sep * (1-etaa0) * ${fparse ic_calc/F} * ${fparse 1/flux_calc}'
    boundary = 'metal_sep_int'
    outputs = exodus
  []

[]

[Postprocessors]

  execute_on = 'initial timestep_end final'
  
  [dt]
    type = TimestepSize
    outputs = all
    execute_on = '${execute_on}'
  []

  [computation_time]
   type = PerfGraphData
   section_name = "Root"
   data_type = total
   execute_on = timestep_end
   outputs = all
  []

  [horizontal_pos]
    type = FindValueOnLine
    target = 0.5
    v = etaa0
    start_point = '${fparse -longitudnal_length/2/characteristic_length} 0 0'
    end_point = '${fparse initial_metal_sei - 200} 0 0'
    tol = 1e-8
    execute_on = '${execute_on}'
  []

  [vertical_pos_negative_y]
    type = FindValueOnLine
    target = 0.5
    v = etaa0
    start_point = '${fparse initial_metal_sei - 000} 0 0'
    end_point = '${fparse initial_metal_sei - 000} ${fparse -lateral_length/2/characteristic_length} 0'
    tol = 1e-8
    execute_on = '${execute_on}'
  [] 

  [vertical_pos_positive_y]
    type = FindValueOnLine
    target = 0.5
    v = etaa0
    start_point = '${fparse initial_metal_sei - 000} 0 0'
    end_point = '${fparse initial_metal_sei - 000} ${fparse lateral_length/2/characteristic_length} 0'
    tol = 1e-8
    execute_on = '${execute_on}'
  [] 

  [pp_xNa]
    type = ElementAverageValue
    variable = xLi
    execute_on = ${execute_on}
    outputs = all
    block = '0'
  []

  [applied_flux_at_interface]
    type = SideDiffusiveFluxAverage
    variable = w
    boundary = 'metal_sep_int'
    diffusivity = L_BB
    execute_on = ${execute_on}
    outputs = all
  []

  [applied_current]
    type = SideDiffusiveFluxAverage
    variable = v_sep
    boundary = right
    diffusivity = epsilon_sep
    execute_on = ${execute_on}
    outputs = all
  []

  #This is the applied current in A/m2 calculated 
  #from the applied current using the Side diffusive flux average

  [applied_current_in_dim]
    type = ParsedPostprocessor
    function = 'applied_current * ${fparse ic_calc}'
    pp_names = 'applied_current'
    execute_on = ${execute_on}
    outputs = none
  []

  [applied_current_at_interface_dim]
    type = ParsedPostprocessor
    function = 'applied_flux_at_interface * ${fparse flux_calc * F}'
    pp_names = 'applied_flux_at_interface'
    execute_on = ${execute_on}
    outputs = none
  [] 

[]


#[VectorPostprocessors]

  #[./fields]
  #   type = LineValueSampler
  #   start_point = '${fparse -longitudnal_length/2/characteristic_length} 0 0'
  #   end_point = '${fparse longitudnal_length/2/characteristic_length} 0 0'
  #   variable = 'v_sep v_metal current_density_0'
  #   num_points = ${fparse int(longitudnal_length/grid_spacing) + 1}
  #   sort_by =  id
  #   execute_on =  'final'
  #   outputs = all
  #[../]

#[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 10
    coloring_algorithm = jp
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]



[Executioner]

  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'

  #petsc_options = '-info -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor -snes_monitor -ksp_monitor_true_residual '
  #petsc_options = '-snes_converged_reason'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu mumps'

  #petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm       lu          nonzero'


  l_max_its = 30
  dtmin = 1e-4
  l_tol = 1e-6
  nl_max_its = 25
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1e-9

  start_time = 0
  #num_steps = 30
  end_time = ${fparse real_time/characteristic_time}

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    growth_factor = 1.5
    cutback_factor = 0.8
    iteration_window = 2
    optimal_iterations = 25
  []

[]

[Outputs]

  file_base = ${env OUTPUTFILE_DIR}

  exodus = true
  checkpoint = true

  [./csv]
    type = CSV
    execute_on = 'final' #Dump the VectorPostprocessors at only initial and final
    execute_postprocessors_on = 'initial timestep_end final'
  [../]

  #[./console]
  #  type = Console
  #  output_file = true
  #  outlier_variable_norms = false
  #[../]

[]

[Debug]
  show_var_residual_norms = false
[]
