module constants

const g = 9.81             #: Gravitational constant, m/s^2
const Cp = 1004.0          #: Specific heat of dry air at constant pressure, J/kg
const L = 2.25e6           #: Latent heat of condensation, J/kg
const rho_w = 1e3          #: Density of water, kg/m^3
const R = 8.314            #: Universal gas constant, J/(mol K)
const Mw = 18.0/1e3        #: Molecular weight of water, kg/mol
const Ma = 28.9/1e3        #: Molecular weight of dry air, kg/mol
const Rd = R/Ma            #: Gas constant for dry air, J/(kg K)
const Rv = R/Mw            #: Gas constant for water vapor, J/(kg K)
const Dv = 3.e-5           #: Diffusivity of water vapor in air, m^2/s
const ac = 1.0             #: condensation constant
const Ka = 2.e-2           #: Thermal conductivity of air, J/m/s/K
const at = 0.96            #: thermal accomodation coefficient
const epsilon = 0.622      #: molecular weight of water / molecular weight of dry air

# Additional fixed model parameters
const N_STATE_VARS = 7
# const STATE_VARS = ["z", "P", "T", "wv", "wc", "wi", "S"]
# const STATE_VAR_MAP = {var: i for i, var in enumerate(STATE_VARS)}

end # module constants
