#
# bp_codutividade: Script utilizado para validar a condutividade do BP
#
# Referência: 
# Liu, Chao, et al. "Surface plasmon resonances between silver nanoribbons and 
# anisotropic black phosphorus to light confinement." 
# Journal of Physics D: Applied Physics 54.22 (2021): 225202.
#
# %% [ Imports ]--------------------------------------------------------------
from numpy import pi, log, exp, linspace, real, imag 
import matplotlib.pyplot as plt
import scipy.constants as constantes # https://docs.scipy.org/doc/scipy/reference/constants.html

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %% [ Funções ]--------------------------------------------------------------
eV    = lambda x: x*1.602176634e-19 # Converte x de J para eletronvolt
Joule = lambda x: x/1.602176634e-19 # Converte x para de eV para J


# %% [ Variáveis e constantes ]-----------------------------------------------
hbar = constantes.value(u"reduced Planck constant")     #  Js: 
m0   = constantes.value(u"electron mass")               #  kg: 9.1093837015e-31
kBT  = constantes.value(u"Boltzmann constant")*300      #   J: 
e    = constantes.value(u"elementary charge")           #   C: 1.602176634e-19
c    = constantes.value(u"speed of light in vacuum")    # m/s: 299792458.0

# %% [ Cálculo da condutividade ]---------------------------------------------
Delta  = Joule(2.0)            # eV
a      = 0.223*1e-9     # m
gamma1 = 4 * a/pi       # eVm
eta    = 10*1e-3        # eV

# ----------[ Massas ]----------#
eta_c  = hbar**2 / (0.4*m0)
nu_c   = hbar**2 / (1.4*m0)

mAC    = hbar**2 / (2*gamma1**2 / Delta + eta_c)
mZZ    = hbar**2 / (2*nu_c)

# ----------[ nc: Densidade de portadores ]----------#
m1D   = (mAC*mZZ)**(0.5)  # kg
EF_EC = 0.6               # eV: EF - EC
nc   = (m1D*kBT / (pi*hbar**2)) * log(1+exp(EF_EC/kBT))

# ----------[ DAC: Peso de Drude AC ]----------#
DAC = pi * e**2 * nc / mAC

# ----------[ DZZ: Peso de Drude ZZ ]----------#
DZZ = pi * e**2 * nc / mZZ

# ----------[ freq e lda ]----------#
lda   = linspace(4.0,20.0,100)*1e-6   # m
freq  = c/lda                         # Hz
omega = 2*pi*freq

# ----------[ sigmaAC: Condutividade AC ]----------#
eta_hbar = eta/hbar
sigAC = 1j*DAC / (pi*(omega + 1j*eta_hbar))

# ----------[ sigmaZZ: Condutividade ZZ ]----------#
sigZZ = 1j*DZZ / (pi*(omega + 1j*eta_hbar))

# %% [ Plot ]-----------------------------------------------------------------
fig_eps = make_subplots()
fig_eps.add_trace(go.Scatter(x=lda, y=real(sigAC), name='Re{sigACC}' ))
fig_eps.add_trace(go.Scatter(x=lda, y=imag(sigAC), name='Imag{sigAC}'))

fig_eps.update_yaxes(type="log")

fig_eps.show()