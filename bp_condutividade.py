#
# bp_codutividade: Script utilizado para validar a condutividade do BP
#
# Referência: 
# Liu, Chao, et al. "Surface plasmon resonances between silver nanoribbons and 
# anisotropic black phosphorus to light confinement." 
# Journal of Physics D: Applied Physics 54.22 (2021): 225202.
#
#%% [ Imports ]--------------------------------------------------------------
from numpy import pi, log, exp, linspace, real, imag 
import matplotlib.pyplot as plt
import scipy.constants as constantes # https://docs.scipy.org/doc/scipy/reference/constants.html

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from warnings import warn

# [ Função para a condutividade ]---------------------------------------------

def bp_sigma(range_lda=None, range_freq=None, pot_quimico=None):

    # ----------[ Funções ]---------- #
    eV    = lambda x: x/1.602176634e-19 # Converte x J para eV
    Joule = lambda x: x*1.602176634e-19 # Converte x eV para J

    # ----------[ Variáveis e constantes ]---------- #
    hbar = constantes.value(u"reduced Planck constant")     #  Js: 
    m0   = constantes.value(u"electron mass")               #  kg: 9.1093837015e-31
    kBT  = constantes.value(u"Boltzmann constant")*300      #   J: 
    e    = constantes.value(u"elementary charge")           #   C: 1.602176634e-19
    c    = constantes.value(u"speed of light in vacuum")    # m/s: 299792458.0

    Delta  = Joule(2.0)       # J
    a      = 0.223*1e-9       # m
    gamma1 = Joule(4 * a/pi)  # Jm
    eta    = Joule(10*1e-3)   # J

    # ----------[ Testes de sanidade ]---------- #
    if range_lda is not None:
        lda = range_lda
        freq = c/lda
    else:
        if range_freq is not None:
            freq = range_freq
        else:
            raise Exception("Faltando parâmetros. range_freq ou range_lda devem estar definidos.")
        
    if pot_quimico is not None:
        mu_c = pot_quimico
    else:
        mu_c = 0.6
        warn('pot_quimico não está definido. Usando o valor padrão de 0.6eV ')
    # Testa se o usuário inseriu um valor

    # ----------[ Massas ]---------- #
    eta_c  = hbar**2 / (0.4*m0)
    nu_c   = hbar**2 / (1.4*m0)

    mAC    = hbar**2 / (2*gamma1**2 / Delta + eta_c)
    mZZ    = hbar**2 / (2*nu_c)

    # ----------[ Densidade de portadores ]---------- #
    m1D   = (mAC*mZZ)**(0.5)  # kg
    EF_EC = Joule(mu_c)       # J: EF - EC (Potencial qímico)

    nc   = (m1D*kBT / (pi*hbar**2)) * log(1+exp(EF_EC/kBT))

    # ----------[ Condutividade ]---------- #
    omega = 2*pi*freq
    eta_hbar = eta/hbar
    
    DAC = pi * e**2 * nc / mAC
    sigAC = 1j*DAC / (pi*(omega + 1j*eta_hbar))

    DZZ = pi * e**2 * nc / mZZ
    sigZZ = 1j*DZZ / (pi*(omega + 1j*eta_hbar))

    return {'xx':sigAC, 'yy':sigZZ}


# [ Plot ]------------------------------------------------------------------

llda   = linspace(4.0,20.0,100)*1e-6   # m
mmu = 0.6
sigma = bp_sigma(range_lda=llda, range_freq=ffreq, pot_quimico=0.6)

x  = llda
y1 = sigma['xx']
y1 = sigma['yy']

fig_eps = make_subplots()
fig_eps.add_trace(go.Scatter(x=x, y=real(y1), name='Re{sig_ACC}' ))
fig_eps.add_trace(go.Scatter(x=x, y=imag(y1), name='Imag{sig_AC}'))

fig_eps.update_yaxes(type="log")

fig_eps.show()