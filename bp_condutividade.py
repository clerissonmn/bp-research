#
# bp_codutividade: Script utilizado para validar a condutividade do BP
#
# Referência: 
# Liu, Chao, et al. "Surface plasmon resonances between silver nanoribbons and 
# anisotropic black phosphorus to light confinement." 
# Journal of Physics D: Applied Physics 54.22 (2021): 225202.
#
#%% [ Imports ]--------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constantes # https://docs.scipy.org/doc/scipy/reference/constants.html

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from warnings import warn

# [ Função para a condutividade ]---------------------------------------------

def bp_sigma(range_lda=None, range_freq=None, pot_quimico=None):
    """Calcula a condutividade para uma camada de fosforeno com base em 
    Liu, Chao, et al. Journal of Physics D: Applied Physics 54.22 (2021): 225202.

    Todos os parâmetros foram convertidos para assegurar estarem no SI. O range de 

    de frequências ou de  comprimentos de onda são obrigatórios mas mútuamente exclusivos, 
    ou seja, se forem passados os dois, a prioridade é comprimento de onda. O potencial 
    químico deve ser dado em eV (a conversão para o SI será feita automática quando Soportuno).

    Args:
        range_lda (list, tuple ou escalar, obrigatório): Range de comprimento de ondas em m.Padrão é None.
        range_freq (list, tuple ou escalar): Range de frequências em Hz. Padrão é None.
        pot_quimico (float, optional): Potêncial químico em eV.Padrão é None.

    Returns:
        dict de complexos: dict no formato {'xx':[a+jb,...], 'yy':[a+jb,...]} com as 
                            condutividades AC e ZZ respectivamente.
    """
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
    gamma1 = Joule(4 * a/np.pi)  # Jm
    eta    = Joule(10*1e-3)   # J

    # ----------[ Testes de sanidade ]---------- #
    if range_lda is not None:
        lda = np.asarray(range_lda)
        freq = c/lda
    else:
        if range_freq is not None:
            freq = np.asarray(range_freq)
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

    nc   = (m1D*kBT / (np.pi*hbar**2)) * np.log(1+np.exp(EF_EC/kBT))

    # ----------[ Condutividade ]---------- #
    omega = 2*np.pi*freq
    eta_hbar = eta/hbar
    
    DAC = np.pi * e**2 * nc / mAC
    sigAC = 1j*DAC / (np.pi*(omega + 1j*eta_hbar))

    DZZ = np.pi * e**2 * nc / mZZ
    sigZZ = 1j*DZZ / (np.pi*(omega + 1j*eta_hbar))

    return {'xx':sigAC, 'yy':sigZZ}
    
# [ Plot ]------------------------------------------------------------------

llda   = np.linspace(4.0,20.0,100)*1e-6   # m
mmu = 0.6
sigma = bp_sigma(range_lda=llda, pot_quimico=0.6)

x  = llda
y1 = sigma['xx']
y1 = sigma['yy']

fig_eps = make_subplots()
fig_eps.add_trace(go.Scatter(x=x, y=np.real(y1), name='Re{sig_ACC}' ))
fig_eps.add_trace(go.Scatter(x=x, y=np.imag(y1), name='Imag{sig_AC}'))

fig_eps.update_yaxes(type="log")

fig_eps.show()

# %%
