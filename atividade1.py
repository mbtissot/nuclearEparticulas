import numpy as np

def valores(T, Zx, Za, Mx, Ma):
	K = 8.617e-5 #[ev K⁻¹]
	mp = 938 # [MeV]

	Mx = mp*Mx # [MeV]
	Ma = mp*Ma # [MeV]
	mr = Mx*Ma/(Mx+Ma) # [MeV]
	densp = 4.45e31 # [m⁻³]

	Eg = 986.4*(mr/mp)*(Zx*Za)**2 # [keV]
	kt = K*T/1000 # [keV]
	E0 = (Eg*(kt**2)/4)**(1/3) # [keV]
	deltaE = (2*Eg*kt**5)**(1/6) * 3**(-1/2) # [keV]
	nu = -2/3 + (Eg/(4*kt))**(1/3)
	c = 3e8
	R = c*(densp*densp/2) * (2**(13/6)) * (3**(-1/2)) * (3.8e-50) * (Eg**(1/6)) * ((1000*mr)**(-1/2)) * (kt**(-2/3)) * np.exp(-3*(Eg/(4*kt))**(1/3))
	tau = densp/(2*R) # [s]
	return Eg, kt, E0, deltaE, nu, R, tau

print("------ Ciclo ppI ------")
valppi = [1.6e7, 1, 1, 1, 1]
PPI = valores(*valppi)
print(f"""Eg: {PPI[0]:.2f} keV
KT: {PPI[1]:.2f} keV
E0: {PPI[2]:.2f} keV
deltaE: {PPI[3]:.2f} keV
Nu: {PPI[4]:.2f}
Rpp: {PPI[5]}
Tau: {PPI[6]:.0f} seconds = {PPI[6]/(3.154e7):.2e} anos = {PPI[6]/(3.154e16):.2f} bilhões de anos""")


print("\n------ Ciclo CNO ------")
valcno = [2.5e7, 7, 1, 14, 1]
CNO = valores(*valcno)
print(f"""Eg: {CNO[0]:.2f} keV = {CNO[0]/1000:.2f} MeV
KT: {CNO[1]:.2f} keV
E0: {CNO[2]:.2f} keV
deltaE: {CNO[3]:.2f} keV
Nu: {CNO[4]:.2f}""")