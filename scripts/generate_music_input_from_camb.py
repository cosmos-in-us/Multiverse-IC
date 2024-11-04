import os
import json
import numpy as np
import camb
from camb import model

def initialize_camb_params(cosmo_params, z):
    h = cosmo_params['h']
    Omega_b = cosmo_params['Omega_b']
    Omega_c = cosmo_params['Omega_c']
    Omega_k = cosmo_params['Omega_k']
    w0 = cosmo_params['w0']
    wa = cosmo_params['wa']
    As = cosmo_params['As']
    ns = cosmo_params['ns']
    kmax = cosmo_params['kmax']
    k_per_logint = cosmo_params['k_per_logint']
    
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=100*h, ombh2=Omega_b*h*h, omch2=Omega_c*h*h, omk=Omega_k)
    pars.DarkEnergy = camb.dark_energy.DarkEnergyPPF(w=w0, wa=wa)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_matter_power(redshifts=z, kmax=kmax, k_per_logint=k_per_logint)
    pars.NonLinear = model.NonLinear_none
    
    return pars

def adjust_As_at_z0(pars, cosmo_params):
    z0 = pars.Transfer.PK_redshifts[-1]  # CAMB sorts redshifts (earliest first)
    if z0 != 0:
        print("Redshifts should include 0 here, as we normalize sigma8 at z=0.")
        return None

    results = camb.get_results(pars)
    sigma8_calc = results.get_sigma8()[-1]  # sigma8 at z=0
    sigma8_goal = cosmo_params['sigma8_z0_WMAP5']  # target sigma8 at z=0

    print(f"Target sigma8 = {sigma8_goal} at z = {z0}")
    print(f"Initial sigma8 = {sigma8_calc}, As = {pars.InitPower.As} (trial)")

    As = pars.InitPower.As * (sigma8_goal / sigma8_calc)**2  # new normalization
    pars.InitPower.set_params(As=As, ns=pars.InitPower.ns)  # set new parameterization

    results = camb.get_results(pars)
    sigma8 = results.get_sigma8()[-1]  # new sigma8 at z=0
    print(f"Normalized sigma8 = {sigma8}, diff={sigma8 - sigma8_goal}, As = {pars.InitPower.As} (new)")
    
    return pars

def main():
    # Load cosmological parameters from a JSON file
    with open('cosmo_params.json', 'r') as f:
        cosmo_params = json.load(f)
    
    # Calculate derived parameters
    cosmo_params['sigma8_z0_WMAP5'] = 1 / cosmo_params['bias8_z0_WMAP5']
    cosmo_params['Omega_c'] = cosmo_params['Omega_m'] - cosmo_params['Omega_b']
    cosmo_params['Omega_l'] = 1 - cosmo_params['Omega_m']
    
    z = [cosmo_params['z_start'], 150, 100, 20, 10, 1, 0]
    z.sort(reverse=True)

    pars    = initialize_camb_params(cosmo_params, z=z) # initialize
    pars    = adjust_As_at_z0(pars, cosmo_params) # normalization
    cosmo_params['As'] = pars.InitPower.As  # update the As value
    
    # Generate results
    results = camb.get_results(pars)
    transfers = results.get_matter_transfer_data()
    T = transfers.transfer_data  # Shape: [var, k, z]

    # Prepare output file name
    iz = z.index(cosmo_params['z_start'])
    plus_w0 = '+' if cosmo_params['w0'] >= 0 else ''
    plus_wa = '+' if cosmo_params['wa'] >= 0 else ''
    fname = f"camb_transfer_z{z[iz]:03.0f}_cpl{plus_w0}{cosmo_params['w0']:.1f}{plus_wa}{cosmo_params['wa']:.1f}_lmin{cosmo_params['levelmin']:02d}.txt"
    
    # Save transfer function data
    header = "\n".join([f"{key}: {value}" for key, value in cosmo_params.items()])
    np.savetxt(f"./{fname}", T[:, :, iz].T, fmt="%.10e", header=header)
    
    print(f"MUSIC input data created: ./{fname}")
    print("The number of k values =", T.shape[1])

if __name__ == "__main__":
    main()
