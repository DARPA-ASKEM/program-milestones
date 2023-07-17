def SEIRHD_Model(y, t, N, beta, r_I_to_R, r_I_to_H, r_E_to_I, r_H_to_R, r_H_to_D, p_I_to_H, p_I_to_R, p_H_to_D, p_H_to_R):
    S, E, I, R, H, D = y

    dSdt = -beta * I * S / N
    dEdt = beta* I * S / N - r_E_to_I * E
    dIdt = r_E_to_I * E - (r_I_to_H * p_I_to_H) * I - (r_I_to_R * p_I_to_R * I)
    dRdt = (r_I_to_R * p_I_to_R * I) + (r_H_to_R * p_H_to_R * H)
    dHdt = (r_I_to_H * p_I_to_H * I) - (r_H_to_D * p_H_to_D * H)  - (r_H_to_R * p_H_to_R * H)
    dDdt = r_H_to_D * p_H_to_D * H
    return dSdt, dEdt, dIdt, dRdt, dHdt, dDdt