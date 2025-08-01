def compute_depth_averaged_vorticity_budget(u, v, w, z_r, pm, pn, h, zeta):
    # Constants
    dx = 1 / pm  # shape: (eta_rho, xi_rho)
    dy = 1 / pn

    # Step 1: Interpolate u and v to rho points (use CROCO grid tools)
    u_rho = u2rho_3d(u)
    v_rho = v2rho_3d(v)

    # Step 2: Compute zeta (vertical vorticity) at each level
    dzeta_dx = d_dx(v_rho, dx)  # central difference
    dzeta_dy = d_dy(u_rho, dy)
    zeta_3d = dzeta_dx - dzeta_dy

    # Step 3: Compute horizontal gradients of zeta (for advection term)
    dzetadx = d_dx(zeta_3d, dx)
    dzetady = d_dy(zeta_3d, dy)
    dzetadz = d_dz(zeta_3d, z_r)

    # Interpolate w to rho grid
    w_rho = w2rho_3d(w)

    # Advection
    adv = -(u_rho * dzetadx + v_rho * dzetady + w_rho * dzetadz)

    # Stretching term
    div = divergence(u, v, w, dx, dy, z_r)
    stretching = zeta_3d * div

    # Tilting
    dw_dx = d_dx(w_rho, dx)
    dw_dy = d_dy(w_rho, dy)
    dv_dz = d_dz(v_rho, z_r)
    du_dz = d_dz(u_rho, z_r)
    tilting = dw_dx * dv_dz - dw_dy * du_dz

    # Step 4: Integrate in z (depth average)
    dz = np.diff(z_r, axis=0)  # vertical thicknesses
    H = h + zeta

    zeta_dt = depth_average(zeta_3d, dz, H)
    adv_avg = depth_average(adv, dz, H)
    stretching_avg = depth_average(stretching, dz, H)
    tilting_avg = depth_average(tilting, dz, H)

    return zeta_dt, adv_avg, stretching_avg, tilting_avg

