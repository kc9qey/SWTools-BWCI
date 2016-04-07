module data_module

  USE kind_module, ONLY : double
  

  INTEGER :: &
    ncycle, &
    imin, imax, &
    nproc_y, nproc_z, &
    myid
  REAL(KIND=double) :: &
    time
  REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: &
    x_ef, &
    y_ef, &
    z_ef, &
    dx_cf, &
    dy_cf, &
    dz_cf, &
    x_cf, &
    y_cf, &
    z_cf, &
    e_nu_c_bar, &
    f_nu_e_bar
  REAL(kind=double), allocatable, dimension(:,:) :: &
    r_shock, &
    r_shock_mn, &
    r_shock_mx, &
    tau_adv, &
    tau_heat_nu, &
    tau_heat_nuc, &
    r_nse, &
    e_rad, &
    elec_rad
  REAL(kind=double), allocatable, dimension(:,:,:) :: &
    rho_c, &
    t_c, &
    ye_c, &
    u_c, &
    v_c, &
    w_c, &
    pMD, &
    sMD, &
    dudt_nuc, &
    dudt_nu, &
    grav_x_c, &
    grav_y_c, &
    grav_z_c, &
    rsphere_mean, &
    dsphere_mean, &
    tsphere_mean, &
    msphere_mean, &
    esphere_mean, &
    r_gain
  REAL(kind=double), allocatable, dimension(:,:,:) :: &
    unurad, &
    nnurad, &
    nse_c, &
    a_nuc_rep_c, &
    z_nuc_rep_c, &
    be_nuc_rep_c, &
    uburn_c, &
    duesrc
  REAL(kind=double), allocatable, dimension(:,:,:,:) :: &
    unukrad, &
    unujrad, &
    nnukrad, &
    nnujrad, &
    nu_r, &
    nu_rt, &
    nu_rho, &
    nu_rhot, &
    psi0dat, &
    psi1dat, &
    xn_c
  REAL(kind=double), allocatable, dimension(:,:,:,:,:) :: &
    dnurad, &
    psi0_c
    

end module data_module
