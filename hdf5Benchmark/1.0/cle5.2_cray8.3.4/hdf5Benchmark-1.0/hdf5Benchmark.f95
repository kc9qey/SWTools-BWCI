program hdf5Test

  use io_module, only: &
    model_write_hdf5, io_walltime
  use data_module
  
  use mpi
  
  integer :: &
    nproc, &
    mpiError, &
    my_j_ray_dim, &
    my_k_ray_dim, &
    j_ray_min, &
    j_ray_max, &
    k_ray_min, &
    k_ray_max, &
    iCycle
  integer, parameter :: &
    !nx = 300000, &
    nx = 30, &
    ny = 8, &
    nz = 8, &
    nnu = 2, &
    nez = 2, &
    nnc = 2, &
    endCycle = 10
  
  nproc_y = 12
  nproc_z = 1
  ncycle  = 1
  imin    = 1
  imax    = nx
  time    = 0.0
  
  call MPI_INIT(mpiError)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, mpiError)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpiError)
  
  if(nproc_y*nproc_z /= nproc)then
    print*, nproc_y, '*', n_proc_z, '/=',  n_proc
    call MPI_ABORT(MPI_COMM_WORLD, mpiError, mpiError)
  end if
  
  my_j_ray_dim = ny/nproc_y
  my_k_ray_dim = nz/nproc_z
      
  j_ray_min = MOD(myid, nproc_y) * my_j_ray_dim + 1
  k_ray_min = (myid/nproc_y) * my_k_ray_dim + 1
  j_ray_max = MOD(myid, nproc_y) * my_j_ray_dim + my_j_ray_dim
  k_ray_max = (myid/nproc_y) * my_k_ray_dim + my_k_ray_dim
  
  !print*, myid, my_j_ray_dim, my_k_ray_dim
  !print*, myid, j_ray_min, j_ray_max, k_ray_min, k_ray_max
  
  allocate(x_ef(nx+1))
  allocate(y_ef(ny+1))
  allocate(z_ef(nz+1))
  allocate(dx_cf(nx))
  allocate(dy_cf(ny))
  allocate(dz_cf(nz))
  allocate(x_cf(nx))
  allocate(y_cf(ny))
  allocate(z_cf(nz))
  allocate(f_nu_e_bar(nx))
  allocate(e_nu_c_bar(nx))
  
  allocate(r_shock(ny,nz))
  allocate(r_shock_mn(ny,nz))
  allocate(r_shock_mx(ny,nz))
  allocate(tau_adv(ny,nz))
  allocate(tau_heat_nu(ny,nz))
  allocate(tau_heat_nuc(ny,nz))
  allocate(r_nse(ny,nz))
  
  allocate(rho_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(t_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(ye_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(u_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(v_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(w_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(pMD(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(sMD(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(dudt_nuc(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(dudt_nu(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(grav_x_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(grav_y_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(grav_z_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(rsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(dsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(tsphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(msphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(esphere_mean(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(r_gain(nnu+1,my_j_ray_dim,my_k_ray_dim))  
  
  allocate(e_rad(my_j_ray_dim,my_k_ray_dim))
  allocate(elec_rad(my_j_ray_dim,my_k_ray_dim))
  allocate(unurad(nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nnurad(nnu,my_j_ray_dim,my_k_ray_dim))
  
  allocate(nse_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(a_nuc_rep_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(z_nuc_rep_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(be_nuc_rep_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(uburn_c(nx,my_j_ray_dim,my_k_ray_dim))
  allocate(duesrc(nx,my_j_ray_dim,my_k_ray_dim))  
  
  allocate(unukrad(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(unujrad(nx,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nnukrad(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nnujrad(nx,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nu_r(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nu_rt(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nu_rho(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(nu_rhot(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(psi0dat(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(psi1dat(nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(xn_c(nx,nnc,my_j_ray_dim,my_k_ray_dim))
    
  allocate(dnurad(nx,nez,nnu,my_j_ray_dim,my_k_ray_dim))
  allocate(psi0_c(nx,nez,nnu,my_j_ray_dim,my_k_ray_dim))

  call random_number(x_ef)
  call random_number(y_ef)
  call random_number(z_ef)
  call random_number(dx_cf)
  call random_number(dy_cf)
  call random_number(dz_cf)
  call random_number(x_cf)
  call random_number(y_cf)
  call random_number(z_cf)
  call random_number(e_nu_c_bar)
  call random_number(f_nu_e_bar)
  call random_number(r_shock)
  call random_number(r_shock_mn)
  call random_number(r_shock_mx)
  call random_number(tau_adv)
  call random_number(tau_heat_nu)
  call random_number(tau_heat_nuc)
  call random_number(r_nse)
  call random_number(rho_c)
  call random_number(t_c)
  call random_number(ye_c)
  call random_number(u_c)
  call random_number(v_c)
  call random_number(w_c)
  call random_number(pMD)
  call random_number(sMD)
  call random_number(dudt_nuc)
  call random_number(dudt_nu)
  call random_number(grav_x_c)
  call random_number(grav_y_c)
  call random_number(grav_z_c)
  call random_number(rsphere_mean)
  call random_number(dsphere_mean)
  call random_number(tsphere_mean)
  call random_number(msphere_mean)
  call random_number(esphere_mean)
  call random_number(r_gain)
  call random_number(e_rad)
  call random_number(elec_rad)
  call random_number(unurad)
  call random_number(nnurad)
  call random_number(nse_c)
  call random_number(a_nuc_rep_c)
  call random_number(z_nuc_rep_c)
  call random_number(be_nuc_rep_c)
  call random_number(uburn_c)
  call random_number(duesrc)
  call random_number(unukrad)
  call random_number(unujrad)
  call random_number(nnukrad)
  call random_number(nnujrad)
  call random_number(nu_r)
  call random_number(nu_rt)
  call random_number(nu_rho)
  call random_number(nu_rhot)
  call random_number(psi0dat)
  call random_number(psi1dat)
  call random_number(xn_c)
  call random_number(dnurad)
  call random_number(psi0_c)
  
  do iCycle = 1, endCycle
    call model_write_hdf5()
    if(myid == 0 .and. mod(ncycle, 10) == 0)then
      print*, 'Performance:', ncycle*1.0/io_walltime, ' file/sec'
    end if
    ncycle = ncycle+1
    time = time + 0.5
  end do
  
  call MPI_FINALIZE(mpiError)

end program hdf5Test
