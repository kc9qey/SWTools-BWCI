!-----------------------------------------------------------------------
!    File:         IO_MODULE
!    Module:       io_module
!    Type:         Module
!    Author:       Reuben D Budiardja, Dept of Physics, UTK,
!                  Knoxville, TN 37996
!
!    Date:         03/21/08
!
!    Purpose:
!      Contains the subroutines necessary for IO
!
!-----------------------------------------------------------------------


MODULE io_module

  USE kind_module, ONLY : double
  USE numerical_module, ONLY : zero, half, one, epsilon
  USE data_module

  USE HDF5
  USE MPI
  
  
  IMPLICIT none
  
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: filename = 'RadHyd3D_output_'
  
  
  LOGICAL, PRIVATE                    :: io_initialized = .FALSE.
  INTEGER, PRIVATE                    :: nx              ! x-array extent
  INTEGER, PRIVATE                    :: ny              ! y-array extent globaly
  INTEGER, PRIVATE                    :: nz              ! z-array extent globaly
  INTEGER, PRIVATE                    :: nez             ! Number of Neutrino Energy-Zone 
  INTEGER, PRIVATE                    :: nnu             ! neutrino flavor array extent
  INTEGER, PRIVATE                    :: nnc             ! Composition array extent
  
  !-----------------------------------------------------------------------
  !        Variables related to the domain decomposition 
  !
  !  my_j_ray_dim : Number of ray owned by local processor in j direction
  !  my_k_ray_dim : Number of ray owned by local processor in k direction
  !  my_j_ray     : Index of ray in j dimension local to a processor
  !  my_k_ray     : Index of ray in k dimension local to a processor
  !  j_ray        : Global index of j ray
  !  k_ray        : Global index of k ray
  !  j_ray_min    : Min of global index of j ray for this processor
  !  k_ray_min    : Min of global index of k ray for this processor
  !  j_ray_max    : Max of global index of j ray for this processor
  !  k_ray_max    : Max of global index of k ray for this processor
  !  nproc_y      : Number of processors assigned to the y dimension
  !  nproc_z      : Number of processors assigned to the z dimension
  !  
  !-----------------------------------------------------------------------
  
  INTEGER, PRIVATE                    :: my_j_ray_dim
  INTEGER, PRIVATE                    :: my_k_ray_dim
  INTEGER, PRIVATE                    :: my_j_ray
  INTEGER, PRIVATE                    :: my_k_ray
  INTEGER, PRIVATE                    :: j_ray
  INTEGER, PRIVATE                    :: k_ray
  INTEGER, PRIVATE                    :: j_ray_min
  INTEGER, PRIVATE                    :: k_ray_min
  INTEGER, PRIVATE                    :: j_ray_max
  INTEGER, PRIVATE                    :: k_ray_max
  
  INTEGER, PRIVATE                    :: io_count
  REAL(KIND=double), PUBLIC           :: io_walltime  

  
  PUBLIC      :: model_write_hdf5
  
  PRIVATE     :: write_ray_hyperslab
  PRIVATE     :: write_ray_hyperslab_dbl_2d
  PRIVATE     :: write_ray_hyperslab_dbl_3d
  PRIVATE     :: write_ray_hyperslab_dbl_4d
  PRIVATE     :: write_ray_hyperslab_dbl_5d
  PRIVATE     :: write_ray_hyperslab_int_3d
  
  PRIVATE     :: write_1d_slab
  PRIVATE     :: write_1d_slab_int
  PRIVATE     :: write_1d_slab_double
  
  PRIVATE     :: read_ray_hyperslab
  PRIVATE     :: read_ray_hyperslab_dbl_2d
  PRIVATE     :: read_ray_hyperslab_dbl_3d
  PRIVATE     :: read_ray_hyperslab_dbl_4d
  PRIVATE     :: read_ray_hyperslab_dbl_5d
  PRIVATE     :: read_ray_hyperslab_int_3d
  
  PRIVATE     :: read_1d_slab
  PRIVATE     :: read_1d_slab_int
  PRIVATE     :: read_1d_slab_double
  
  INTERFACE write_ray_hyperslab
    MODULE PROCEDURE write_ray_hyperslab_int_3d
    MODULE PROCEDURE write_ray_hyperslab_dbl_2d
    MODULE PROCEDURE write_ray_hyperslab_dbl_3d
    MODULE PROCEDURE write_ray_hyperslab_dbl_4d
    MODULE PROCEDURE write_ray_hyperslab_dbl_5d
  END INTERFACE write_ray_hyperslab
  
  INTERFACE write_1d_slab
    MODULE PROCEDURE write_1d_slab_int
    MODULE PROCEDURE write_1d_slab_double
  END INTERFACE write_1d_slab

  INTERFACE read_ray_hyperslab
    MODULE PROCEDURE read_ray_hyperslab_int_3d
    MODULE PROCEDURE read_ray_hyperslab_dbl_2d
    MODULE PROCEDURE read_ray_hyperslab_dbl_3d
    MODULE PROCEDURE read_ray_hyperslab_dbl_4d
    MODULE PROCEDURE read_ray_hyperslab_dbl_5d
  END INTERFACE read_ray_hyperslab
  
  INTERFACE read_1d_slab
    MODULE PROCEDURE read_1d_slab_int
    MODULE PROCEDURE read_1d_slab_double
  END INTERFACE read_1d_slab


  CONTAINS

  
  SUBROUTINE initialized_io()
    
    INTEGER :: error
    
    !------------------------------------------------------------------------
    !       Initialize Variables related to the decomposition
    !------------------------------------------------------------------------
    
    !-- FIXME: Check if this is the best way to do it
  
    if(.NOT.io_initialized)THEN
      nx = size(x_cf)
      ny = size(y_cf)
      nz = size(z_cf)
      nez = size(psi0_c, dim=2)
      nnu = size(psi0_c, dim=3)
      nnc = size(xn_c, dim=2)

      my_j_ray_dim = ny/nproc_y
      my_k_ray_dim = nz/nproc_z
      
      j_ray_min = MOD(myid, nproc_y) * my_j_ray_dim + 1
      k_ray_min = (myid/nproc_y) * my_k_ray_dim + 1
      j_ray_max = MOD(myid, nproc_y) * my_j_ray_dim + my_j_ray_dim
      k_ray_max = (myid/nproc_y) * my_k_ray_dim + my_k_ray_dim
      
      io_walltime = zero
      io_count    = 0

      io_initialized = .TRUE.
    END IF
  
  END SUBROUTINE initialized_io


  SUBROUTINE model_write_hdf5()
    !-----------------------------------------------------------------------
    !
    !    File:         model_write_hdf5
    !    Type:         Subprogram
    !    Author:       Reuben Budiardja, Dept of Physics, UTK
    !
    !    Date:         3/20/08
    !
    !    Purpose:
    !      To dump the model configuration in HDF5 file.
    !
    !    Subprograms called:
    !        write_ray_hyperslab
    !        write_1d_slab
    !
    !    Input arguments:
    !  nx          : x-array extent
    !  nnu         : neutrino flavor array extent
    !
    !    Output arguments:
    !        none
    !
    !    Include files:
    !  edit_modulee, eos_snc_x_module, nu_dist_module, nu_energy_grid_module,
    !  radial_ray_module
    !
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    !       File, group, dataset, and dataspace Identifier 
    !-----------------------------------------------------------------------
    CHARACTER(LEN=12)                :: suffix
    INTEGER(HID_T)                   :: file_id         ! HDF5 File identifier  
    INTEGER(HID_T)                   :: group_id        ! HDF5 Group identifier
    INTEGER(HID_T)                   :: dataset_id      ! HDF5 dataset identifier
    INTEGER(HID_T)                   :: dataspace_id    ! HDF5 dataspace identifier
    INTEGER(HID_T)                   :: plist_id        ! HDF5 property list   
    
    INTEGER(HSIZE_T)                 :: thresshold    = 524288
    INTEGER(HSIZE_T)                 :: alignment     = 262144
    INTEGER(SIZE_T)                  :: sieve_buffer  = 524288
    
    INTEGER                          :: dset_rank       ! Dataset rank
    INTEGER                          :: error           ! Error Flag
    INTEGER                          :: FILE_INFO_TEMPLATE

    INTEGER(HSIZE_T), dimension(1)   :: datasize1d
    INTEGER(HSIZE_T), dimension(2)   :: datasize2d
    INTEGER(HSIZE_T), dimension(2)   :: mydatasize2d
    INTEGER(HSIZE_T), dimension(2)   :: slab_offset2d
    INTEGER(HSIZE_T), dimension(3)   :: datasize3d
    INTEGER(HSIZE_T), dimension(3)   :: mydatasize3d
    INTEGER(HSIZE_T), dimension(3)   :: slab_offset3d
    INTEGER(HSIZE_T), dimension(4)   :: datasize4d
    INTEGER(HSIZE_T), dimension(4)   :: mydatasize4d
    INTEGER(HSIZE_T), dimension(4)   :: slab_offset4d
    INTEGER(HSIZE_T), dimension(5)   :: datasize5d
    INTEGER(HSIZE_T), dimension(5)   :: mydatasize5d
    INTEGER(HSIZE_T), dimension(5)   :: slab_offset5d
    
    REAL(KIND=double)                :: io_startime
    REAL(KIND=double)                :: io_endtime
    
    CALL initialized_io()    
    
    !------------------------------------------------------------------------
    !       Create and Initialize File using Default Properties
    !------------------------------------------------------------------------

    WRITE(suffix, fmt='(i9.9,a3)') ncycle,'.h5'

    io_startime = MPI_WTIME()
    
    CALL h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_sieve_buf_size_f(plist_id, sieve_buffer, error)
    CALL h5pset_alignment_f(plist_id, thresshold, alignment, error);
    CALL MPI_Info_create(FILE_INFO_TEMPLATE, error)
    CALL MPI_Info_set(FILE_INFO_TEMPLATE, "romio_cb_write", "ENABLE", error)
    CALL MPI_Info_set(FILE_INFO_TEMPLATE, "romio_cb_read", "ENABLE", error) 
    CALL MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "33554432", error)
    CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

    CALL h5fcreate_f(filename//suffix, &
           H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    CALL h5pclose_f(plist_id, error)
    
    !------------------------------------------------------------------------
    !      Create Mesh Group and Write Mesh Associated Data
    !      FIXME: This only needs to be done by processor 0
    !------------------------------------------------------------------------
      
    CALL h5gcreate_f(file_id, '/mesh', group_id, error)

    !-- Write Mesh Metadata
    datasize1d(1) = 3
    CALL write_1d_slab('array_dimensions', (/nx,ny,nz/), group_id, &
           datasize1d, 'array_dimensions')
    
    !-- Write Problem Cycle and Time
    datasize1d(1) = 0
    CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
    CALL h5dcreate_f(group_id, 'time', H5T_NATIVE_DOUBLE, dataspace_id, &
                     dataset_id, error)
    iF(myid == 0) &
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, time, datasize1d, &
                      error)
    CALL h5dclose_f(dataset_id, error)
    CALL h5sclose_f(dataspace_id, error)

    datasize1d(1) = 0
    CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
    CALL h5dcreate_f(group_id, 'cycle', H5T_NATIVE_INTEGER, dataspace_id, &
                     dataset_id, error)
    IF(myid == 0) &
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, ncycle, datasize1d, &
                      error)
    CALL h5dclose_f(dataset_id, error)
    CALL h5sclose_f(dataspace_id, error)

    !-- Write Radial Index bound for MGFLD shifted radial arrays
    datasize1d(1) = 2
    CALL write_1d_slab('radial_index_bound', (/imin,imax/), group_id, &
           datasize1d, 'radial_index_bound')

    !-- Write Zone Face Coordinates
    datasize1d(1) = nx+1
    CALL write_1d_slab('x_ef', x_ef, group_id, &
           datasize1d, 'X Grid Zone Face', 'cm')
    datasize1d(1) = ny+1
    CALL write_1d_slab('y_ef', y_ef, group_id, &
           datasize1d, 'Y Grid Zone Face', 'rad')
    datasize1d(1) = nz+1
    CALL write_1d_slab('z_ef', z_ef, group_id, &
           datasize1d, 'Z Grid Zone Face', 'rad')
    
    !-- Write Zone Midpoint Coordinates
    datasize1d(1) = nx
    CALL write_1d_slab('x_cf', x_cf, group_id, &
           datasize1d, 'X Grid Zone Midpoint', 'cm')
    datasize1d(1) = ny
    CALL write_1d_slab('y_cf', y_cf, group_id, &
           datasize1d, 'Y Grid Zone Midpoint', 'rad')
    datasize1d(1) = nz
    CALL write_1d_slab('z_cf', z_cf, group_id, &
           datasize1d, 'Z Grid Zone Midpoint', 'rad')
    
    !-- Write Zone Width 
    datasize1d(1) = nx
    CALL write_1d_slab('dx_cf', dx_cf, group_id, &
           datasize1d, 'X Zone Width', 'cm')
    datasize1d(1) = ny
    CALL write_1d_slab('dy_cf', dy_cf, group_id, &
           datasize1d, 'Y Zone Width', 'rad')
    datasize1d(1) = nz
    CALL write_1d_slab('dz_cf', dz_cf, group_id, &
           datasize1d, 'Z Zone Width', 'rad')
    
    CALL h5gclose_f(group_id, error)
    
    !------------------------------------------------------------------------
    !      Create Physical Variables group and Write Physical Data
    !      Each processors write its share to a hyperslab of the data 
    !
    !              \\\\\  Model Restart and Plot Files /////
    !
    !------------------------------------------------------------------------
    
    CALL h5gcreate_f(file_id, '/physical_variables', group_id, error)
    
    !-----------------------------------------------------------------------
    !  Independent thermodynamic variables
    !-----------------------------------------------------------------------
    datasize3d = (/nx,ny,nz/)
    mydatasize3d = (/nx, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('rho_c', rho_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Average Density','gm/c^3')
    CALL write_ray_hyperslab('t_c', t_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Average Temperature','K')
    CALL write_ray_hyperslab('ye_c', ye_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Average Electron Fraction','K')
    
    !-----------------------------------------------------------------------
    !  Independent mechanical variables
    !-----------------------------------------------------------------------
    CALL write_ray_hyperslab('u_c', u_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered Average Velocity X Direction','cm/s')
    CALL write_ray_hyperslab('v_c', v_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered Average Velocity Y Direction','cm/s')
    CALL write_ray_hyperslab('w_c', u_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered Average Velocity Z Direction','cm/s')
    
    !-----------------------------------------------------------------------
    !  Independent radiation variables and bookkeeping arrays
    !-----------------------------------------------------------------------
    datasize5d = (/nx,nez,nnu,ny,nz/)
    mydatasize5d = (/nx, nez, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset5d = (/0,0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('psi0_c', psi0_c, group_id, &
           datasize5d, mydatasize5d, slab_offset5d, &
           'Zero Moment of the Neutrino Occupation Distribution')
    CALL write_ray_hyperslab('dnurad', dnurad, group_id, &
           datasize5d, mydatasize5d, slab_offset5d, &
           'Total number of Neutrinos Per Unit Energy That Have ' &
           //'Crossed Outer Boundary Radial Zone', 'MeV^-1')
           
    datasize4d = (/nez,nnu,ny,nz/)
    mydatasize4d = (/nez, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('unukrad', unukrad, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Cumulative Energy Emmitted from the Core For Each '& 
           // 'Neutrinos Type of Each Energy Zone', 'Ergs')
    
    datasize4d = (/nx,nnu,ny,nz/)
    mydatasize4d = (/nx, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('unujrad', unujrad, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Cumulative Energy for Each Neutrino Type '& 
           // 'Transported Across Radial Boundary Zone', 'Ergs')
    
    datasize2d = (/ny,nz/)
    mydatasize2d = (/my_j_ray_dim, my_k_ray_dim/)
    slab_offset2d = (/j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('e_rad', e_rad, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, & 
           'Cumulative Material Energy Entering (-) or Leaving (+) the Grid', &
           'Ergs')
    CALL write_ray_hyperslab('elec_rad', elec_rad, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, & 
           'Net Number of electrons Advected In (-) or Out (+) Of the Grid')
    
    datasize3d = (/nnu,ny,nz/)
    mydatasize3d = (/nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('unurad', unurad, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Cumulative Energy Emitted From the Core for Each Neutrino Type', &
           'Ergs')
    
    datasize4d = (/nez,nnu,ny,nz/)
    mydatasize4d = (/nez, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('nnukrad', nnukrad, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Cumulative Number Emmitted from the Core For Each '& 
           // 'Neutrinos Type of Each Energy Zone')
    
    datasize4d = (/nx,nnu,ny,nz/)
    mydatasize4d = (/nx, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('nnujrad', nnujrad, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Net Number For Each Neutrino Type' &
           // 'Transported Across Radial Boundary Zone')
    
    datasize3d = (/nnu,ny,nz/)
    mydatasize3d = (/nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('nnurad', nnurad, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Cumulative Number Of Each Neutrino Type Emmited by the Core')
           
    datasize4d = (/nez,nnu,ny,nz/)
    mydatasize4d = (/nez, nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('nu_r', nu_r, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Number of Neutrinos for Each Energy Group Radiated Across r_nurad')
    CALL write_ray_hyperslab('nu_rt', nu_rt, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Cumulative Number of Neutrinos for Each Energy Group Radiated ' &
           //'Across r_nurad')
    CALL write_ray_hyperslab('nu_rho', nu_rho, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Number of Neutrinos for Each Energy group Across rho_nurad '&
           // 'in Time dtnuradplot')
    CALL write_ray_hyperslab('nu_rhot', nu_rhot, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Cumulative Number of Neutrinos for Each Energy Group ' &
           //'Radiated Across rho_nurad')
    
    CALL write_ray_hyperslab('psi0dat', psi0dat, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Time Integrated psi0')
    CALL write_ray_hyperslab('psi1dat', psi1dat, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Time Integrated psi1')
    
    datasize4d = (/nx,nnc,ny,nz/)
    mydatasize4d = (/nx, nnc, my_j_ray_dim, my_k_ray_dim/)
    slab_offset4d = (/0,0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('xn_c', xn_c, group_id, &
           datasize4d, mydatasize4d, slab_offset4d, & 
           'Mass Fraction for Each Nucleus')
    
    datasize3d = (/nx,ny,nz/)
    mydatasize3d = (/nx, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('nse_c', nse_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'NSE Flag')
    CALL write_ray_hyperslab('a_nuc_rep_c', a_nuc_rep_c, group_id, &
            datasize3d, mydatasize3d, slab_offset3d, &
           'Mass Number of the Representative Heavy Nucleus')
    CALL write_ray_hyperslab('z_nuc_rep_c', z_nuc_rep_c, &
           group_id, datasize3d, mydatasize3d, slab_offset3d, &
           'Charge Number of the Representative Heavy Nucleus')
    CALL write_ray_hyperslab('be_nuc_rep_c', be_nuc_rep_c, &
           group_id, datasize3d, mydatasize3d, slab_offset3d, &
           'Binding Energy of the Representative Heavy Nucleus', 'MeV')
    CALL write_ray_hyperslab('uburn_c', uburn_c, &
           group_id, datasize3d, mydatasize3d, slab_offset3d, &
           'Cumulative Energy Generated in Zone by Nuclear Reaction', 'Ergs/gm')
    CALL write_ray_hyperslab('duesrc', duesrc, &
           group_id, datasize3d, mydatasize3d, slab_offset3d, &
           'Cumulative Energy Glitches per Unit Mass')
    
    datasize1d(1) = nx
    CALL write_1d_slab('e_nu_c_bar', e_nu_c_bar, group_id, &
           datasize1d, 'Angular Averaged Neutrino Energy Density', &
          'Ergs/cm^3')
    CALL write_1d_slab('f_nu_e_bar', f_nu_e_bar, group_id, &
           datasize1d, 'Angular Averaged Neutrino Energy Flux', &
           'Ergs/cm^2/s')
    
    datasize3d = (/nx,ny,nz/)
    mydatasize3d = (/nx, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('pMD', pMD, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Pressure')
    CALL write_ray_hyperslab('sMD', sMD, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Entropy')
    CALL write_ray_hyperslab('dudt_nuc', dudt_nuc, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Energy Generation Rate By Nuclear Reaction For ' &
           //'The Current Time Step', 'Ergs/gm')
    CALL write_ray_hyperslab('dudt_nu', dudt_nu, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Energy Deposition Rate By All Neutrinos in Radial Zone', &
           'Ergs/gm/s')
    CALL write_ray_hyperslab('grav_x_c', grav_x_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered X-Component of Gravitational Acceleration', &
           'cm/s^2/g')
    CALL write_ray_hyperslab('grav_y_c', grav_y_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered Y-Component of Gravitational Acceleration', &
           'cm/s^2/g')
    CALL write_ray_hyperslab('grav_z_c', grav_z_c, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Zone Centered Z-Component of Gravitational Acceleration', &
           'cm/s^2/g')
    
    datasize2d = (/ny,nz/)
    mydatasize2d = (/my_j_ray_dim, my_k_ray_dim/)
    slab_offset2d = (/j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('r_shock', r_shock, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Radius of Shock Maximum')
    CALL write_ray_hyperslab('r_shock_mn', r_shock_mn, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Minimum Estimated Shock Radius')
    CALL write_ray_hyperslab('r_shock_mx', r_shock_mx, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Maximum Estimated Shock Radius')
    CALL write_ray_hyperslab('tau_adv', tau_adv, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Advection Time Scale', 's')
    CALL write_ray_hyperslab('tau_heat_nu', tau_heat_nu, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Neutrino Heating Time Scale', 's')
    CALL write_ray_hyperslab('tau_heat_nuc', tau_heat_nuc, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Nuclear Heating Time Scale', 's')
    CALL write_ray_hyperslab('r_nse', r_nse, group_id, &
           datasize2d, mydatasize2d, slab_offset2d, &
           'Radius of NSE-nonNSE Boundary')
    
    datasize3d = (/nnu,ny,nz/)
    mydatasize3d = (/nnu, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('rsphere_mean', rsphere_mean, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Mean Neutrinosphere Radius')
    CALL write_ray_hyperslab('dsphere_mean', dsphere_mean, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Mean Neutrinosphere Density')
    CALL write_ray_hyperslab('tsphere_mean', tsphere_mean, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Mean Neutrinosphere Temperature')
    CALL write_ray_hyperslab('msphere_mean', msphere_mean, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Mean Neutrinosphere Enclosed Mass')
    CALL write_ray_hyperslab('esphere_mean', esphere_mean, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Mean Neutrinosphere Energy')
    
    datasize3d = (/nnu+1,ny,nz/)
    mydatasize3d = (/nnu+1, my_j_ray_dim, my_k_ray_dim/)
    slab_offset3d = (/0,j_ray_min-1,k_ray_min-1/)
    CALL write_ray_hyperslab('r_gain', r_gain, group_id, &
           datasize3d, mydatasize3d, slab_offset3d, &
           'Gain Radius')
    
    CALL h5gclose_f(group_id, error)

    !------------------------------------------------------------------------
    !      Cleanup
    !------------------------------------------------------------------------

    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    call MPI_Info_free(FILE_INFO_TEMPLATE, error)
    
    io_endtime = MPI_WTIME()
    io_walltime = io_walltime + (io_endtime-io_startime)
    io_count = io_count+1
    
    IF(myid == 0)THEN
      WRITE(*, '(a30,i9,a10,i5,a20,f10.5)')'*** HDF5 Model dump at cycle ', ncycle, &
        'IO Count:', io_count, 'IO elapsed time: ', io_walltime 
    ENDIF
    
    RETURN

  END SUBROUTINE model_write_hdf5
  
  
  

  SUBROUTINE write_1d_slab_int(name, value, group_id, datasize, &
               desc_option, unit_option)
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(1), INTENT(IN)  :: datasize
    INTEGER, DIMENSION(:), INTENT(IN)           :: value
    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(1, datasize, dataspace_id, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, &
                     dataspace_id, dataset_id, error)
    IF(myid == 0) &
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, &
                      value, datasize, error)
    CALL h5sclose_f(dataspace_id, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)

  END SUBROUTINE write_1d_slab_int

  SUBROUTINE write_1d_slab_double(name, value, group_id, datasize, &
               desc_option, unit_option)
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(1), INTENT(IN)  :: datasize
    REAL(kind=double), DIMENSION(:), INTENT(IN) :: value
    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(1, datasize, dataspace_id, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, &
                     dataspace_id, dataset_id, error)
    IF(myid == 0) &
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, &
                      value, datasize, error)
    CALL h5sclose_f(dataspace_id, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)

  END SUBROUTINE write_1d_slab_double


  SUBROUTINE write_ray_hyperslab_dbl_2d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(2), INTENT(IN)  :: global_datasize
    INTEGER(HSIZE_T), dimension(2), INTENT(IN)  :: local_datasize
    INTEGER(HSIZE_T), dimension(2), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:), INTENT(IN) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: plist_id
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(2, global_datasize, filespace, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, filespace, &
                     dataset_id, error)
    CALL h5sclose_f(filespace, error)
    
    CALL h5screate_simple_f(2, local_datasize, memspace, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           local_datasize, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)
    CALL h5pclose_f(plist_id, error)                     
    
  END SUBROUTINE write_ray_hyperslab_dbl_2d


  SUBROUTINE write_ray_hyperslab_dbl_3d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: global_datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: local_datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:), INTENT(IN) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: plist_id
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(3, global_datasize, filespace, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, filespace, &
                     dataset_id, error)
    CALL h5sclose_f(filespace, error)
    
    CALL h5screate_simple_f(3, local_datasize, memspace, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           local_datasize, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)
    CALL h5pclose_f(plist_id, error)                     
    
  END SUBROUTINE write_ray_hyperslab_dbl_3d


  SUBROUTINE write_ray_hyperslab_dbl_4d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(4), INTENT(IN)  :: global_datasize
    INTEGER(HSIZE_T), dimension(4), INTENT(IN)  :: local_datasize
    INTEGER(HSIZE_T), dimension(4), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:,:), INTENT(IN) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: plist_id
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(4, global_datasize, filespace, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, filespace, &
                     dataset_id, error)
    CALL h5sclose_f(filespace, error)
    
    CALL h5screate_simple_f(4, local_datasize, memspace, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           local_datasize, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)
    CALL h5pclose_f(plist_id, error)                     
    
  END SUBROUTINE write_ray_hyperslab_dbl_4d


  SUBROUTINE write_ray_hyperslab_dbl_5d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(5), INTENT(IN)  :: global_datasize
    INTEGER(HSIZE_T), dimension(5), INTENT(IN)  :: local_datasize
    INTEGER(HSIZE_T), dimension(5), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:,:,:), INTENT(IN) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: plist_id
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(5, global_datasize, filespace, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_DOUBLE, filespace, &
                     dataset_id, error)
    CALL h5sclose_f(filespace, error)
    
    CALL h5screate_simple_f(5, local_datasize, memspace, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           local_datasize, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)
    CALL h5pclose_f(plist_id, error)                     
    
  END SUBROUTINE write_ray_hyperslab_dbl_5d


  SUBROUTINE write_ray_hyperslab_int_3d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    CHARACTER(*), INTENT(IN)                    :: name
    CHARACTER(*), INTENT(IN), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(IN), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: global_datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: local_datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: slab_offset
    INTEGER, DIMENSION(:,:,:), INTENT(IN)       :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: plist_id
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
    INTEGER                                     :: error
    
    CALL h5screate_simple_f(3, global_datasize, filespace, error)
    CALL h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, filespace, &
                     dataset_id, error)
    CALL h5sclose_f(filespace, error)
    
    CALL h5screate_simple_f(3, local_datasize, memspace, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           local_datasize, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    IF(present(desc_option))THEN
      attr_len = len(desc_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Desc', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    IF(present(unit_option))THEN
      attr_len = len(unit_option)
      CALL h5screate_simple_f(1, adims, dataspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attr_len, error)
      CALL h5acreate_f(dataset_id, 'Unit', atype_id, dataspace_id, &
             attr_id, error)
      IF(myid == 0) &
        CALL h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END IF
    CALL h5dclose_f(dataset_id, error)
    CALL h5pclose_f(plist_id, error)                     
    
  END SUBROUTINE write_ray_hyperslab_int_3d
  
  
  SUBROUTINE read_1d_slab_int(name, value, group_id, datasize)
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(1), INTENT(IN)  :: datasize
    INTEGER, DIMENSION(:), INTENT(OUT)          :: value
    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER                                     :: error
    
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, &
                    value, datasize, error)
    CALL h5dclose_f(dataset_id, error)

  END SUBROUTINE read_1d_slab_int

  
  SUBROUTINE read_1d_slab_double(name, value, group_id, datasize)
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(1), INTENT(IN)  :: datasize
    REAL(kind=double), DIMENSION(:), INTENT(OUT):: value
    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER                                     :: error
    
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, &
                    value, datasize, error)
    CALL h5dclose_f(dataset_id, error)

  END SUBROUTINE read_1d_slab_double


  SUBROUTINE read_ray_hyperslab_dbl_2d(name, value, group_id, &
               datasize, slab_offset)
    
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(2), INTENT(IN)  :: datasize
    INTEGER(HSIZE_T), dimension(2), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:), INTENT(OUT) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HSIZE_T), dimension(2)              :: null_offset
    INTEGER                                     :: error
    
    null_offset = 0
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           datasize, error)
    CALL h5screate_simple_f(2, datasize, memspace, error)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, &
           datasize, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dataset_id, error)
    
  END SUBROUTINE read_ray_hyperslab_dbl_2d


  SUBROUTINE read_ray_hyperslab_dbl_3d(name, value, group_id, &
               datasize, slab_offset)
    
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:), INTENT(OUT) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HSIZE_T), dimension(3)              :: null_offset
    INTEGER                                     :: error
    
    null_offset = 0
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           datasize, error)
    CALL h5screate_simple_f(3, datasize, memspace, error)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, &
           datasize, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dataset_id, error)
    
  END SUBROUTINE read_ray_hyperslab_dbl_3d


  SUBROUTINE read_ray_hyperslab_dbl_4d(name, value, group_id, &
               datasize, slab_offset)
    
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(4), INTENT(IN)  :: datasize
    INTEGER(HSIZE_T), dimension(4), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:,:), INTENT(OUT) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HSIZE_T), dimension(4)              :: null_offset
    INTEGER                                     :: error
    
    null_offset = 0
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           datasize, error)
    CALL h5screate_simple_f(4, datasize, memspace, error)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, &
           datasize, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dataset_id, error)
    
  END SUBROUTINE read_ray_hyperslab_dbl_4d


  SUBROUTINE read_ray_hyperslab_dbl_5d(name, value, group_id, &
               datasize, slab_offset)
    
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(5), INTENT(IN)  :: datasize
    INTEGER(HSIZE_T), dimension(5), INTENT(IN)  :: slab_offset
    REAL(kind=double), DIMENSION(:,:,:,:,:), INTENT(OUT) :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HSIZE_T), dimension(5)              :: null_offset
    INTEGER                                     :: error
    
    null_offset = 0
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           datasize, error)
    CALL h5screate_simple_f(4, datasize, memspace, error)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, &
           datasize, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dataset_id, error)
    
  END SUBROUTINE read_ray_hyperslab_dbl_5d


  SUBROUTINE read_ray_hyperslab_int_3d(name, value, group_id, &
               datasize, slab_offset)
    
    CHARACTER(*), INTENT(IN)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: datasize
    INTEGER(HSIZE_T), dimension(3), INTENT(IN)  :: slab_offset
    INTEGER, DIMENSION(:,:,:), INTENT(OUT)       :: value
    
    INTEGER(HID_T)                              :: filespace    
    INTEGER(HID_T)                              :: memspace    
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HSIZE_T), dimension(3)              :: null_offset
    INTEGER                                     :: error
    
    null_offset = 0
    CALL h5dopen_f(group_id, name, dataset_id, error)
    CALL h5dget_space_f(dataset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, &
           datasize, error)
    CALL h5screate_simple_f(3, datasize, memspace, error)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, &
           datasize, error)
    CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dataset_id, error)
    
    
  END SUBROUTINE read_ray_hyperslab_int_3d
  
  
END MODULE io_module
