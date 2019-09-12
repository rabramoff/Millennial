! main program start
  PROGRAM Millennial
! History
! Xiaofeng Xu created this code program to play with Millennial model structure (ICOS workshop Mar 14-16, 2016 in Boulder, CO)
! The code is created in May - June 2016, solely by Xiaofeng XU (xxu@mail.sdsu.edu)
! This is a toy verion of the Millennial model (C only version, N P will be added in future updates)

  implicit none
  integer,parameter   :: r8 = selected_real_kind(12)  !8 byte real
  integer         :: nr, i, n, flag_output  , flag_annual     
  
 !  for initile file
  character(len = 256) :: initialfile
  character(len = 256) :: soilparafile

! for output file
  character(len = 256) :: outputfile
! end of defining output file

! the input data: driving forces
  real(r8), dimension(:), allocatable :: forc_st
  real(r8), dimension(:), allocatable :: forc_sw
! end of driving forces

! key variables to drive this model : semi-driving forces
  real(r8), dimension(:), allocatable :: forc_npp
  real(r8), dimension(:), allocatable :: forc_roots
  real(r8), dimension(:), allocatable :: forc_exoenzyme !if modified this could be calcuated based on biomass and limitation of N or P
! end of key variables to drive this model : semi-driving forces

!!  key variables to track the system over time
! pools 
  real(r8), dimension(:), allocatable :: LMWC
  real(r8), dimension(:), allocatable :: POM
  real(r8), dimension(:), allocatable :: MB
  real(r8), dimension(:), allocatable :: MINERAL
  real(r8), dimension(:), allocatable :: SOILAGG
! end of pools

! flux
  real(r8), dimension(:), allocatable :: f_LM_leaching    ! LM leaching
  real(r8), dimension(:), allocatable :: f_MI_LM_des    ! Mineral to LMWC
  real(r8), dimension(:), allocatable :: f_LM_MI_sor    ! LMWC to Mineral C
  real(r8), dimension(:), allocatable :: f_LM_MB_uptake ! LMWC to microbial C
  real(r8), dimension(:), allocatable :: f_PO_LM_dep    ! PO 
  real(r8), dimension(:), allocatable :: f_MB_MI_sor
  real(r8), dimension(:), allocatable :: f_PO_SO_agg
  real(r8), dimension(:), allocatable :: f_MI_SO_agg
  real(r8), dimension(:), allocatable :: f_SO_PO_break
  real(r8), dimension(:), allocatable :: f_SO_MI_break
  real(r8), dimension(:), allocatable :: f_MB_atm
! end of flux


! soil properties over time
  real(r8), dimension(:), allocatable :: psi_real
! end of flux

! soil properties
  real      :: sand
  real      :: clay
  real      :: silt
  real      :: maxpsi
  real      :: vwc
  real      :: vwcsat
  real      :: smp_l
  real      :: psisat
  real      :: organic
  real      :: psi
  
  real(r8)    :: k_leaching    
  real(r8)    :: klmc_min     
  real(r8)    :: Qmax   
  real(r8)    :: klmc   
  real(r8)    :: kes    
  real(r8)    :: CUEref
  real(r8)    :: CUET   
  real(r8)    :: Taeref 
  real(r8)    :: Vpom_lmc   
  real(r8)    :: kpom   
  real(r8)    :: k_POMes  
  real(r8)    :: kmic_min
  real(r8)    :: kmic
  real(r8)    :: Vpom_agg
  real(r8)    :: kpom_agg 
  real(r8)    :: Vmin_agg
  real(r8)    :: kmin_agg
  real(r8)    :: AGGmax
  real(r8)    :: kagg
! end

  real(r8)      :: initial_pom
  real(r8)      :: initial_lmwc
  real(r8)      :: initial_mb
  real(r8)      :: initial_mineral
  real(r8)      :: initial_soilagg
! end of key variables

  integer, parameter    :: soil_par_num = 25
! character(len=256)      :: soil_par_f = './soilpara_in'   ! local file name
  integer           :: ier                      ! error code
  character(len=40)     :: soil_par_name(soil_par_num)  ! parameter name
  real(r8)          :: dummy(soil_par_num)

  common  /global/ &
    k_leaching, &  
    klmc_min, &   
    Qmax, & 
    klmc, & 
    kes, &
    CUEref, &
    CUET, & 
    Taeref, &
    Vpom_lmc, & 
    kpom, & 
    k_POMes, &    
    kmic_min, & 
    kmic, & 
    Vpom_agg, &
    kpom_agg, &   
    Vmin_agg, &
    kmin_agg, &
    AGGmax, &
    kagg
  
  write(*,*) "This is the toy version of the millienum model at a daily time step"
  
  write(*,*) "Pleae enter the number for total simulation steps"
  read(*,*) nr

  write(*,*) "please enter the name of parameter file:"
  read(*,*) soilparafile
  
  write(*,*) "Do you want to save model output! 1 for YES, 0 for NO"
  read(*,*) flag_output
  write(*,*) "annaul output or daily? 1 for annual, 0 for daily"
  read(*,*) flag_annual
  if(flag_output == 1) then 
  write(*,*) "please enter the name of file for saving model output!"
  read(*,*) outputfile
  end if
  
! allocate space for key input data
  allocate(forc_st(1:nr))
  allocate(forc_sw(1:nr))
  allocate(forc_npp(1:nr))
  allocate(forc_roots(1:nr))
  allocate(forc_exoenzyme(1:nr))
  allocate(psi_real(1:nr))
  
  allocate(LMWC(1:nr))
  allocate(POM(1:nr))
  allocate(MB(1:nr))
  allocate(MINERAL(1:nr))
  allocate(SOILAGG(1:nr))
  
  allocate(f_LM_leaching(1:nr))
  allocate(f_MI_LM_des(1:nr))
  allocate(f_LM_MI_sor(1:nr))
  allocate(f_LM_MB_uptake(1:nr))
  allocate(f_PO_LM_dep(1:nr))
  allocate(f_MB_MI_sor(1:nr))
  allocate(f_PO_SO_agg(1:nr))
  allocate(f_MI_SO_agg(1:nr))
  allocate(f_SO_PO_break(1:nr))
  allocate(f_SO_MI_break(1:nr))
  allocate(f_MB_atm(1:nr))
! end of the allocation
  
  write(*,*) 'Attempting to read soil parameters .....'
  open(unit = 10, file=soilparafile)
  do i = 1, soil_par_num
  read (10,*,iostat=ier) soil_par_name(i), dummy(i)
  print *, dummy(i)
  if (ier /= 0) then
  write(*,*)'soilpara: error in reading in soilpara_in'
  end if
  end do
  close(10)

! Assign values
  i = 1
  clay        = dummy(i); i = i + 1
  sand        = dummy(i); i = i + 1
  silt        = dummy(i); i = i + 1
  maxpsi      = dummy(i); i = i + 1
  vwcsat      = dummy(i); i = i + 1
  organic     = dummy(i); i = i + 1
  k_leaching      = dummy(i); i = i + 1
  klmc_min      = dummy(i); i = i + 1
  Qmax      = dummy(i); i = i + 1
  klmc        = dummy(i); i = i + 1
  kes       = dummy(i); i = i + 1
  CUEref      = dummy(i); i = i + 1
  CUET      = dummy(i); i = i + 1
  Taeref      = dummy(i); i = i + 1
  Vpom_lmc      = dummy(i); i = i + 1
  kpom      = dummy(i); i = i + 1
  k_POMes     = dummy(i); i = i + 1
  kmic_min      = dummy(i); i = i + 1
  kmic        = dummy(i); i = i + 1
  Vpom_agg      = dummy(i); i = i + 1
  kpom_agg      = dummy(i); i = i + 1
  Vmin_agg      = dummy(i); i = i + 1
  kmin_agg      = dummy(i); i = i + 1
  AGGmax      = dummy(i); i = i +1
  kagg        = dummy(i)

  AGGmax = AGGmax * (0.0265 * clay * 100.0 + 0.1351)
  write(*,*) "Model inializing! "
  write(*,*) "please enter the name of file for initilizing the model"
  read(*,*) initialfile
    
  open(unit = 11, file=initialfile)
    
  read (11,*,iostat=ier) initial_pom, initial_lmwc, initial_mb, initial_mineral, initial_soilagg

  if (ier /= 0) then
  write(*,*)'model inializing failed !'
  else
  write(*,*) "model inialization finished !"
  end if
  close(11)

  print *, 'read data start'
  call readdata(nr, forc_st, forc_sw, forc_npp)
  print *, 'read data end'
  
  LMWC(1)=initial_lmwc
  POM(1)=initial_pom
  MB(1)=initial_mb
  MINERAL(1)=initial_mineral
  SOILAGG(1)=initial_soilagg
  
  do n = 1, nr
  vwc = forc_sw(n)
call soilpsi(sand, clay, silt, vwc, vwcsat, organic, psisat, psi, smp_l)
  psi_real(n) = psi

call decomp(forc_st(n), forc_sw(n), psi_real(n), forc_npp(n), forc_roots(n), &
    forc_exoenzyme(n), clay, LMWC(n), POM(n), MB(n), MINERAL(n), SOILAGG(n), f_LM_leaching(n), f_MI_LM_des(n),&
    f_LM_MI_sor(n), f_LM_MB_uptake(n),f_PO_LM_dep(n), f_MB_MI_sor(n), f_PO_SO_agg(n), f_MI_SO_agg(n),&
    f_SO_PO_break(n), f_SO_MI_break(n),f_MB_atm(n))
 
! upading the pool after each iteration 
  if(n < nr) then
  LMWC(n+1)=LMWC(n)
  POM(n+1)=POM(n)
  MB(n+1)=MB(n)
  MINERAL(n+1)=MINERAL(n)
  SOILAGG(n+1)=SOILAGG(n)
  endif
 
  end do
  
  if(flag_output ==1) then
call writeoutput(flag_annual, nr, forc_st, forc_sw, forc_npp, forc_roots, &
    forc_exoenzyme, LMWC, POM, MB, MINERAL, SOILAGG, f_LM_leaching, f_MI_LM_des,&
    f_LM_MI_sor, f_LM_MB_uptake, f_PO_LM_dep, f_MB_MI_sor, f_PO_SO_agg, f_MI_SO_agg,&
    f_SO_PO_break, f_SO_MI_break, f_MB_atm, outputfile)
  end if

  deallocate(forc_st)
  deallocate(forc_sw)
  deallocate(forc_npp)
  deallocate(forc_roots)
  deallocate(forc_exoenzyme)
  deallocate(psi_real)
  
  deallocate(LMWC)
  deallocate(POM)
  deallocate(MB)
  deallocate(MINERAL)
  deallocate(SOILAGG)
  
  deallocate(f_LM_leaching)
  deallocate(f_MI_LM_des)
  deallocate(f_LM_MI_sor)
  deallocate(f_LM_MB_uptake)
  deallocate(f_PO_LM_dep)
  deallocate(f_MB_MI_sor)
  deallocate(f_PO_SO_agg)
  deallocate(f_MI_SO_agg)
  deallocate(f_SO_PO_break)
  deallocate(f_SO_MI_break)
  deallocate(f_MB_atm)
! end of the allocation
  
  stop
END PROGRAM Millennial
! main program end


! read data subroutine start
subroutine readdata(nr, forc_st, forc_sw, forc_npp)
  implicit none
  integer,parameter     :: r8 = selected_real_kind(12) ! 8 byte real

  integer, intent(in)   :: nr
  real(r8), intent(out)   :: forc_st(1:nr)
  real(r8), intent(out)   :: forc_sw(1:nr)
  real(r8), intent(out)   :: forc_npp(1:nr)

  integer :: n, ier, i
  character(len=256)    :: filename
  print *, "please enter the filename for input data:"
  read (*,*) filename
  open (1001, file=filename, IOSTAT=ier)
  if(ier /= 0) then
  write (*,*) filename, 'file does not exist!'
  end if

  do n = 1, 365 !nr !xiaofeng xu made this change to recycle the data, avoiding read in large dataset
  read(1001,*,iostat=ier) forc_st(n), forc_sw(n), forc_npp(n)

  if (ier /= 0) then
  write(*,*) 'error in reading input data'
  end if
  end do
  close(1001)
  print *, "reading forcing data finished"
  
  do n = 366, nr
  i = mod(n-1, 365) + 1
  forc_st(n) = forc_st(i)
  forc_sw(n) = forc_sw(i)
  forc_npp(n) = forc_npp(i)
  end do
! read
end subroutine readdata
! read data subroutine end

  
subroutine writeoutput(flag_annual, nr, forc_st, forc_sw, forc_npp, forc_roots, &
    forc_exoenzyme, LMWC, POM, MB, MINERAL, SOILAGG, f_LM_leaching, f_MI_LM_des, &
    f_LM_MI_sor, f_LM_MB_uptake, f_PO_LM_dep, f_MB_MI_sor,f_PO_SO_agg, f_MI_SO_agg, &
    f_SO_PO_break, f_SO_MI_break, f_MB_atm, outputfile)
  
  implicit none
  integer,parameter     :: r8 = selected_real_kind(12) ! 8 byte real

  integer, intent(in)   :: flag_annual
  integer, intent(in)   :: nr
  real(r8), intent(in)    :: forc_st(1:nr)
  real(r8), intent(in)    :: forc_sw(1:nr)
  real(r8), intent(in)    :: forc_npp(1:nr)
  real(r8), intent(in)    :: forc_roots(1:nr)
  real(r8), intent(in)    :: forc_exoenzyme(1:nr)
  real(r8), intent(in)    :: LMWC(1:nr)
  real(r8), intent(in)    :: POM(1:nr)
  real(r8), intent(in)    :: MB(1:nr)
  real(r8), intent(in)    :: MINERAL(1:nr)
  real(r8), intent(in)    :: SOILAGG(1:nr)
  real(r8), intent(in)    :: f_LM_leaching(1:nr)
  real(r8), intent(in)    :: f_MI_LM_des(1:nr)
  real(r8), intent(in)    :: f_LM_MI_sor(1:nr)
  real(r8), intent(in)    :: f_LM_MB_uptake(1:nr)
  real(r8), intent(in)    :: f_PO_LM_dep(1:nr)
  real(r8), intent(in)    :: f_MB_MI_sor(1:nr)
  real(r8), intent(in)    :: f_PO_SO_agg(1:nr)
  real(r8), intent(in)    :: f_MI_SO_agg(1:nr)
  real(r8), intent(in)    :: f_SO_PO_break(1:nr)
  real(r8), intent(in)    :: f_SO_MI_break(1:nr)
  real(r8), intent(in)    :: f_MB_atm(1:nr)

  
  character(len = 256), intent(in)  :: outputfile
  
  integer         :: n, ier, j
  real(r8)        :: clmw, cpom, cmb, cmoa, cagg
  
  open (1000, FILE=outputfile)
! annual output
if(flag_annual == 1) then
do n = 1, nr / 365
  clmw = 0.0
  cpom = 0.0
  cmb = 0.0
  cmoa = 0.0
  cagg = 0.0
  
do j = 1, 365
  clmw = clmw + LMWC((n - 1) * 365 + j)
  cpom = cpom + POM((n - 1) * 365 + j)
  cmb = cmb + MB((n - 1) * 365 + j)
  cmoa = cmoa + MINERAL((n - 1) * 365 + j)
  cagg = cagg + SOILAGG((n - 1) * 365 + j)
end do
  write(1000,*,iostat=ier)  n, &
  clmw / 365., cpom/365., cmb/365., cmoa/365., cagg/365.
  
  if (ier /= 0) then
  write(*,*) 'error in writing output'
  end if
end do
! annual output

else
! daily output
  do n = 1, nr
  write(1000,*,iostat=ier)  n, &
! forc_st(n), forc_sw(n), forc_npp(n), forc_roots(n), forc_exoenzyme(n), &
  LMWC(n), POM(n), MB(n), MINERAL(n), SOILAGG(n)!, &
! f_LM_leaching(n), f_MI_LM_des(n),f_LM_MI_sor(n), f_LM_MB_uptake(n),f_PO_LM_dep(n), &
! f_MB_MI_sor(n),f_PO_SO_agg(n), f_MI_SO_agg(n),f_SO_PO_break(n), f_SO_MI_break(n), f_MB_atm(n)!
  if (ier /= 0) then
  write(*,*) 'error in writing output'
  end if
  end do
! daily output
end if
  close(1000)
end subroutine writeoutput
! write output subroutine end

! decomposition subroutine start
subroutine decomp(forc_st, forc_sw, psi, forc_npp, forc_roots, &
    forc_exoenzyme, clay, LMWC, POM, MB, MINERAL, SOILAGG, f_LM_leaching, f_MI_LM_des,&
    f_LM_MI_sor, f_LM_MB_uptake,f_PO_LM_dep, f_MB_MI_sor,f_PO_SO_agg, f_MI_SO_agg,&
    f_SO_PO_break, f_SO_MI_break, f_MB_atm)
    
  implicit none
  integer,parameter :: r8 = selected_real_kind(12)  ! 8 byte real
  integer,parameter :: r6 = selected_real_kind(8)   ! 8 byte real
  real(r8), intent(in) :: forc_st           ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
  real(r8), intent(in) :: forc_sw           ! soil moisture (fraction)
  real(r8), intent(in) :: psi               ! soil water potential at saturation for CN code (MPa)
  real(r8), intent(in) :: forc_npp
  real(r8), intent(in) :: forc_roots
! real(r8), intent(in) :: pH

  real(r8),intent(inout)  :: forc_exoenzyme 
  real, intent(inout) :: clay 
  real(r8),intent(inout)  :: LMWC   
  real(r8),intent(inout)  :: POM      
  real(r8),intent(inout)  :: MB      
  real(r8),intent(inout)  :: MINERAL     
  real(r8),intent(inout)  :: SOILAGG     
  real(r8),intent(inout)  :: f_LM_leaching     
  real(r8),intent(inout)  :: f_MI_LM_des     
  real(r8),intent(inout)  :: f_LM_MI_sor     
  real(r8),intent(inout)  :: f_LM_MB_uptake      
  real(r8),intent(inout)  :: f_PO_LM_dep    
  real(r8),intent(inout)  :: f_MB_MI_sor           
  real(r8),intent(inout)  :: f_PO_SO_agg           
  real(r8),intent(inout)  :: f_MI_SO_agg           
  real(r8),intent(inout)  :: f_SO_PO_break    
  real(r8),intent(inout)  :: f_SO_MI_break
  real(r8),intent(inout)  :: f_MB_atm
  
  real(r8) :: k_leaching 
  real(r8) :: klmc_min      
  real(r8) :: Qmax    
  real(r8) :: klmc    
  real(r8) :: kes   
  real(r8) :: CUEref
  real(r8) :: CUET    
  real(r8) :: Taeref  
  real(r8) :: Vpom_lmc  
  real(r8) :: kpom    
  real(r8) :: k_POMes 
  real(r8) :: kmic_min
  real(r8) :: kmic      
  real(r8) :: Vpom_agg
  real(r8) :: kpom_agg
  real(r8) :: Vmin_agg
  real(r8) :: kmin_agg
  real(r8) :: AGGmax
  real(r8) :: kagg

  ! local pointers to implicit out scalars
  !
  ! !OTHER LOCAL VARIABLES:

  real    :: temp, temp2, temp3 ! temporary variables
  real    :: psi_tem1, psi_tem2
  real    :: k_sorption             ! temporar variable for k of sorption
! real    :: Qmax       ! maximum sorption capacity  mg / kg (mayes et al, 2012, SSSAJ)
  real(r8)  :: t_scalar           ! soil temperature scalar for decomp
  real(r8)  :: t_scalar_mb        ! soil temperature scalar for decomp
  real(r8)  :: minpsi, maxpsi       ! limits for soil water scalar for decomp
! real    :: psi                        ! temporary soilpsi for water scalar
  real    :: w_scalar           ! soil water scalar for decomp
  real    :: rate_scalar        ! combined rate scalar for decomp
  real    :: pH
  real  :: f_SO_break
  real(r8)  :: t_scalar_reverse       ! soil temperature scalar for decomp
  real(r8)  :: w_scalar_reverse       ! soil temperature scalar for decomp
  
  !~ !-----------------------------------------------------------------------

  common  /global/ &
    k_leaching, & 
    klmc_min, &   
    Qmax, & 
    klmc, & 
    kes, &
    CUEref, &
    CUET, & 
    Taeref, &
    Vpom_lmc, & 
    kpom, & 
    k_POMes, &    
    kmic_min, & 
    kmic, &   
    Vpom_agg, &
    kpom_agg, &
    Vmin_agg, &
    kmin_agg, &
    AGGmax, &
    kagg

  t_scalar = 0._r8
  t_scalar_reverse = 0._r8
  temp = (forc_st - 15._r8) / 10._r8
  t_scalar = t_scalar + 2. **(temp)
  t_scalar_reverse = t_scalar_reverse + 0.5**(temp)
  
  t_scalar_mb = 0._r8
  temp = (forc_st - 15._r8) / 10._r8
  t_scalar_mb = t_scalar_mb + 2.**(temp)  
  
  minpsi = -10.0_r8
  w_scalar = 0._r8
  maxpsi = -0.01_r8
  pH = 7.0
  
! xiaofeng replaced above codes with following
  !if (psi > minpsi) then
  !w_scalar = w_scalar + (psi-minpsi)*(psi-maxpsi)/((psi-minpsi)*(psi-maxpsi) - &
  !  (psi-(maxpsi-(maxpsi-minpsi)/3.))*(psi-(maxpsi-(maxpsi-minpsi)/3.)))
  !end if
  !w_scalar = w_scalar ** 0.5

  !#century temperature function
  t_scalar = (11.75 + (29.7 / 3.1415926) * ATAN(real(3.1415926*0.031*(forc_st - 15.4)))) / &
  (11.75 + (29.7 / 3.1415926) * ATAN(real(3.1415926 * 0.031 *(30.0 - 15.4))))
  t_scalar_mb = t_scalar
  
  !#century water function
  w_scalar = 1.0 / (1.0 + 30. * EXP(real(-9.0 * forc_sw)))

  ! LMWC -> out of sysem LWMMWC leaching
  if (LMWC > 0._r8) then
        f_LM_leaching = LMWC * k_leaching * t_scalar !* w_scalar ! Xiaofeng removed water impact, after review at GBC June,2017
  end if

  ! LMWC -> MINERAL: This desorption function is from Mayes 2012, SSAJ
  klmc_min = (10.0 ** (-0.186 * pH - 0.216)) / 24.0

  Qmax = 10.0 ** (0.297 * log(clay * 100.0) + 2.355 + 0.50) ! 1.35 is bulk density to convert Q from mg/kg to g/m2

  temp = (klmc_min * Qmax * LMWC ) / (2. + klmc_min * LMWC) - MINERAL

  f_LM_MI_sor = (temp / Qmax + 0.0015) * LMWC / 50. * t_scalar * w_scalar !* t_scalar * w_scalar !* (LMWC / 200) * (LMWC / 200)

  if (f_LM_MI_sor < (LMWC * 0.9)) then
  f_LM_MI_sor = f_LM_MI_sor 
  else
  f_LM_MI_sor = LMWC * 0.9
  end if
  
  ! LMWC -> MB
  if (LMWC > 0._r8) then
  f_LM_MB_uptake = LMWC * klmc * t_scalar * w_scalar * MB / (MB + kes) * LMWC / (20. + LMWC)
  temp2 = f_LM_MB_uptake * (1. - (CUEref + CUET * (forc_st - Taeref)))
  if(temp2 < 0._r8) then
  temp2 = 0_r8
  end if
  f_LM_MB_uptake = f_LM_MB_uptake - temp2
  end if

  ! POM -> LMWC
  if (POM > 0._r8) then
        f_PO_LM_dep = Vpom_lmc * POM / (POM + kpom) * t_scalar * w_scalar !* (1. - MB / (MB + k_POMes)) 
  end if

  if(f_PO_LM_dep > (0.9 * POM)) then
  f_PO_LM_dep = 0.9 * POM
  end if
  
  if (MB > 0._r8 .and. MINERAL < Qmax) then
  f_MB_MI_sor = MB * kmic * 0.15 * t_scalar_mb * w_scalar  !* (MB / 200) * (MB / 200)
  else
  f_MB_MI_sor = 0.
  end if
  
  if(f_MB_MI_sor > 0.9 * MB) then
  f_MB_MI_sor = 0.9 * MB
  end if
  if(f_MB_MI_sor < 0.) then
  f_MB_MI_sor = 0.
  end if
  
  ! MB -> ATM
  if (MB > 0._r8) then
        f_MB_atm = temp2 + MB * kmic * t_scalar_mb * w_scalar 
  end if
  
  ! POM -> SOILAGG
  if (POM > 0._r8) then
        f_PO_SO_agg = Vpom_agg * POM / (kpom_agg + POM) * (1. - SOILAGG / AGGmax) * t_scalar * w_scalar
  end if
  
  if(f_PO_SO_agg > 0.9 * POM) then
  f_PO_SO_agg = 0.9 * POM
  end if

  ! MINERAL -> SOILAGG
  if (MINERAL > 0._r8) then
        f_MI_SO_agg = Vmin_agg * MINERAL / (kmin_agg + MINERAL) * (1. - SOILAGG / AGGmax) !* t_scalar * w_scalar
  end if

  if(f_MI_SO_agg>0.9 * MINERAL) then
  f_MI_SO_agg = 0.9 * MINERAL
  end if
  
  ! SOILAGG -> MINERAL
  if (SOILAGG > 0._r8) then
        f_SO_break = SOILAGG * kagg * t_scalar * w_scalar
  f_SO_PO_break = f_SO_break * 1.5 / 3.
  f_SO_MI_break = f_SO_break * 1.5 / 3.
  end if

  if((f_PO_LM_dep + f_PO_SO_agg) > POM) then
  temp3 = POM / (f_PO_LM_dep + f_PO_SO_agg)
  f_PO_LM_dep = f_PO_LM_dep * temp3
  f_PO_SO_agg = f_PO_SO_agg * temp3
  end if
  
  LMWC = LMWC + (f_PO_LM_dep + f_MI_LM_des - f_LM_leaching - f_LM_MI_sor - f_LM_MB_uptake - temp2) + forc_npp / 3.
  
  POM = POM + (f_SO_PO_break - f_PO_LM_dep - f_PO_SO_agg) + forc_npp * 2. / 3.
  
  MB = MB + (f_LM_MB_uptake - f_MB_MI_sor - f_MB_atm)
  
  MINERAL = MINERAL + (f_LM_MI_sor + f_MB_MI_sor + f_SO_MI_break - f_MI_LM_des - f_MI_SO_agg)
  
  SOILAGG = SOILAGG + (f_PO_SO_agg + f_MI_SO_agg - f_SO_PO_break - f_SO_MI_break)
  
end subroutine decomp
  ! decomposition subroutine end

! hydrological properties start
subroutine soilpsi(sand, clay, silt, vwc, vwcsat, organic, psisat, psi, smp_l)
implicit none
  integer,parameter   :: r8 = selected_real_kind(12)    ! 8 byte real
  real, intent(in)    :: sand     ! percentage
  real, intent(in)    :: clay     ! percentage
  real, intent(in)    :: silt     ! percentage
  real, intent(in)    :: vwc
  real, intent(in)    :: vwcsat
  real, intent(in)    :: organic    ! read-in the organic matter content g / m3
  real, intent(out)     :: psisat
  real, intent(out)     :: psi
  real, intent(out)   :: smp_l    ! soil matric potential (mm)
  real          :: bsw    ! Clapp and Hornberger "b"
  real          :: bsw2     ! Clapp and Hornberger "b" for CN module
  real          :: smp    ! msoil matrix potential  it seems this is exactly same as smp_l
  real          :: sucsat     ! minimum soil suction
  real          :: s_node     ! soil wetness
  real          :: smpmin   ! restriction for min of soil potential
  real          :: om_frac  ! organic matter fraction
  real          :: om_b     ! Clapp Hornberger parameter for organic soil (letts, 2000) Line 188 of iniTimeConst.F90
  real          :: organic_max  ! orgnaic matter hwere oil is assumed to act like peat
  real          :: om_sucsat  ! saturated suction for organic matter (Lets, 2000)
 
  om_sucsat = 10.3_r8           ! saturated suction for organic matter (Lets, 2000)
  smpmin = -1._r8             ! restriction for min of soil potential line 750 of iniTimeConst.F90
  organic_max = 130._r8
  om_b = 2.7_r8
  if(vwc > vwcsat) then
  write(*,*) 'vwcsat is less than vwc'
  end if
 
  om_frac = min(organic / organic_max, 1._r8)
  sucsat = 10. * ( 10.**(1.88-0.0131*sand) )
  bsw = (1.-om_frac)*(2.91 + 0.159*clay) + om_frac*om_b   
  sucsat = (1.-om_frac)*sucsat + om_sucsat*om_frac  
  s_node = min(1.0, max(vwc / vwcsat, 0.01))      
  smp = max(smpmin, (-sucsat * (s_node ** (-bsw))))
  smp_l = smp
  bsw2 = -(3.1 + 0.157 * clay - 0.003 * sand)
  psisat = -(exp((1.54 - 0.0095*sand + 0.0063*(100.0-sand-clay))*log(10.0))*9.8e-5_r8)
  psi = psisat * ((vwc/vwcsat)**bsw2)
end subroutine soilpsi
! hydrological properties end