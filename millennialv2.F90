module soilType

  implicit none
  integer,parameter   :: p8 = selected_real_kind(12)  !8 byte real

! set parameters
type, public :: soil_type
  real(p8)    :: param_pi
  real(p8)    :: param_pa 
  real(p8)    :: kaff_pl    
  real(p8)    :: alpha_pl   
  real(p8)    :: eact_pl   
  real(p8)    :: gas_const
  real(p8)    :: rate_pa
  real(p8)    :: rate_break
  real(p8)    :: rate_leach
  real(p8)    :: param_pad
  real(p8)    :: param_p1
  real(p8)    :: param_p2
  real(p8)    :: param_pH
  real(p8)    :: param_bulkd
  real(p8)    :: param_c1
  real(p8)    :: param_c2
  real(p8)    :: param_clay
  real(p8)    :: kaff_lb
  real(p8)    :: alpha_lb
  real(p8)    :: eact_lb
  real(p8)    :: rate_bd
  real(p8)    :: rate_ma

  real(p8)    :: cue_ref
  real(p8)    :: cue_t
  real(p8)    :: tae_ref

  real(p8)    :: kaff_lm
  real(p8)    :: kaff_ld
  real(p8)    :: kaff_bm
! soil properties
  real(p8)      :: matpot ! soil matric potential kPa
  real(p8)      :: lambda ! dependence of respiration on soil matric potential kPa
  real(p8)      :: porosity ! total porosity
  real(p8)      :: kamin ! minimum relative rate in saturated soil
end type soil_type

contains
! read data subroutine start
subroutine readdata(nr, forc_st, forc_sw, forc_npp)
  implicit none
  integer, parameter     :: r8 = selected_real_kind(12) ! 8 byte real

  integer,  intent(in)   :: nr
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

  
subroutine writeoutput(flag_annual, nr, &
    scalar_wd, scalar_wb, kaff_eq, param_qmax, param_pb, vmax_pl, vmax_lb, &
    LMWC, POM, MIC, MAOM, AGG, f_AG_break, f_PO_AG,&
    f_MA_AG, f_PO_LM, f_LM_leach, f_LM_MA, f_LM_MB, f_MB_turn,&
    f_MA_LM, f_MB_atm, outputfile)
  
  implicit none
  integer, parameter      :: r8 = selected_real_kind(12) ! 8 byte real
  integer,  intent(in)    :: flag_annual
  integer,  intent(in)    :: nr
  real(r8), intent(in)    :: scalar_wd(1:nr)
  real(r8), intent(in)    :: scalar_wb(1:nr)
  real(r8), intent(in)    :: kaff_eq(1:nr)
  real(r8), intent(in)    :: param_qmax(1:nr)
  real(r8), intent(in)    :: param_pb(1:nr)
  real(r8), intent(in)    :: vmax_pl(1:nr)
  real(r8), intent(in)    :: vmax_lb(1:nr)
  real(r8), intent(in)    :: LMWC(1:nr)
  real(r8), intent(in)    :: POM(1:nr)
  real(r8), intent(in)    :: MIC(1:nr)
  real(r8), intent(in)    :: MAOM(1:nr)
  real(r8), intent(in)    :: AGG(1:nr)
  real(r8), intent(in)    :: f_AG_break(1:nr)
  real(r8), intent(in)    :: f_PO_AG(1:nr)
  real(r8), intent(in)    :: f_MA_AG(1:nr)
  real(r8), intent(in)    :: f_PO_LM(1:nr)
  real(r8), intent(in)    :: f_LM_leach(1:nr)
  real(r8), intent(in)    :: f_LM_MA(1:nr)
  real(r8), intent(in)    :: f_LM_MB(1:nr)
  real(r8), intent(in)    :: f_MB_turn(1:nr)
  real(r8), intent(in)    :: f_MA_LM(1:nr)
  real(r8), intent(in)    :: f_MB_atm(1:nr)

  
  character(len = 256), intent(in)  :: outputfile
  
  integer         :: n, ier, j
  real(r8)        :: clmw, cpom, cmic, cmoa, cagg, cswd, cswb, ckeq, cqmx, cppb
  real(r8)        :: cvpl, cvlb, cfab, cfpa, cfma, cfpl, cfll, cflm, cflb
  real(r8)        :: cfmt, cfml, cfmo 
  
  open (1000, FILE=outputfile)

! annual output
if(flag_annual == 1) then
  do n = 1, nr / 365
    clmw = 0._r8; cpom = 0._r8; cmic = 0._r8; cmoa = 0._r8; cagg = 0._r8
    cswd = 0._r8; cswb = 0._r8; ckeq = 0._r8; cqmx = 0._r8; cppb = 0._r8; 
    cvpl = 0._r8; cvlb = 0._r8; cfab = 0._r8; cfpa = 0._r8; cfma = 0._r8; 
    cfpl = 0._r8; cfll = 0._r8; cflm = 0._r8; cflb = 0._r8; cfmt = 0._r8; 
    cfml = 0._r8; cfmo = 0._r8
    
    do j = 1, 365
      clmw = clmw + LMWC((n - 1) * 365 + j)
      cpom = cpom + POM((n - 1) * 365 + j)
      cmic = cmic + MIC((n - 1) * 365 + j)
      cmoa = cmoa + MAOM((n - 1) * 365 + j)
      cagg = cagg + AGG((n - 1) * 365 + j)
      cswd = cswd + scalar_wd((n - 1) * 365 + j)
      cswb = cswb + scalar_wb((n - 1) * 365 + j)
      ckeq = ckeq + kaff_eq((n - 1) * 365 + j)
      cqmx = cqmx + param_qmax((n - 1) * 365 + j)
      cppb = cppb + param_pb((n - 1) * 365 + j)
      cvpl = cvpl + vmax_pl((n - 1) * 365 + j)
      cvlb = cvlb + vmax_lb((n - 1) * 365 + j)
      cfab = cfab + f_AG_break((n - 1) * 365 + j)
      cfpa = cfpa + f_PO_AG((n - 1) * 365 + j)
      cfma = cfma + f_MA_AG((n - 1) * 365 + j)
      cfpl = cfpl + f_PO_LM((n - 1) * 365 + j)
      cfll = cfll + f_LM_leach((n - 1) * 365 + j)
      cflm = cflm + f_LM_MA((n - 1) * 365 + j)
      cflb = cflb + f_LM_MB((n - 1) * 365 + j)
      cfmt = cfmt + f_MB_turn((n - 1) * 365 + j)
      cfml = cfml + f_MA_LM((n - 1) * 365 + j)
      cfmo = cfmo + f_MB_atm((n - 1) * 365 + j)

    end do
  write(1000,*,iostat=ier)  n, &
  clmw/365., cpom/365., cmic/365., cmoa/365., cagg/365., cswd/365., &
  cswb/365., ckeq/365., cqmx/365., cppb/365., cvpl/365., cvlb/365., cfab/365., &
  cfpa/365., cfma/365., cfpl/365., cfll/365., cflm/365., cflb/365., &
  cfmt/365., cfml/365., cfmo/365.
    
    if (ier /= 0) then
    write(*,*) 'error in writing output'
    end if
  end do
  ! annual output

  else

  ! daily output
    do n = 1, nr
    write(1000,*,iostat=ier)  n, &
    LMWC(n), POM(n), MIC(n), MAOM(n), AGG(n), &
    scalar_wd(n), scalar_wb(n), kaff_eq(n), param_qmax(n), param_pb(n), vmax_pl(n), vmax_lb(n), &
    f_AG_break(n), f_PO_AG(n),f_MA_AG(n), f_PO_LM(n), f_LM_leach(n), f_LM_MA(n), &
    f_LM_MB(n), f_MB_turn(n), f_MA_LM(n), f_MB_atm(n) 
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
subroutine decomp(this, forc_st, forc_npp, &
    scalar_wd, scalar_wb, kaff_eq, param_qmax, param_pb, vmax_pl, vmax_lb, &
    LMWC, POM, MIC, MAOM, AGG, f_AG_break, f_PO_AG,&
    f_MA_AG, f_PO_LM, f_LM_leach, f_LM_MA, f_LM_MB, f_MB_turn,&
    f_MA_LM, f_MB_atm)

  implicit none
  class(soil_type), intent(inout) :: this
  integer,parameter :: r8 = selected_real_kind(12)  ! 8 byte real
  real(r8), intent(in) :: forc_st           ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
  real(r8), intent(in) :: forc_npp
  real(r8), intent(in) :: scalar_wd
  real(r8), intent(in) :: scalar_wb
  real(r8), intent(in) :: kaff_eq 
  real(r8), intent(in) :: param_qmax
  real(r8), intent(inout) :: param_pb
  real(r8), intent(inout) :: vmax_pl
  real(r8), intent(inout) :: vmax_lb

  real(r8),intent(inout)  :: LMWC   
  real(r8),intent(inout)  :: POM      
  real(r8),intent(inout)  :: MIC      
  real(r8),intent(inout)  :: MAOM     
  real(r8),intent(inout)  :: AGG     
  real(r8),intent(inout)  :: f_AG_break    
  real(r8),intent(inout)  :: f_PO_AG    
  real(r8),intent(inout)  :: f_MA_AG     
  real(r8),intent(inout)  :: f_PO_LM    
  real(r8),intent(inout)  :: f_LM_leach 
  real(r8),intent(inout)  :: f_LM_MA         
  real(r8),intent(inout)  :: f_LM_MB          
  real(r8),intent(inout)  :: f_MB_turn         
  real(r8),intent(inout)  :: f_MA_LM    
  real(r8),intent(inout)  :: f_MB_atm

!Equation 3
  vmax_pl = this%alpha_pl * exp(-this%eact_pl / (this%gas_const * (forc_st + 273.15_r8)))

!Equation 2
  ! POM -> LMWC
  if (POM > 0._r8) then
  f_PO_LM = POM * vmax_pl * MIC * scalar_wd / (this%kaff_pl + POM + MIC) 
  endif
!Equation 6
  ! POM -> AGG
  if (POM > 0._r8) then
  f_PO_AG = this%rate_pa * scalar_wd * POM
  endif
!Equation 7
  ! AGG -> MAOM
  if (AGG > 0._r8) then
  f_AG_break = this%rate_break * scalar_wd * AGG
  endif
!Equation 9
  ! LMWC -> out of system leaching
  if (LMWC > 0._r8) then
  f_LM_leach = this%rate_leach * scalar_wd * LMWC 
  endif
!Equation 11
  this%kaff_lm = kaff_eq * this%param_pad !WORKING

!Equation 15
  this%kaff_ld = kaff_eq * (1 - this%param_pad)

!Equation 23
  this%kaff_bm = kaff_eq * this%param_pad

!Equation 10
  ! LMWC -> MAOM
  if (LMWC > 0._r8) then
  f_LM_MA = scalar_wd * this%kaff_lm * LMWC * (param_qmax - MAOM)
  endif
!Equation 14
  ! MAOM -> LMWC
  if (MAOM > 0._r8) then
  f_MA_LM = this%kaff_ld * MAOM
  endif
!Equation 17
  vmax_lb = this%alpha_lb * exp(-this%eact_lb / (this%gas_const * (forc_st + 273.15_r8)))

!Equation 16
  ! LMWC -> MIC
  if (LMWC > 0._r8) then
  f_LM_MB = vmax_lb * scalar_wb * LMWC * MIC / (this%kaff_lb + MIC + LMWC)
  endif
!Equation 18
  ! MIC -> MAOM/LMWC
  if (MIC > 0._r8) then
  f_MB_turn = this%rate_bd * MIC**2.0_r8
  endif
!Equation 20
  ! MAOM -> AGG
  if (MAOM > 0._r8) then
  f_MA_AG = this%rate_ma * scalar_wd * MAOM
  endif
!Equation 22
    if (MIC > 0._r8) then
  param_pb = scalar_wd * this%kaff_bm * MIC * (param_qmax - MAOM) 
    endif
!Equation 25
  ! microbial growth flux, but is not used in mass balance
!WORKING

!Equation 26
  ! MIC -> atmosphere
  if (MIC > 0._r8) then
  f_MB_atm = f_LM_MB * (this%cue_ref - this%cue_t * (forc_st - this%tae_ref))
  endif
!Equation 1
  POM = POM + forc_npp * this%param_pi + f_AG_break * this%param_pa - f_PO_AG - f_PO_LM

!Equation 8
  LMWC = LMWC + forc_npp * (1._r8 - this%param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + &
  f_MB_turn * (1._r8 - param_pb) + f_MA_LM

!Equation 19
  AGG = AGG + f_MA_AG + f_PO_AG - f_AG_break

!Equation 21
  MAOM = MAOM + f_LM_MB - f_MA_LM + f_MB_turn * param_pb - f_MA_AG + f_AG_break * (1._r8 - this%param_pa)
  
!Equation 24
  MIC = MIC + f_LM_MB - f_MB_turn - f_MB_atm

!Equation 27
!WORKING
  
end subroutine decomp
  ! decomposition subroutine end

! hydrological properties start
subroutine soilwater(this, forc_sw, scalar_wd, scalar_wb)
implicit none
  integer,parameter   :: r8 = selected_real_kind(12)    ! 8 byte real
  class(soil_type), intent(in) :: this
  real(r8), intent(in)    :: forc_sw    
  real(r8), intent(out)   :: scalar_wd
  real(r8), intent(out)   :: scalar_wb

 !Equation 5
 scalar_wd = (forc_sw/this%porosity)**0.5_r8

 !Equation 4
 scalar_wb = exp(this%lambda * this%matpot) * (this%kamin + (1._r8 - this%kamin) * &
  ((this%porosity - forc_sw) / this%porosity)**0.5_r8) * scalar_wd

end subroutine soilwater
! hydrological properties end

end module soilType

! main program start
PROGRAM Millennial
! History
! Xiaofeng Xu created this code program to play with Millennial model structure (ICOS workshop Mar 14-16, 2016 in Boulder, CO)
! The code is created in May - June 2016, solely by Xiaofeng XU (xxu@mail.sdsu.edu)
! This is a toy verion of the Millennial model (C only version, N P will be added in future updates)
  
  use soilType
  implicit none

  type(soil_type) :: this
  integer,parameter     :: r8 = selected_real_kind(12)  !8 byte real
  integer               :: nr, i, n, flag_output, flag_annual     
  
 !  for initile file
  character(len = 256) :: initialfile
  character(len = 256) :: soilparafile

! for output file
  character(len = 256) :: outputfile
! end of defining output file

! the input data: driving forces
  real(r8), dimension(:), allocatable :: forc_st
  real(r8), dimension(:), allocatable :: forc_sw ! volumetric water content mm3/mm3
  real(r8), dimension(:), allocatable :: forc_npp
! end of driving forces

!!  key variables to track the system over time
! pools 
  real(r8), dimension(:), allocatable :: LMWC
  real(r8), dimension(:), allocatable :: POM
  real(r8), dimension(:), allocatable :: MIC
  real(r8), dimension(:), allocatable :: MAOM
  real(r8), dimension(:), allocatable :: AGG
! end of pools

! flux
  real(r8), dimension(:), allocatable :: f_AG_break     ! AGG -> POM + MAOM;  Aggregate breakdown
  real(r8), dimension(:), allocatable :: f_PO_AG        ! AGG -> POM;         Aggregation
  real(r8), dimension(:), allocatable :: f_MA_AG        ! AGG -> MAOM;        Aggregation
  real(r8), dimension(:), allocatable :: f_PO_LM        ! POM -> LMWC;        Decomposition
  real(r8), dimension(:), allocatable :: f_LM_leach     ! LMWC -> out;        Leaching
  real(r8), dimension(:), allocatable :: f_LM_MA        ! LMWC -> MAOM;       Adsorption
  real(r8), dimension(:), allocatable :: f_LM_MB        ! LMWC -> MIC;        Uptake
  real(r8), dimension(:), allocatable :: f_MB_turn      ! MIC -> LMWC + MAOM; Microbial turnover
  real(r8), dimension(:), allocatable :: f_MA_LM        ! MAOM -> LMWC;       Desoprtion
  real(r8), dimension(:), allocatable :: f_MB_atm       ! MIC -> out;         Microbial maintenance
! end of flux

! time-dependent parameters
  real(r8), dimension(:), allocatable :: scalar_wd
  real(r8), dimension(:), allocatable :: scalar_wb
  real(r8), dimension(:), allocatable :: kaff_eq
  real(r8), dimension(:), allocatable :: param_qmax
  real(r8), dimension(:), allocatable :: param_pb
  real(r8), dimension(:), allocatable :: vmax_pl
  real(r8), dimension(:), allocatable :: vmax_lb
! end of time-dependent parameters

  real(r8)      :: initial_pom
  real(r8)      :: initial_lmwc
  real(r8)      :: initial_mic
  real(r8)      :: initial_maom
  real(r8)      :: initial_agg
! end of key variables

  integer, parameter    :: soil_par_num = 28
! character(len=256)      :: soil_par_f = './soilpara_in'   ! local file name
  integer           :: ier                      ! error code
  character(len=40)     :: soil_par_name(soil_par_num)  ! parameter name
  real(r8)          :: dummy(soil_par_num)
  
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

  allocate(scalar_wd(1:nr)) ; scalar_wd(1:nr) = 0._r8
  allocate(scalar_wb(1:nr)) ; scalar_wb(1:nr) = 0._r8
  allocate(kaff_eq(1:nr))   ; kaff_eq(1:nr) = 0._r8
  allocate(param_qmax(1:nr)); param_qmax(1:nr) = 0._r8
  allocate(param_pb(1:nr))  ; param_pb(1:nr) = 0._r8
  allocate(vmax_pl(1:nr))   ; vmax_pl(1:nr) = 0._r8
  allocate(vmax_lb(1:nr))   ; vmax_lb(1:nr) = 0._r8
  
  allocate(LMWC(1:nr))
  allocate(POM(1:nr))
  allocate(MIC(1:nr))
  allocate(MAOM(1:nr))
  allocate(AGG(1:nr))
  
  allocate(f_AG_break(1:nr))
  allocate(f_PO_AG (1:nr))
  allocate(f_MA_AG(1:nr))
  allocate(f_PO_LM (1:nr))
  allocate(f_LM_leach(1:nr))
  allocate(f_LM_MA(1:nr))
  allocate(f_LM_MB (1:nr))
  allocate(f_MB_turn(1:nr))
  allocate(f_MA_LM  (1:nr))
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
  this%param_pi       = dummy(i); i = i + 1
  this%param_pa       = dummy(i); i = i + 1
  this%kaff_pl        = dummy(i); i = i + 1
  this%alpha_pl       = dummy(i); i = i + 1
  this%eact_pl        = dummy(i); i = i + 1
  this%rate_pa        = dummy(i); i = i + 1
  this%rate_break     = dummy(i); i = i + 1
  this%rate_leach     = dummy(i); i = i + 1
  this%param_pad      = dummy(i); i = i + 1
  this%param_p1       = dummy(i); i = i + 1
  this%param_p2       = dummy(i); i = i + 1
  this%param_pH       = dummy(i); i = i + 1
  this%param_bulkd    = dummy(i); i = i + 1
  this%param_c1       = dummy(i); i = i + 1
  this%param_c2       = dummy(i); i = i + 1
  this%param_clay     = dummy(i); i = i + 1
  this%kaff_lb        = dummy(i); i = i + 1
  this%alpha_lb       = dummy(i); i = i + 1
  this%eact_lb        = dummy(i); i = i + 1
  this%rate_bd        = dummy(i); i = i + 1
  this%rate_ma        = dummy(i); i = i + 1
  this%cue_ref        = dummy(i); i = i + 1
  this%cue_t          = dummy(i); i = i + 1
  this%tae_ref        = dummy(i); i = i + 1
  this%matpot         = dummy(i); i = i + 1
  this%lambda         = dummy(i); i = i + 1
  this%porosity       = dummy(i); i = i + 1
  this%kamin          = dummy(i)
  this%gas_const      = 8.31446_r8 ! Universal gas constant (J/K/mol)

  write(*,*) "Model inializing! "
  write(*,*) "please enter the name of file for initilizing the model"
  read(*,*) initialfile
    
  open(unit = 11, file=initialfile)
    
  read (11,*,iostat=ier) initial_pom, initial_lmwc, initial_mic, initial_maom, initial_agg

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
  MIC(1)=initial_mic
  MAOM(1)=initial_maom
  AGG(1)=initial_agg

  do n = 1, nr

!These should be outside the loop and written as parameters, but might want to leave option of time dependence
!Equation 12
  !Mayes 2012, SSAJ
  kaff_eq(n) =  exp(-this%param_p1 * this%param_pH - this%param_p2)

!Equation 13
  !Mayes 2012, SSAJ
  param_qmax(n) = this%param_bulkd * exp(this%param_c1 * log(this%param_clay) + this%param_c2)

call soilwater(this, forc_sw(n), scalar_wd(n), scalar_wb(n))

call decomp(this, forc_st(n), forc_npp(n), &
    scalar_wd(n), scalar_wb(n), kaff_eq(n), param_qmax(n), param_pb(n), vmax_pl(n), vmax_lb(n), &
    LMWC(n), POM(n), MIC(n), MAOM(n), AGG(n), f_AG_break(n), f_PO_AG(n),&
    f_MA_AG(n), f_PO_LM(n), f_LM_leach(n), f_LM_MA(n), f_LM_MB(n), f_MB_turn(n),&
    f_MA_LM(n), f_MB_atm(n))
 
! updating the pool after each iteration 
  if(n < nr) then
  LMWC(n+1)=LMWC(n)
  POM(n+1)=POM(n)
  MIC(n+1)=MIC(n)
  MAOM(n+1)=MAOM(n)
  AGG(n+1)=AGG(n)
  !print *, 'LMWC:', LMWC(n), ' POM:', POM(n), ' MIC:', MIC(n), ' MAOM:', MAOM(n), ' AGG:', AGG(n)
  endif

  print *, 'finished time step', n 

  end do
  
  if(flag_output ==1) then
call writeoutput(flag_annual, nr, &
    scalar_wd, scalar_wb, kaff_eq, param_qmax, param_pb, vmax_pl, vmax_lb, &
    LMWC, POM, MIC, MAOM, AGG, f_AG_break, f_PO_AG,&
    f_MA_AG, f_PO_LM, f_LM_leach, f_LM_MA, f_LM_MB, f_MB_turn,&
    f_MA_LM, f_MB_atm, outputfile)
  end if

    print *, 'write output finished !'

  deallocate(forc_st)
  deallocate(forc_sw)
  deallocate(forc_npp)
  
  deallocate(LMWC)
  deallocate(POM)
  deallocate(MIC)
  deallocate(MAOM)
  deallocate(AGG)
  
  deallocate(f_AG_break)
  deallocate(f_PO_AG)
  deallocate(f_MA_AG)
  deallocate(f_PO_LM)
  deallocate(f_LM_leach)
  deallocate(f_LM_MA)
  deallocate(f_LM_MB)
  deallocate(f_MB_turn)
  deallocate(f_MA_LM)
  deallocate(f_MB_atm)

  deallocate(scalar_wd)
  deallocate(scalar_wb)
  deallocate(kaff_eq)
  deallocate(param_qmax)
  deallocate(param_pb)
  deallocate(vmax_pl)
  deallocate(vmax_lb)
! end of the allocation
  
  stop
END PROGRAM Millennial
! main program end