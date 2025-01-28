module datm_datamode_atmmesh_mod

  use ESMF             , only : ESMF_State, ESMF_LOGMSG_INFO, ESMF_LogWrite, ESMF_SUCCESS
  use NUOPC            , only : NUOPC_Advertise
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod    , only : shr_const_TkFrz, shr_const_pi, shr_const_ocn_ref_sal
  use shr_sys_mod      , only : shr_sys_abort
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, chkerr
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write
  use dshr_strdata_mod , only : shr_strdata_type

  implicit none
  private ! except

  public  :: datm_datamode_atmmesh_advertise
  public  :: datm_datamode_atmmesh_init_pointers
  public  :: datm_datamode_atmmesh_advance

  ! export fields
  real(r8), pointer :: Sa_pslv   (:) => null() ! sea level pressure
  real(r8), pointer :: Sa_u10m   (:) => null() ! zonal wind height 10m
  real(r8), pointer :: Sa_v10m   (:) => null() ! meridional wind height 10m
  real(r8), pointer :: Sa_t2m    (:) => null() ! surface temperature height 2m
  real(r8), pointer :: Sa_q2m    (:) => null() ! surface humidity height 2m
  real(r8), pointer :: Faxa_lwdn (:) => null() ! downward longwave radiation
  real(r8), pointer :: Faxa_swnet(:) => null() ! net shortwave radiation
  real(r8), pointer :: Faxa_rain (:) => null() ! total precipitation

  character(*) , parameter :: nullstr = 'null'
  character(*) , parameter :: rpfile  = 'rpointer.atm'
  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine datm_datamode_atmmesh_advertise(exportState, fldsexport, flds_scalar_name, rc)

    ! input/output variables
    type(esmf_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Advertise export fields
    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'Sa_pslv')
    call dshr_fldList_add(fldsExport, 'Sa_u10m')
    call dshr_fldList_add(fldsExport, 'Sa_v10m')
    call dshr_fldList_add(fldsExport, 'Sa_t2m' )
    call dshr_fldList_add(fldsExport, 'Sa_q2m' )
    call dshr_fldList_add(fldsExport, 'Faxa_lwdn')
    call dshr_fldList_add(fldsExport, 'Faxa_swnet')
    call dshr_fldList_add(fldsExport, 'Faxa_rain')

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(datm_comp_advertise): Fr_atm'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine datm_datamode_atmmesh_advertise

  !===============================================================================
  subroutine datm_datamode_atmmesh_init_pointers(exportState, sdat, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(shr_strdata_type) , intent(in)    :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname='(datm_init_pointers): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize pointers to export fields
    call dshr_state_getfldptr(exportState, 'Sa_pslv'   , fldptr1=Sa_pslv   , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_u10m'   , fldptr1=Sa_u10m   , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_v10m'   , fldptr1=Sa_v10m   , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_t2m'    , fldptr1=Sa_t2m    , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_q2m'    , fldptr1=Sa_q2m    , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_lwdn' , fldptr1=Faxa_lwdn , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swnet', fldptr1=Faxa_swnet, allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_rain' , fldptr1=Faxa_rain , allowNullReturn=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (associated(Sa_pslv   )) Sa_pslv   (:) = 0.0_r8
    if (associated(Sa_u10m   )) Sa_u10m   (:) = 0.0_r8
    if (associated(Sa_v10m   )) Sa_v10m   (:) = 0.0_r8
    if (associated(Sa_t2m    )) Sa_t2m    (:) = 0.0_r8
    if (associated(Sa_q2m    )) Sa_q2m    (:) = 0.0_r8
    if (associated(Faxa_lwdn )) Faxa_lwdn (:) = 0.0_r8
    if (associated(Faxa_swnet)) Faxa_swnet(:) = 0.0_r8
    if (associated(Faxa_rain )) Faxa_rain (:) = 0.0_r8

  end subroutine datm_datamode_atmmesh_init_pointers

  !===============================================================================
  subroutine datm_datamode_atmmesh_advance(mainproc, logunit, mpicom, rc)

    ! input/output variables
    logical                , intent(in)    :: mainproc
    integer                , intent(in)    :: logunit
    integer                , intent(in)    :: mpicom
    integer                , intent(out)   :: rc

    ! local variables
    integer  :: n                   ! indices
    integer  :: lsize               ! size of attr vect
    character(len=*), parameter :: subname='(datm_datamode_atmmesh_advance): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !lsize = size(Sa_pslv)

    !do n = 1,lsize
    !   !--- bottom layer height ---
    !   if (.not. associated(strm_z)) Sa_z(n) = 30.0_r8
    !end do

  end subroutine datm_datamode_atmmesh_advance

end module datm_datamode_atmmesh_mod
