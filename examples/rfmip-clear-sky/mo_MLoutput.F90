!Write gas concetrations, water vapour content, temperature, pressure
!and optical properties to use as training/evaluation data for Neural
!Networks !
!
module mo_MLoutput 
use mo_gas_concentrations, only: ty_gas_concs
use mo_optical_props,      only: ty_optical_props,ty_optical_props_arry
use mo_gas_optics,      only: ty_gas_optics
use mo_rte_kind, only: wp
implicit none
public :: init_MLoutput,write_MLoutput_tau
character(len=10) :: Filename = "output.txt"


contains


subroutine init_MLoutput(gas_names,ngas)
  character(len=*), dimension(:), intent(in) :: gas_names
  integer, intent(in) :: ngas
  open(666,file=Filename,status='replace')
  write(666,'(21A8)') gas_names(:ngas),'p_lay','t_lay','Taus -->'
  close(666)
end subroutine init_MLoutput


subroutine write_MLoutput_tau(ngas,ncol,nlay,gas_desc,gas_names,&
                              p_lay,t_lay,optical_props)
  integer, intent(in) :: ngas,ncol,nlay
  integer :: igas,wi,wj
  class(ty_optical_props_arry),intent(in) :: optical_props
  real(wp),allocatable,dimension(:,:,:) :: vmr
  real(wp),dimension(ncol,nlay),intent(in) :: p_lay,t_lay
  character(len=*), dimension(:), intent(in) :: gas_names
  character(len=128)          :: error_msg
  type(ty_gas_concs), intent(in) :: gas_desc
  allocate(vmr(ncol,nlay,ngas))
  do igas = 1,ngas
    error_msg = gas_desc%get_vmr(gas_names(igas),vmr(:,:,igas))
  enddo
  
  open(666,file=Filename,position='append')
  do wi = 1,ncol
    do wj = 1,nlay
       write(666,'(276E15.7)') vmr(wi,wj,:),p_lay(wi,wj),t_lay(wi,wj),optical_props%tau(wi,wj,:)
    enddo
  enddo
  close(666)

  deallocate(vmr)
end subroutine write_MLoutput_tau




end module

































