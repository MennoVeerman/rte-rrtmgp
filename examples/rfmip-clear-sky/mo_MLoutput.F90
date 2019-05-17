!Write gas concetrations, water vapour content, temperature, pressure
!and optical properties to use as training/evaluation data for Neural
!Networks !
!
module mo_MLoutput 
use mo_gas_concentrations, only: ty_gas_concs
use mo_optical_props,      only: ty_optical_props,ty_optical_props_arry
use mo_source_functions,only: ty_source_func_lw
use mo_gas_optics,      only: ty_gas_optics
use mo_rte_kind, only: wp
implicit none
public :: init_MLoutput,write_MLoutput_tau
character(len=11) :: FilenameT = "outputT.txt"
character(len=11) :: FilenameP = "outputP.txt"

logical::write_gases = .false.

contains


subroutine init_MLoutput(gas_names,ngas)
  character(len=*), dimension(:), intent(in) :: gas_names
  character(len=7), dimension(:), allocatable:: tau_names,lay_names,lvi_names,lvd_names
  character(6) :: string1,string2,string3,string4
  integer, intent(in) :: ngas
  integer :: nbands = 256,iband
  allocate(tau_names(nbands))
  allocate(lay_names(nbands))
  allocate(lvi_names(nbands))
  allocate(lvd_names(nbands))

  do iband = 1,nbands
    write(string1,'(A,I3.3)') 'tau',iband
    write(string2,'(A,I3.3)') 'lay',iband
    write(string3,'(A,I3.3)') 'lvi',iband
    write(string4,'(A,I3.3)') 'lvd',iband

    tau_names(iband) = string1    
    lay_names(iband) = string2
    lvi_names(iband) = string3
    lvd_names(iband) = string4

  enddo
  open(666,file=FilenameT,status='replace')
  if (write_gases) then
    write(666,'(17A8)',advance='no') gas_names(2:ngas)
  endif
  write(666,'(261A8)') gas_names(1),gas_names(2),gas_names(3),'p_lay','t_lay',tau_names(:)

  close(666)

  open(667,file=FilenameP,status='replace')
  if (write_gases) then
    write(667,'(17A8)',advance='no') gas_names(2:ngas)
  endif
  write(667,'(775A8)') gas_names(1),gas_names(2),gas_names(3),'p_lay','t_lay','t_levB','t_levT',lay_names(:),lvi_names(:),lvd_names(:)

  close(667)

end subroutine init_MLoutput


subroutine write_MLoutput_tau(ngas,ncol,nlay,gas_desc,gas_names,&
                              p_lay,p_lev,t_lay,t_lev,optical_props)
  integer, intent(in) :: ngas,ncol,nlay
  integer :: igas,wi,wj
  class(ty_optical_props_arry),intent(in) :: optical_props
  real(wp),allocatable,dimension(:,:,:) :: vmr
  real(wp),dimension(ncol,nlay),intent(in) :: p_lay,t_lay
  real(wp),dimension(ncol,nlay+1),intent(in) :: p_lev,t_lev
  character(len=*), dimension(:), intent(in) :: gas_names
  character(len=128)          :: error_msg
  type(ty_gas_concs), intent(in) :: gas_desc
  allocate(vmr(ncol,nlay,ngas))
  do igas = 1,ngas
    error_msg = gas_desc%get_vmr(gas_names(igas),vmr(:,:,igas))
  enddo
  
  open(666,file=FilenameT,position='append')
  do wi = 1,ncol  
    do wj = 1,nlay
       if (write_gases) then
           write(666,'(17E15.7)',advance='no') vmr(wi,wj,2:ngas)
       endif
           write(666,'(261E15.7)')vmr(wi,wj,1),vmr(wi,wj,2),vmr(wi,wj,3),p_lay(wi,wj),&
                      t_lay(wi,wj),optical_props%tau(wi,wj,:)/ (p_lev(wi,wj+1)-p_lev(wi,wj))
      

     enddo
  enddo
  close(666)
  deallocate(vmr)
end subroutine write_MLoutput_tau

subroutine write_MLoutput_planck(ngas,ncol,nlay,gas_desc,gas_names,&
                              p_lay,t_lay,t_lev,source)
  integer, intent(in) :: ngas,ncol,nlay
  integer :: igas,wi,wj
  class(ty_source_func_lw),intent(in) :: source
  real(wp),allocatable,dimension(:,:,:) :: vmr
  real(wp),dimension(ncol,nlay),intent(in) :: p_lay,t_lay
  real(wp),dimension(ncol,nlay+1),intent(in) :: t_lev
  character(len=*), dimension(:), intent(in) :: gas_names
  character(len=128)          :: error_msg
  type(ty_gas_concs), intent(in) :: gas_desc
  allocate(vmr(ncol,nlay,ngas))
  do igas = 1,ngas
    error_msg = gas_desc%get_vmr(gas_names(igas),vmr(:,:,igas))
  enddo

  open(666,file=FilenameP,position='append')
  do wi = 1,ncol
    do wj = 1,nlay
       if (write_gases) then
         write(666,'(17E15.7)',advance='no') vmr(wi,wj,2:ngas)
       endif
         write(666,'(775E15.7)')vmr(wi,wj,1),vmr(wi,wj,2),vmr(wi,wj,3),p_lay(wi,wj),&
                      t_lay(wi,wj), t_lev(wi,wj), t_lev(wi,wj+1), source%lay_source(wi,wj,:),&
                      source%lev_source_inc(wi,wj,:),source%lev_source_dec(wi,wj,:)


   enddo
  enddo
  close(666)
  deallocate(vmr)
end subroutine write_MLoutput_planck



end module

































