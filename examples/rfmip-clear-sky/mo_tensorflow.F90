module mo_tensorflow
use mo_rte_kind, only: wp
use mo_gas_concentrations, only: ty_gas_concs
use mo_optical_props,      only: ty_optical_props,ty_optical_props_arry

implicit none

public :: write_predict!tau_pred,python_tensorflow

character(len=19) :: prednameT   = 'predictdata_tau.csv'
character(len=22) :: prednameP   = 'predictdata_planck.csv'

character(len=17) :: ptauname   = 'predicted_tau.txt'
character(len=20) :: planckname = 'predicted_planck.txt'
character(len=10) :: pythpath   = 'pythonpath'
character(len=10) :: tenspath   = 'tensorpath'

contains

subroutine write_predict(ngas,ncol,nlay,gas_desc,gas_names,&
                        p_lay,t_lay,p_lev,t_lev,optical_props)
integer                                    :: igas,icol,ilay
integer,intent(in)                         :: ngas,ncol,nlay
real(wp),allocatable,dimension(:,:,:)      :: vmr
real(wp),dimension(ncol,nlay),intent(in)   :: p_lay,t_lay
real(wp),dimension(ncol,nlay+1),intent(in) :: p_lev,t_lev
character(len=*), dimension(:), intent(in) :: gas_names
character(len=128)          :: error_msg
type(ty_gas_concs), intent(in) :: gas_desc
class(ty_optical_props_arry),intent(inout) :: optical_props


allocate(vmr(ncol,nlay,ngas))
do igas = 1,ngas
    error_msg = gas_desc%get_vmr(gas_names(igas),vmr(:,:,igas))
enddo

open(333,file=prednameT,status='replace')
do ilay =  1,nlay
    do icol = 1,ncol
        write(333,'(5(E15.7))') vmr(icol,ilay,1),vmr(icol,ilay,2),vmr(icol,ilay,3),p_lay(icol,ilay),t_lay(icol,ilay)
    enddo
enddo
close(333)

open(333,file=prednameP,status='replace')
do ilay =  1,nlay
    do icol = 1,ncol
        write(333,'(5(E15.7))') vmr(icol,ilay,1),vmr(icol,ilay,2),vmr(icol,ilay,3),p_lay(icol,ilay),t_lay(icol,ilay), t_lev(icol,ilay), t_lev(icol, ilay+1)
    enddo
enddo
close(333)


call python_tensorflow
call read_predict_tau(ncol,nlay,p_lev,p_lay,optical_props)
!call read_predict_planck(ncol,nlay,p_lev,p_lay,optical_props)
end subroutine write_predict

subroutine python_tensorflow()
call system("source ~/virtualenv/firstCNN_intel_CPU/bin/activate")
call system("~/virtualenv/firstCNN_intel_CPU/bin/python3 ~/rte-rrtmgp/examples/rfmip-clear-sky/neuralnetwork/tau_predict.py")
!call system("~/virtualenv/firstCNN_intel_CPU/bin/python3 ~/rte-rrtmgp/examples/rfmip-clear-sky/neuralnetwork/plk_predict.py")
call system("~/virtualenv/firstCNN_intel_CPU/bin/deactivate")
end subroutine python_tensorflow

subroutine read_predict_tau(ncol,nlay,p_lev,p_lay,optical_props)
real(wp),allocatable,dimension(:,:)        :: tau_pred,tautemp
real(wp),dimension(ncol,nlay+1),intent(in) :: p_lev,p_lay
integer,intent(in)                         :: ncol,nlay
integer                                    :: icol,ilay,N,Nn,itau
integer                                    :: ntau = 256, nfeat = 6
class(ty_optical_props_arry),intent(inout) :: optical_props
real(wp),dimension(262)                    :: M0,M1,S0,S1,M,S
logical                                    :: isparallel
isparallel = .true.
allocate(tau_pred(ncol*nlay,ntau))
allocate(tautemp(ntau,ncol*nlay))

open(666,file=ptauname,status='old')
read(666,*) tautemp
tau_pred = transpose(tautemp)
close(666)

if (isparallel == .true.) then
    open(666,file='means0_tau.txt',status='old')
    read(666,*) M0
    close(666)
    open(666,file='means1_tau.txt',status='old')
    read(666,*) M1
    close(666)
    open(666,file='stdev0_tau.txt',status='old')
    read(666,*) S0
    close(666)
    open(666,file='stdev1_tau.txt',status='old')
    read(666,*) S1
    close(666)
else if (isparallel == .false.) then
    open(666,file='means_tau.txt',status='old')
    read(666,*) M
    close(666)
    open(666,file='stdev_tau.txt',status='old')
    read(666,*) S
    close(666)
endif
do itau = 1,ntau  
    N = 1
    do ilay =  1,nlay
        do icol = 1,ncol 
            if (isparallel == .true. ) then
                if (p_lay(icol,ilay) < 9948.43156 .and. .not. any((/177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,&
                                                                    225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240/)==itau)) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S0(itau+nfeat)*tau_pred(N,itau)+M0(itau+nfeat)) * (p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)             
                else if (p_lay(icol,ilay) >= 9948.43156) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S1(itau+nfeat)*tau_pred(N,itau)+M1(itau+nfeat)) * (p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)
                endif           
            else if (isparallel ==.false.) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S(itau+nfeat)*tau_pred(N,itau)+M(itau+nfeat)) *(p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)
            endif
            N = N + 1
        enddo
    enddo
enddo
deallocate(tau_pred)
end subroutine read_predict_tau

subroutine read_predict_planck(ncol,nlay,p_lev,p_lay,optical_props)
real(wp),allocatable,dimension(:,:)        :: tau_pred,tautemp
real(wp),dimension(ncol,nlay+1),intent(in) :: p_lev,p_lay
integer,intent(in)                         :: ncol,nlay
integer                                    :: icol,ilay,N,Nn,itau
integer                                    :: ntau = 256, nfeat = 6
class(ty_optical_props_arry),intent(inout) :: optical_props
real(wp),dimension(262)                    :: M0,M1,S0,S1,M,S
logical                                    :: isparallel
isparallel = .true.
allocate(tau_pred(ncol*nlay,ntau))
allocate(tautemp(ntau,ncol*nlay))

open(666,file=planckname,status='old')
read(666,*) tautemp
tau_pred = transpose(tautemp)
close(666)

if (isparallel == .true.) then
    open(666,file='means0_planck.txt',status='old')
    read(666,*) M0
    close(666)
    open(666,file='means1_planck.txt',status='old')
    read(666,*) M1
    close(666)
    open(666,file='stdev0_planck.txt',status='old')
    read(666,*) S0
    close(666)
    open(666,file='stdev1_planck.txt',status='old')
    read(666,*) S1
    close(666)
else if (isparallel == .false.) then
    open(666,file='means_planck.txt',status='old')
    read(666,*) M
    close(666)
    open(666,file='stdev_planck.txt',status='old')
    read(666,*) S
    close(666)
endif

do itau = 1,ntau
    N = 1
    do ilay =  1,nlay
        do icol = 1,ncol
            if (isparallel == .true. .and. icol > 50) then
                if (p_lay(icol,ilay) < 9948.43156 .and. .not. any((/177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,&
                                                                    225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240/)==itau)) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S0(itau+nfeat)*tau_pred(N,itau)+M0(itau+nfeat)) * (p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)
                else if (p_lay(icol,ilay) >= 9948.43156) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S1(itau+nfeat)*tau_pred(N,itau)+M1(itau+nfeat))*(p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)
                endif

           else if (isparallel ==.false.) then
                    optical_props%tau(icol,ilay,itau) = max(exp(S(itau+nfeat)*tau_pred(N,itau)+M(itau+nfeat))*(p_lev(icol,ilay+1)-p_lev(icol,ilay)) , 0.)
            endif
            N = N + 1
        enddo
    enddo
enddo
deallocate(tau_pred)
end subroutine read_predict_planck

end module

