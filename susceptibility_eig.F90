#include "fintrf.h"
!======================================================================
!
#if 0
!
!     convec.F
!     .F file needs to be preprocessed to generate .for equivalent
!
#endif

subroutine convec(Nk,Nb,Nf,q,kv,fr,lambda,phsus,ppsus)
    implicit none
    integer(kind=8) ::  Nk,Nb,Nf
    integer(kind=4) :: o1,o2,o3,o4,b,bp,m,n
    integer(kind=4) :: ik

    real(kind=8) :: q(2),lambda,fr(Nf,2),kv(Nk,3),mu=0.d0
    real(kind=8) :: e1,e2,deltaE=1d-8;
    real*8 vk(2),wk

    ! complex(kind=8) :: Hpk(Nb,Nb),Hmk(Nb,Nb),Hkq(Nb,Nb)
    complex(kind=8) :: uk(Nb,Nb),umk(Nb,Nb),ukq(Nb,Nb)
    real(kind=8) :: ek(Nb),emk(Nb),ekq(Nb)
    complex(kind=8) :: fk(Nf)
    complex(kind=8) :: wfcoe,Lindhard,tmp1,tmp2

    complex(kind=8) xph(Nb**2,Nb**2)
    complex(kind=8) xpp(Nb**2,Nb**2)
    complex(kind=8) phsus(Nf*Nb**2,Nf*Nb**2)
    complex(kind=8) ppsus(Nf*Nb**2,Nf*Nb**2)
    
    real(kind=8),external :: fermi,dfermi
    ! complex(kind=8),external :: hk

    character*100 strout
    complex(kind=8) , parameter :: one=cmplx(0.d0,1.d0)
    real(kind=8) , parameter ::pi=asin(1.d0)*2
    
    phsus=cmplx(0.00,0.00)
    ppsus=cmplx(0.00,0.00)

    ! write(strout,*) "Nb:",Nb
    ! call mexPrintf(strout//achar(10))

    do ik=1,Nk
        vk=kv(ik,1:2); wk=kv(ik,3)
        ! uk=hk(vk,mu)
        ! umk=hk(-vk,mu)
        ! ukq=hk(vk+q,mu)

        ! write(strout,*) "flag1"
        ! call mexPrintf(strout//achar(10))

        call hk(vk,mu,uk)
        call hk(-vk,mu,umk)
        call hk(vk+q,mu,ukq)

        ! write(strout,*) "flag2"
        ! call mexPrintf(strout//achar(10))

        call ZHEIGEN(Nb,uk,ek)
        call ZHEIGEN(Nb,umk,emk)
        call ZHEIGEN(Nb,ukq,ekq)

        ! if (ik==3) then
        !     write(strout,*) "x1:",real(ek(1))
        !     call mexPrintf(strout//achar(10))
        !     write(strout,*) "x2:",real(ek(2))
        !     call mexPrintf(strout//achar(10))
        !     write(strout,*) "x3:",real(ek(3))
        !     call mexPrintf(strout//achar(10))
        ! endif


        ! write(strout,*) "flag3"
        ! call mexPrintf(strout//achar(10))

        xph=cmplx(0.00,0.00)
        xpp=cmplx(0.00,0.00)
        do o1=1,Nb
            do o2=1,Nb
                do o3=1,Nb
                    do o4=1,Nb
                        do b=1,Nb
                            do bp=1,Nb
                                
                                !ph susceptibility
                                wfcoe=ukq(o1,b) *conjg(ukq(o3,b)) &
                                &    *uk(o4,bp)*conjg(uk(o2,bp))
                                e1=real(ekq(b))
                                e2=real(ek(bp))
                                if(abs( e1-e2)> deltaE) then
                                    Lindhard=(fermi(e1,lambda)-fermi(e2,lambda))/(e1-e2)
                                else
                                    Lindhard=dfermi((e1+e2)/2,lambda)
                                end if
                                xph(o1+(o2-1)*Nb,o3+(o4-1)*Nb)= xph(o1+(o2-1)*Nb,o3+(o4-1)*Nb) + wfcoe*Lindhard

                                !pp susceptibility
                                wfcoe=ukq(o1,b) *conjg(ukq(o3,b)) &
                                &    *umk(o2,bp)*conjg(umk(o4,bp))
                                e1=-real(ekq(b))
                                e2=real(emk(bp))
                                if(abs( e1-e2)> deltaE) then
                                    Lindhard=(fermi(e1,lambda)-fermi(e2,lambda))/(e1-e2)
                                else
                                    Lindhard=dfermi((e1+e2)/2,lambda)
                                end if
                                xpp(o1+(o2-1)*Nb,o3+(o4-1)*Nb)= xpp(o1+(o2-1)*Nb,o3+(o4-1)*Nb) - wfcoe*Lindhard

                            end do
                        end do
                    end do
                end do 
            end do
        end do

        do m=1,Nf
            fk(m)=exp(one*( fr(m,1)*vk(1)+fr(m,2)*vk(2) ))
        end do

        do m=1,Nf
            do n=1,Nf
                do o1=1,Nb**2
                    do o2=1,Nb**2
                        phsus(m+(o1-1)*Nf,n+(o2-1)*Nf)=phsus(m+(o1-1)*Nf,n+(o2-1)*Nf) &
                        &         +fk(m)*xph(o1,o2)*conjg(fk(n))*wk
                        ppsus(m+(o1-1)*Nf,n+(o2-1)*Nf)=ppsus(m+(o1-1)*Nf,n+(o2-1)*Nf) &
                        &         +fk(m)*xpp(o1,o2)*conjg(fk(n))*wk
                    end do
                end do
            end do
        end do

    end do
    
        


    return
    end



function fermi(e,T)
    implicit none
    real(kind=8) :: T
    real(kind=8) :: fermi,e

    fermi= e/(4* (cosh(e/2/T))**2 * T**2)
    return
    end

function dfermi(e,T)
    implicit none
    real(kind=8) :: T
    real(kind=8) :: dfermi,e

    dfermi=(T- e* tanh(e/2/T))/(4* (cosh(e/2/T))**2 * T**3)
    return
    end

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    implicit none
    !mexFunction arguments:
    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs

    !Function declarations:
    mwPointer mxGetPr, mxGetPi
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    !Array information:
    mwPointer mx, my, nx , ny
    mwSize row
    mwSize kv_size,fr_size,Nchi

    real(kind=8) ::lambda
    integer(kind=8) :: Nk,Nb,Nf

    real(kind=8) :: q(2)
    real(kind=8),allocatable :: kv(:,:),fr(:,:)
    complex(kind=8),allocatable :: phsus(:,:)
    complex(kind=8),allocatable :: ppsus(:,:)
    character*200 strout1


    !Check for proper number of arguments.
    if (nrhs /= 4) then
        call mexErrMsgIdAndTxt ('MATLAB:convec:nInput', &
    &                           'nine inputs required.')
    elseif (nlhs > 2) then
        call mexErrMsgIdAndTxt ('MATLAB:convec:nOutput', &
    &                           'Too many output arguments.')
    endif

    call mxCopyPtrToReal8(mxGetPr(prhs(1)),q,2)
    call mxCopyPtrToReal8(mxGetPr(prhs(4)),lambda,1)

    mx = mxGetM(prhs(2))
    Nk=int(mx)

    nx=mxGetM(prhs(3))
    Nf=int(nx)

    Nb=3

    kv_size=mx*3
    fr_size=nx*2
    Nchi=Nf*Nb**2

    ! write(strout1,*) "Nk:",Nk
    ! call mexPrintf(strout1//achar(10))
    ! write(strout1,*) "Nf:",Nf
    ! call mexPrintf(strout1//achar(10))
    ! write(strout1,*) "Nb:",Nb
    ! call mexPrintf(strout1//achar(10))
    ! write(strout1,*) "Nchi:",Nchi
    ! call mexPrintf(strout1//achar(10))

    ALLOCATE(kv(Nk,3))
    ALLOCATE(fr(Nf,2))
    ALLOCATE(phsus(Nchi,Nchi))
    ALLOCATE(ppsus(Nchi,Nchi))

    ! call mxCopyPtrToComplex16(mxGetPr(prhs(2)), mxGetPi(prhs(2)),kv,kv_size)
    ! call mxCopyPtrToComplex16(mxGetPr(prhs(3)), mxGetPi(prhs(3)),fr,fr_size)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),kv,kv_size)
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),fr,fr_size)

    plhs(1) = mxCreateDoubleMatrix(Nchi,Nchi, 1)
    plhs(2) = mxCreateDoubleMatrix(Nchi,Nchi, 1)

    call convec(Nk,Nb,Nf,q,kv,fr,lambda,phsus,ppsus)
    
    call mxCopyComplex16ToPtr(phsus,mxGetPr(plhs(1)), mxGetPi(plhs(1)),Nchi**2)
    call mxCopyComplex16ToPtr(ppsus,mxGetPr(plhs(2)), mxGetPi(plhs(2)),Nchi**2)

    ! write(strout1,*) "flag5"
    ! call mexPrintf(strout1//achar(10))

    DEALLOCATE(kv)
    DEALLOCATE(fr)
    DEALLOCATE(phsus)
    DEALLOCATE(ppsus)

    return
    end

subroutine hk(k,mu,H)
    implicit none
    integer(kind=4):: i,j
    real(kind=8) :: k(2) ,mu, k1 ,k2 
    complex(kind=8) :: H(3,3),one=cmplx(0,1)
    H=0
    do i=1,3
        H(i,i)=H(i,i)-mu
    end do
      
    k1=k(1)
    k2=k(1)/2 + sqrt(3.0)/2.0*k(2)
      
    H(1,2)=-(1+exp(-one*k1))
    H(1,3)=-(1+exp(-one*k2))
    H(2,3)=-(1+exp(one*(k1-k2)))
        
    do j=1,3; 
        do i=1,j-1
            H(j,i)=conjg( H(i,j) )
        end do
    end do
    return
end

SUBROUTINE ZHEIGEN(N,A,EVAL)
    IMPLICIT NONE
    INTEGER N,LWORK,INFO
    DOUBLE PRECISION RWORK(N*3),EVAL(N)
    COMPLEX*16 WORK(N*2),A(N,N)
  
    IF(N<=0)STOP 'INVALID N @ ZHEIGEN'
    IF(N==1)THEN; EVAL(1)=real(A(1,1)); A=1; RETURN; END IF
  
    LWORK=N*2
    !CALL CHEEV('V','U',N,A,N,EVAL,WORK,LWORK,RWORK,INFO)
    CALL ZHEEV('V','U',N,A,N,EVAL,WORK,LWORK,RWORK,INFO)
    IF(INFO/=0)STOP 'ZHEIGEN FAILED'
    RETURN
  END SUBROUTINE ZHEIGEN


! function hk(k,mu)
!     implicit none
!     integer(kind=4) :: i,j
!     real(kind=8) :: k(2) ,mu, k1 ,k2 
!     complex(kind=8) :: hk(3,3),one=cmplx(0,1)

!     ! hk=cmplx(0,0)
!     hk=0
!     do i=1,3
!         hk(i,i)=hk(i,i)-mu
!     end do

!     k1=k(1)
!     k2=k(1)/2 + sqrt(3.0)/2.0*k(2)

!     hk(1,2)=-(1+exp(-one*k1))
!     hk(1,3)=-(1+exp(-one*k2))
!     hk(2,3)=-(1+exp(one*(k1-k2)))
    
!     do j=1,3; 
!         do i=1,j-1
!             hk(j,i)=conjg( hk(i,j) )
!         end do
!     end do

!     return
!     end       