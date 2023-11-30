#include "fintrf.h"
!======================================================================
!
#if 0
!
!
#endif

subroutine convec(Nk,Nb,Nf,q,uk,ek,fk,lambda,phsus,ppsus)
    implicit none
    integer(kind=4) ::  Nk,Nb,Nf
    integer(kind=4) :: o1,o2,o3,o4,b,bp,m,n
    integer(kind=4) :: mk1,nk1,mq,nq,mkq,nkq,mmk,nmk

    real(kind=8) :: q(2),lambda,Lindhard
    real(kind=8) :: e1,e2,deltaE=1d-6,eps=1d-16;

    complex(kind=8) :: uk(Nk,Nk,Nb,Nb)
    complex(kind=8) :: ek(Nk,Nk,Nb)
    complex(kind=8) :: fk(Nk,Nk,Nf)
    complex(kind=8) :: wfcoe,tmp1,tmp2

    complex(kind=8) xph(Nb**2,Nb**2)
    complex(kind=8) xpp(Nb**2,Nb**2)
    complex(kind=8) phsus(Nf*Nb**2,Nf*Nb**2)
    complex(kind=8) ppsus(Nf*Nb**2,Nf*Nb**2)
    
    real(kind=8),external :: fermi,dfermi

    character*200 strout1,strout2,strout3,strout4,strout5

    phsus=cmplx(0.00,0.00)
    ppsus=cmplx(0.00,0.00)

    mq=anint(q(1)*Nk)
    nq=anint(q(2)*Nk)
    
    do mk1=1,Nk
        mmk=-mk1+Nk
        if (mmk==0) then
            mmk=Nk
        end if

        mkq=modulo(mk1+mq,Nk)
        if (mkq==0) then
            mkq=Nk
        end if

        do nk1=1,Nk
            nmk=-nk1+Nk
            if (nmk==0) then
                nmk=Nk
            end if
            
            nkq=modulo(nk1+nq,Nk)
            if (nkq==0) then
                nkq=Nk
            end if

            xph=cmplx(0.00,0.00)
            xpp=cmplx(0.00,0.00)
            do o1=1,Nb
                do o2=1,Nb
                    do o3=1,Nb
                        do o4=1,Nb
                            do b=1,Nb
                                do bp=1,Nb
                                    
                                    !ph susceptibility
                                    wfcoe=uk(mkq,nkq,o1,b) *conjg(uk(mkq,nkq,o3,b)) &
                                    &    *uk(mk1,nk1,o4,bp)*conjg(uk(mk1,nk1,o2,bp))
                                    e1=real(ek(mkq,nkq,b))
                                    e2=real(ek(mk1,nk1,bp))
                                    if( abs(e1-e2) > deltaE ) then
                                        Lindhard=(fermi(e1,lambda)-fermi(e2,lambda))/(e1-e2+eps)
                                    else
                                        Lindhard=dfermi((e1+e2)/2,lambda)
                                    end if
                                    xph(o1+(o2-1)*Nb,o3+(o4-1)*Nb)= xph(o1+(o2-1)*Nb,o3+(o4-1)*Nb) + wfcoe*Lindhard

                                    !pp susceptibility
                                    wfcoe=uk(mkq,nkq,o1,b) *conjg(uk(mkq,nkq,o3,b)) &
                                    &    *uk(mmk,nmk,o2,bp)*conjg(uk(mmk,nmk,o4,bp))
                                    e1=-real(ek(mkq,nkq,b))
                                    e2=real(ek(mmk,nmk,bp))
                                    if( abs(e1-e2) > deltaE ) then
                                        Lindhard=(fermi(e1,lambda)-fermi(e2,lambda))/(e1-e2+eps)
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
                do n=1,Nf
                    do o1=1,Nb**2
                        do o2=1,Nb**2
                            phsus(m+(o1-1)*Nf,n+(o2-1)*Nf)=phsus(m+(o1-1)*Nf,n+(o2-1)*Nf) &
                            &         +fk(mk1,nk1,m)*xph(o1,o2)*conjg(fk(mk1,nk1,n))/Nk**2
                            ppsus(m+(o1-1)*Nf,n+(o2-1)*Nf)=ppsus(m+(o1-1)*Nf,n+(o2-1)*Nf) &
                            &         +fk(mk1,nk1,m)*xpp(o1,o2)*conjg(fk(mk1,nk1,n))/Nk**2
                        end do
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
    mwSize uk_size,ek_size,fk_size,Nchi

    real(kind=8) ::lambda
    integer(kind=4) :: Nk,Nb,Nf

    real(kind=8) :: q(2)
    complex(kind=8),allocatable :: uk(:,:,:,:),ek(:,:,:),fk(:,:,:)
    complex(kind=8),allocatable :: phsus(:,:)
    complex(kind=8),allocatable :: ppsus(:,:)
    character*200 strout1
    integer(kind=4) :: o1,o2,o3,o4,b,bp,m,n


    !Check for proper number of arguments.
    if (nrhs /= 5) then
        call mexErrMsgIdAndTxt ('MATLAB:convec:nInput', &
    &                           'nine inputs required.')
    elseif (nlhs > 2) then
        call mexErrMsgIdAndTxt ('MATLAB:convec:nOutput', &
    &                           'Too many output arguments.')
    endif

    call mxCopyPtrToReal8(mxGetPr(prhs(1)),q,2)
    call mxCopyPtrToReal8(mxGetPr(prhs(5)),lambda,1)

    mx = mxGetM(prhs(2))
    my = mxGetN(prhs(2))
    nx=my/mx*1.0;
    Nk=int(mx)
    Nb=int(sqrt(real(nx)))

    ny=mxGetN(prhs(4))
    Nf=int(ny/mx*1.0)

    uk_size=mx*mx*Nb*Nb
    ek_size=mx*mx*Nb
    fk_size=mx*mx*Nf
    Nchi=Nf*Nb**2

    ALLOCATE(uk(mx,mx,Nb,Nb))
    ALLOCATE(ek(mx,mx,Nb))
    ALLOCATE(fk(mx,mx,Nf))
    ALLOCATE(phsus(Nchi,Nchi))
    ALLOCATE(ppsus(Nchi,Nchi))

    call mxCopyPtrToComplex16(mxGetPr(prhs(2)), mxGetPi(prhs(2)),uk,uk_size)
    call mxCopyPtrToComplex16(mxGetPr(prhs(3)), mxGetPi(prhs(3)),ek,ek_size)
    call mxCopyPtrToComplex16(mxGetPr(prhs(4)), mxGetPi(prhs(4)),fk,fk_size)

    plhs(1) = mxCreateDoubleMatrix(Nchi,Nchi, 1)
    plhs(2) = mxCreateDoubleMatrix(Nchi,Nchi, 1)

    call convec(Nk,Nb,Nf,q,uk,ek,fk,lambda,phsus,ppsus)
    
    call mxCopyComplex16ToPtr(phsus,mxGetPr(plhs(1)), mxGetPi(plhs(1)),Nchi**2)
    call mxCopyComplex16ToPtr(ppsus,mxGetPr(plhs(2)), mxGetPi(plhs(2)),Nchi**2)

    DEALLOCATE(uk)
    DEALLOCATE(ek)
    DEALLOCATE(fk)
    DEALLOCATE(phsus)
    DEALLOCATE(ppsus)

    return
    end