    !
    !**** new module ************************************************************
    !
    module mTransporte
    !
    implicit none
    !
    real*8,  allocatable :: satElem(:), satElemAnt(:)
    real*8,  allocatable :: satElemL(:), satElemL0(:), satElem0(:)
    !
    public  :: transport, satcor, satinit, calc_prod
    private :: maxder, courant, ktdd_rt_STRANG, ktdd_rt_GODUNOV, ktdd_rt_CORREA, fkt, nbflux, derivau, udd_rt
    private :: xmax_speed, ymax_speed, zmax_speed, f, g, h, df, dg, dh
    private :: dif1o, mass
    contains
    !
    !=======================================================================
    !
    subroutine transport(v, GEOTIME)
    !
    !     Non-Linear Hyperbolic Equation
    !
    !     s,t + div(F(s)) = 0
    !
    !     F(s)= [f g]^t
    !
    !     Scheme: KT

    use mPropGeoFisica,    only: hx, hy, hz, phi, permkx
    use mMalha,            only: nsd, numelReserv, nen, numLadosElem, numLadosReserv
    use mMalha,            only: conecLadaisElem, local, listaDosElemsPorFace
    use mGlobaisEscalares, only: dtBlocoTransp,SPLITT
    use mHidrodinamicaRT,  only: ndofV
    use mGlobaisEscalares, only: tempoTotalTransporte
    use mPropGeoFisica,    only: tprt_prod,np_rand_prod, dtprt_prod
    !
    implicit none
    !
    real*8,  intent(in)    :: v(ndofV,numLadosReserv)
    !
    real(8), allocatable   :: uf(:)
    integer, dimension(100):: ntime
    !
    integer :: nt,nprt
    real(8) :: cr,crf,GEOTIME
    real(8) :: dt,vxmax,vymax,vzmax,vmax

    real*8 :: t1,t2
    !
    call timing(t1)
    !
    print*, "calculando o TRANSPORTE:"

    allocate(uf (numelReserv)); uf = 0.0d0
    ntime=0
    nprt =0
    !
    !      pi=4.d0*datan(1.d0)
    !
    vxmax=0.d00;vymax=0.d00;vzmax=0.d00;vmax=0.d00
    !
    !.... courant prescrito
    !
    cr = 0.125D0
    !
    !.... verificando a maior velocidade
    !
    call maxder(ndofV,numelReserv,vxmax,vymax,vzmax,vmax,phi,permkx,v)
    !
    if (nsd==2) write(*,*) "velocidades maximas",vxmax,vymax
    if (nsd==3) write(*,*) "velocidades maximas",vxmax,vymax,vzmax
    !
    !.... calculando numero de passos no tempo para o Courant fixo
    !
    if (nsd==2) nt = int(dtBlocoTransp*(vxmax/hx+vymax/hy)/cr)+1
    if (nsd==3) nt = int(dtBlocoTransp*(vxmax/hx+vymax/hy+vzmax/hz)/cr)+1
    !
    dt=dtBlocoTransp/nt
    print*, "Calculando para ", dtBlocoTransp, " dias"
    !
    write(*,*) 'dt=',dt,cr
    print*, "calculando saturacao com cr=", cr
    !
    !.... calculando o numero de Courant resultante
    !
    crf = courant(dt,hx,hy,hz,vxmax,vymax,vzmax)
    !
    write(*,*) 'Courant=',crf,' com ', nt,' passos para o transporte'
    !
    !.... Kurganov-Tadmor DxD
    !
    IF (SPLITT.EQ.'STRANG') THEN
        call ktdd_rt_STRANG(nsd,nt,nen, &
            conecLadaisElem,listaDosElemsPorFace,&
            dt,uf,v, numLadosElem,GEOTIME)
    ELSE
        IF(SPLITT.EQ.'CORREA')THEN
            call ktdd_rt_CORREA(nsd,nt, &
                conecLadaisElem,listaDosElemsPorFace,&
                dt,uf,v, numLadosElem,GEOTIME)
            write(*,*)"=================================================="
            write(*,*)"PERFORMING CORREA'S CORRECTION STEP FOR SATURATION"
            write(*,*)"=================================================="
        ELSE
            call ktdd_rt_GODUNOV(nsd,nt, &
                conecLadaisElem,listaDosElemsPorFace,&
                dt,uf,v, numLadosElem,GEOTIME)
        END IF
    ENDIF
    !
    tprt_prod=tprt_prod+dtprt_prod
    np_rand_prod = 0

    call timing(t2)
    tempoTotalTransporte=tempoTotalTransporte+(t2-t1)
    !
    end subroutine transport
    !
    !=======================================================================
    !
    subroutine maxder(ndof,numel,dxmax,dymax,dzmax,dmax,phi,perm,v)
    !
    use mMalha, only: conecLadaisElem, nsd, local
    !
    implicit none
    !
    integer :: numel,ndof
    real(8), dimension(*)    :: phi,perm
    real(8), dimension(ndof,*) :: v
    !
    real(8) :: srw,sro,ss,smin,smax,ds,dd,vv,phin,xk,yk,zk
    real(8) :: dmax,dxmax,dymax,dzmax,dxmin,dymin,dzmin
    real(8) :: uu
    integer :: nstep,n,nel
    !
    dmax =0.d0
    dxmax=0.d0
    dymax=0.d0
    if (nsd==3) dzmax=0.d0
    dxmin=1.d+10
    dymin=1.d+10
    if (nsd==3) dzmin=1.d+10
    !
    ss   = 0.0D0
    smin = srw
    smax = 1.d0-sro
    !
    nstep=100
    ds=(smax-smin)/nstep
    !
    do n=1,nstep+1
        ss=smin+(n-1)*ds
        !
        do nel=1,numel
            phin= phi(nel)
            !
            uu  = satElem(nel)
            xk  = perm(nel)
            yk  = perm(nel)
            if (nsd==3) zk = perm(nel)
            !
            !.... direcao x
            !
            vv = dabs(v(1,conecLadaisElem(1,nel)))
            !old     dd=dabs(dfx(ss,vv,xk))/phin
            dd = dabs(DF(ss,vv,xk))/phin
            if (dd.gt.dxmax) dxmax = dd
            !
            vv = dabs(v(1,conecLadaisElem(3,nel)))
            dd = dabs(DF(ss,vv,xk))/phin
            if (dd.gt.dxmax) dxmax = dd
            !
            !.... direcao y
            !
            vv = dabs(v(1,conecLadaisElem(2,nel)))
            !old      dd=dabs(dfy(ss,vv,xk))/phin
            dd = dabs(DG(ss,vv,xk))/phin
            if (dd.gt.dymax) dymax = dd
            !
            vv = dabs(v(1,conecLadaisElem(4,nel)))
            !old      dd=dabs(dfy(ss,vv,xk))/phin
            dd = dabs(DG(ss,vv,xk))/phin
            if (dd.gt.dymax) dymax = dd
            !
            !.... direcao z
            !
            if (nsd==3) then
                vv = dabs(v(1,conecLadaisElem(5,nel)))
                !old         dd = dabs(dfz(ss,vv,xk))/phin
                dd = dabs(DH(ss,vv,xk))/phin
                if (dd.gt.dzmax) dzmax = dd
                !
                vv = dabs(v(1,conecLadaisElem(6,nel)))
                !old         dd = dabs(dfz(ss,vv,xk))/phin
                dd = dabs(DH(ss,vv,xk))/phin
                if (dd.gt.dzmax) dzmax = dd
            endif
            !
        end do ! nel
        !
    end do ! n

    contains
    !
    !----------------------------------------------------------------------
    !
    function dfx(u,v,xk)
    implicit none
    !
    real(8) :: dfx,u,v,xk
    !
    dfx=v
    !
    end function
    !
    !----------------------------------------------------------------------
    !
    function dfy(u,v,xk)
    implicit none
    !
    real(8) :: dfy,u,v,xk
    !
    dfy=v
    !
    end function
    !
    !----------------------------------------------------------------------
    !
    function dfz(u,v,xk)
    implicit none

    real(8) :: dfz,u,v,xk
    !
    dfz=v
    !
    end  function
    !
    end subroutine
    !
    !=======================================================================
    !
    function courant(dt,hx,hy,hz,vxmax,vymax,vzmax)
    use mMalha, only: nsd
    !
    implicit none
    !
    real(8) :: dt,hx,hy,hz,courant,vxmax,vymax,vzmax
    !
    courant=dt*(vxmax/hx+vymax/hy)
    if(nsd==3) courant=courant + dt*(vzmax/hz)
    !
    end function
    !
    !=======================================================================
    !
    subroutine ktdd_rt_GODUNOV(nsd,nt,&
        conecLadaisElem,listaDosElemsPorFace,dt,uf,&
        v,numLadosElem, GEOTIME)
    !
    use mGlobaisEscalares, only: ns
    use mHidrodinamicaRT,  only: ndofV
    use mMalha,            only: numLadosReserv, numelReserv
    use mPropGeoFisica
    !
    implicit none
    !
    integer :: nsd,nt,numLadosElem
    integer :: conecLadaisElem(numLadosElem,numelReserv)
    integer :: listaDosElemsPorFace(numLadosElem, numLadosReserv)
    real(8) :: dt
    real(8), dimension(numelReserv)     :: uf
    real(8)   :: v(1,numLadosReserv)
    REAL*8  :: GEOTIME
    !
    integer :: nel
    real(8) :: umin,z,umax
    real(8), allocatable          :: du(:,:)
    real(8) :: t1
    !
    integer :: tid
    integer :: numThreads
    integer :: numPrimeiroElemento, numUltimoElemento
    integer :: inicioSol, fimSol
    integer, parameter :: meuNumMaxThreads=32
    real(8),dimension(meuNumMaxThreads)  :: uminAux,umaxAux
    !
    real(8) :: dttotal,dtt
    real(8), dimension(numelReserv)      :: phiINTERP
    !
    !
    !.... Mensagem: KT DxD
    !
    z=0.d0
    !
    if (allocated(du) .eqv. .false.) allocate(du(nsd,numelReserv))
    du = 0.0d0
    !
    !.... identifica os volumes e nohs internos e do contorno
    !
    uf=0.d0
    call nbflux(numelReserv,uf)
    !
    tid=1
    numThreads=1
    uminAux = 1.d+9
    umaxAux =-1.d+9
    umin = 1.d+9
    umax =-1.d+9
    !
    if(tid==1) print*, "numThreads",numThreads
    !
    !.... divisao do Trabalho para a atualizacao da solucao
    !
    numPrimeiroElemento = 1
    numUltimoElemento   = numelReserv
    call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, &
        numThreads, tid-1, inicioSol, fimSol)
    !
    call timing(t1)
    !
    do ns=1, nt
        !
        !.... Interpolacao da porosidade
        dttotal = real(ns)*dt
        dtt     = real(ns-1)*dt
        call interpPHI(phiINTERP,phi,phi0,dtt,GEOTIME,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....    calcula as derivadas
        !
        call derivau(nsd,numLadosElem,conecLadaisElem, &
            listaDosElemsPorFace,du,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....  Calcula os fluxos e evolui no tempo
        !
        call fkt(nsd,ndofV,dt,hx,hy,hz,uf,du,phiINTERP, &
            v,numLadosElem,permkx,permky,permkz,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....  Atualiza a solucao
        !
        do nel=inicioSol,fimSol
            satElem(nel)=satElemAnt(nel)
            if (satElem(nel).gt.umaxAux(tid)) umaxAux(tid) = satElem(nel)
            if (satElem(nel).lt.uminAux(tid)) uminAux(tid) = satElem(nel)
        end do
        !
        !$OMP BARRIER
        !
        !$OMP MASTER
        !
        !....  SATURATION UPDATE WITH POWER LAW FOR TIME FRACTION
        !
        CALL SATCOR(numelReserv,phi,phi0,satElem,dt,GEOTIME)
        !
        !      write(*,'(a,i0,2x,a,3f12.5)') "passo", ns,'massa global = ',sum(satElem), maxval(satElem), minval(satElem)
        !
        !$OMP END MASTER
        !
    end do ! ns
    !
    !     calcula a derivada para a reconstrucao linear
    !
    call derivau(nsd,numLadosElem,conecLadaisElem, &
        listaDosElemsPorFace,du,inicioSol,fimSol)
    !
    ! final da Paralelização
    !$OMP END PARALLEL
    !
    write(*,"(//,'Smax: ',f15.7,2x,'Smin: ',f15.7,//)") maxval(umaxAux), minval(uminAux)
    !
    end subroutine ktdd_rt_GODUNOV
    !
    !=======================================================================
    !
    subroutine ktdd_rt_STRANG(nsd,nt,nen,&
        conecLadaisElem,listaDosElemsPorFace,dt,uf,&
        v,numLadosElem, GEOTIME)
    !
    use mGlobaisEscalares, only: ns
    use mHidrodinamicaRT, only: ndofV
    use mMalha,            only: numLadosReserv, numelReserv
    use mPropGeoFisica
    !
    implicit none
    !
    integer :: nsd,nen,nt,numLadosElem
    integer, dimension(nen,*) :: conecLadaisElem
    integer, dimension(nen,*) :: listaDosElemsPorFace
    real(8) :: dt
    real(8), dimension(numelReserv)     :: uf
    real(8)   :: v(1,numLadosReserv), GEOTIME
    !
    integer :: nel
    real(8) :: umin,z,umax
    real(8), allocatable          :: du(:,:)
    real(8) :: t1
    !
    integer :: tid
    integer :: numThreads
    integer :: numPrimeiroElemento, numUltimoElemento
    integer :: inicioSol, fimSol
    integer, parameter :: meuNumMaxThreads=32
    real(8),dimension(meuNumMaxThreads)  :: uminAux,umaxAux
    !
    real(8) :: dttotal,dtt
    real(8), dimension(numelReserv)      :: phiINTERP
    !
    logical :: strang
    !.... STRANG FRACTIONAL TIME
    !
    REAL(8) :: FRACDT
    !
    !     Mensagem: KT DxD
    !
    z=0.d0
    !
    if (allocated(du) .eqv. .false.) allocate(du(nsd,numelReserv))
    du = 0.0d0
    !
    !.... identifica os volumes e nohs internos e do contorno
    !
    uf=0.d0
    call nbflux(numelReserv,uf)
    !
    tid=1
    numThreads=1
    uminAux = 1.d+9
    umaxAux =-1.d+9
    umin = 1.d+9
    umax =-1.d+9
    !
    strang=.false.
    !
    if(strang.eqv..false.) print*, "em transporte: GODUNOV"
    if(strang.eqv..true.)  print*, "em transporte: STRANG"
    !
    if(tid==1) print*, "numThreads",numThreads
    !
    !...  divisao do Trabalho para a atualizacao da solucao
    !
    numPrimeiroElemento = 1
    numUltimoElemento   = numelReserv
    call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, &
        numThreads, tid-1, inicioSol, fimSol)
    !
    call timing(t1)

    FRACDT = 0.50D0*dt
    !
    do ns=1, nt
        !
        !.... Interpolacao da porosidade
        dttotal = real(ns)*FRACDT
        dtt     = real(ns-1)*FRACDT
        call interpPHI(phiINTERP,phi,phi0,dtt,GEOTIME,inicioSol,fimSol)
        !
        !....   calcula as derivadas
        !
        call derivau(nsd,numLadosElem,conecLadaisElem,&
            listaDosElemsPorFace,du,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !     Calcula os fluxos e evolui no tempo
        !
        !! ONLY STRANG 27/Nov/2014
        call fkt(nsd,ndofV,FRACDT,hx,hy,hz,uf,du,phiINTERP,v, &
            numLadosElem,permkx,permky,permkz,inicioSol,fimSol)
        !$OMP SINGLE
        !
        !$OMP END SINGLE

        ! $OMP BARRIER
        !
        !....  Atualiza a solucao
        !
        do nel=inicioSol,fimSol
            satElem(nel)=satElemAnt(nel)
            if (satElem(nel).gt.umaxAux(tid)) umaxAux(tid) = satElem(nel)
            if (satElem(nel).lt.uminAux(tid)) uminAux(tid) = satElem(nel)
        end do

        !$OMP BARRIER
        !$OMP MASTER
        !
        !  SATURATION UPDATE WITH POWER LAW FOR TIME FRACTION
        !
        CALL SATCOR(numelReserv,phi,phi0,satElem,dt,GEOTIME)
        !$OMP END MASTER
        !$OMP BARRIER
        !.... Interpolacao da porosidade
        dttotal = real(ns+1)*FRACDT
        dtt     = real(ns)*FRACDT
        call interpPHI(phiINTERP,phi,phi0,dtt,GEOTIME,inicioSol,fimSol)
        !$OMP BARRIER
        call fkt(nsd,ndofV,FRACDT,hx,hy,hz,uf,du,phiINTERP,v,&
            numLadosElem,permkx,permky,permkz,inicioSol,fimSol)
        !$OMP BARRIER
        do nel=inicioSol,fimSol
            satElem(nel)=satElemAnt(nel)
            if (satElem(nel).gt.umaxAux(tid)) umaxAux(tid) = satElem(nel)
            if (satElem(nel).lt.uminAux(tid)) uminAux(tid) = satElem(nel)
        end do
        !
    end do ! ns
    !$OMP BARRIER
    !.... calcula a derivada para a reconstrucao linear
    !
    call derivau(nsd,numLadosElem,conecLadaisElem,&
        listaDosElemsPorFace,du,inicioSol,fimSol)
    !
    ! final da Paralelização
    !$OMP END PARALLEL
    !
    write(*,"(//,'Smax: ',f15.7,2x,'Smin: ',f15.7,//)") maxval(umaxAux), minval(uminAux)
    !
    end subroutine ktdd_rt_STRANG
    !
    !=======================================================================
    !
    subroutine fkt(nsd,ndof,dt,hx,hy,hz, &
        &     uf,du,phi,v,numLadosElem,permx,permy,permz,inicio,fim)
    !
    implicit none
    !
    integer, intent(in) :: nsd,numLadosElem,ndof
    real(8), intent(in) :: dt,hx,hy,hz
    real(8), intent(in)       :: v(ndof,*)
    real(8), dimension(*)     :: uf
    real(8), dimension(*), intent(in) :: permx,permy,permz,phi
    real(8), dimension(nsd,*),intent(in) :: du
    !
    integer :: k,n(numLadosElem),nel,posDer
    real(8) :: fluxo(2,numLadosElem)
    real(8),dimension(numLadosElem) :: uel,vel,unb,xkb, p
    real(8) :: pc,xk,yk,zk,uc

    real(8) :: eps=1.d-8
    integer :: inicio, fim

    xk=0.0; yk=0.0; zk=0.0
    uel=0.0; vel=0.0; unb=0.0; p=0.0
    pc=0.0d0
    fluxo=0.0d00

    do nel=inicio,fim
        !
        call vizinhanca(uel,vel,unb,xkb,p,pc,numLadosElem)

        if (n(1)==nel) then
            if(vel(1).lt.eps) then
                fluxo(1,1)=f(unb(1),vel(1),xkb(1))
                fluxo(2,1)=f(uel(1),vel(1),xk)
            else
                fluxo(1,1)=f(uf(nel),vel(1),xk)
                fluxo(2,1)=f(uf(nel),vel(1),xk)
            end if
        else
            fluxo(1,1)=f(unb(1),vel(1),xkb(1))
            fluxo(2,1)=f(uel(1),vel(1),xk)
        end if
        !
        if (n(2)==nel) then
            if(vel(2).lt.eps) then
                fluxo(1,2)=g(unb(2),vel(2),xkb(2))
                fluxo(2,2)=g(uel(2),vel(2),yk)
            else
                fluxo(1,2)=g(uf(nel),vel(2),yk)
                fluxo(2,2)=g(uf(nel),vel(2),yk)
            end if
        else
            fluxo(1,2)=g(unb(2),vel(2),xkb(2))
            fluxo(2,2)=g(uel(2),vel(2),yk)
        endif
        !
        if (n(3)==nel) then
            if(vel(3).gt.eps) then
                fluxo(1,3)=f(unb(3),vel(3),xkb(3))
                fluxo(2,3)=f(uel(3),vel(3),xk)
            else
                fluxo(1,3)=f(uf(nel),vel(3),xk)
                fluxo(2,3)=f(uf(nel),vel(3),xk)
            end if
        else
            fluxo(1,3)=f(uel(3),vel(3),xk)
            fluxo(2,3)=f(unb(3),vel(3),xkb(3))
        end if

        if (n(4)==nel) then
            if(vel(4).gt.eps) then
                fluxo(1,4)=g(unb(4),vel(4),xkb(4))
                fluxo(2,4)=g(uel(4),vel(4),yk)
            else
                fluxo(1,4)=g(uf(nel),vel(4),yk)
                fluxo(2,4)=g(uf(nel),vel(4),yk)
            end if
        else
            fluxo(1,4)=g(uel(4),vel(4),yk)
            fluxo(2,4)=g(unb(4),vel(4),xkb(4))
        endif
        !
        if(nsd==3) then
            if (n(5)==nel) then
                if(vel(5).lt.eps) then
                    fluxo(1,5)=h(unb(5),vel(5),xkb(5))
                    fluxo(2,5)=h(uel(5),vel(5),zk)
                else
                    fluxo(1,5)=h(uf(nel),vel(5),zk)
                    fluxo(2,5)=h(uf(nel),vel(5),zk)
                end if
            else
                fluxo(1,5)=h(unb(5),vel(5),xkb(5))
                fluxo(2,5)=h(uel(5),vel(5),zk)
            end if
            !
            if (n(6)==nel) then
                if(vel(6).gt.eps) then
                    fluxo(1,6)=h(unb(6),vel(6),xkb(6))
                    fluxo(2,6)=h(uel(6),vel(6),zk)
                else
                    fluxo(1,6)=h(uf(nel),vel(6),zk)
                    fluxo(2,6)=h(uf(nel),vel(6),zk)
                end if
            else
                fluxo(1,6)=h(uel(6),vel(6),zk)
                fluxo(2,6)=h(unb(6),vel(6),xkb(6))
            end if
        endif

        satElemAnt(nel)=udd_rt(dt,hx,hy,hz,uc,uel,unb,vel,fluxo,p,pc,xk,yk,zk,xkb)
        !
    end do
    !
    contains

    subroutine vizinhanca(uel,vel,unb,xkb,p,pc,numLadosElem)

    use mMalha,            only: local, numelReserv, conecNodaisElem
    use mMalha,            only: conecLadaisElem, listaDosElemsPorFace

    implicit none

    integer, intent(in) :: numLadosElem
    real(8), intent(inout) :: pc
    real(8), dimension(numLadosElem) :: uel,vel,unb,xkb, p
    integer :: i


    n(1)= listaDosElemsPorFace(2,conecLadaisElem(4,nel))
    n(2)= listaDosElemsPorFace(3,conecLadaisElem(1,nel))
    n(3)= listaDosElemsPorFace(4,conecLadaisElem(2,nel))
    n(4)= listaDosElemsPorFace(1,conecLadaisElem(3,nel))
    if(nsd==3)n(5)= listaDosElemsPorFace(5,conecLadaisElem(5,nel))
    if(nsd==3)n(6)= listaDosElemsPorFace(6,conecLadaisElem(6,nel))

    do i=1, numLadosElem
        if(n(i)>numelReserv) n(i)=0
    enddo

    if(n(1)==0) n(1)=nel
    if(n(2)==0) n(2)=nel
    if(n(3)==0) n(3)=nel
    if(n(4)==0) n(4)=nel
    if(nsd==3) then
        if(n(5)==0) n(5)=nel !z
        if(n(6)==0) n(6)=nel !z
    endif
    !
    !     porosidade nos elementos vizinhos
    !
    pc  =phi(nel)
    p(1)=phi(n(1))
    p(2)=phi(n(2))
    p(3)=phi(n(3))
    p(4)=phi(n(4))
    if(nsd==3) then
        p(5)=phi(n(5)) !z
        p(6)=phi(n(6)) !z
    endif
    !
    !     permeabilidades
    !
    xk    =permx(nel)
    yk    =permy(nel)
    if(nsd==3) zk    =permz(nel)

    xkb(1)=permx(n(1))
    xkb(2)=permy(n(2))
    xkb(3)=permx(n(3))
    xkb(4)=permy(n(4))
    if(nsd==3) then
        xkb(5)=permz(n(5)) !z
        xkb(6)=permz(n(6)) !z
    endif
    !
    !     solucao no elemento
    !
    uc=satElem(nel)
    !
    !     solucoes nos lados: reconstrucao linear
    !
    uel(1)=uc-du(1,nel)/2.d0
    uel(2)=uc-du(2,nel)/2.d0
    uel(3)=uc+du(1,nel)/2.d0
    uel(4)=uc+du(2,nel)/2.d0
    if(nsd==3) then
        uel(5)=uc-du(3,nel)/2.d0 !z
        uel(6)=uc+du(3,nel)/2.d0 !z
    endif
    !
    !     velocidades nos lados do elemento
    !
    vel(1)=v(1,conecLadaisElem(4,nel))
    vel(2)=v(1,conecLadaisElem(1,nel))
    vel(3)=v(1,conecLadaisElem(2,nel))
    vel(4)=v(1,conecLadaisElem(3,nel))
    if(nsd==3) then
        vel(5)=v(1,conecLadaisElem(5,nel)) !z
        vel(6)=v(1,conecLadaisElem(6,nel)) !z
    endif

    !     solucoes nos elementos vizinhos: reconstrucao linear
    !
    do k=1, numLadosELem

        if(n(k)==0)then
            unb(k)=uel(k)
        else
            if(k==1.or.k==3)posDer=1
            if(k==2.or.k==4)posDer=2
            if(k==5.or.k==6)posDer=3

            if(k==1.or.k==2.or.k==5) then
                unb(k)=satElem(n(k))+du(posDer,n(k))/2.d0
            else
                unb(k)=satElem(n(k))-du(posDer,n(k))/2.d0
            end if

        end if
    end do

    end subroutine

    END SUBROUTINE
    !
    !=======================================================================
    !
    subroutine nbflux(numel, uf)
    !
    !     identifica os volumes do contorno com fluxo prescrito
    !
    use mMalha,         only: xc, nsd
    use mPropGeoFisica, only: hx, hy, hz, dimx, dimy, dimz, sinj
    !
    implicit none
    !
    integer :: numel
    real(8), dimension(*)   :: uf
    !
    integer :: nel
    real(8) :: xx,yy,zz
    !
    !       dimx=dabs(x(1,numnp)-x(1,1))
    !       dimy=dabs(x(2,numnp)-x(2,1))
    !       if(nsd==3)dimz=dabs(x(3,numnp)-x(3,1))
    !       hx  =dimx/nelx
    !       hy  =dimy/nely
    !       if(nsd==3)hz  =dimz/nelz

    do nel=1,numel
        !
        uf(nel) =0.d0
        xx=xc(1,nel)
        yy=xc(2,nel)
        if(nsd==3)zz=xc(3,nel)
        !
        !     injecao continua no bordo esquerdo
        !
        if(xx.lt.hx) uf(nel)=sinj
        !
        !     injecao continua no bordo direito
        !
        if(xx.gt.(dimx-hx)) uf(nel)=sinj
        !
        !     injecao continua no bordo frente
        !
        if(yy.lt.hy) uf(nel)=sinj
        !
        !     injecao continua no bordo atras
        !
        if(yy.gt.(dimy-hy)) uf(nel)=sinj

        if(nsd==3) then
            !
            !        injecao continua no bordo inferior
            !
            if(zz.lt.hz) uf(nel)=sinj
            !
            !        injecao continua no bordo superior
            !
            if(zz.gt.(dimz-hz)) uf(nel)=sinj
            !
        endif
        !
        !
    end do
    !
    end subroutine
    !
    !=======================================================================
    !
    subroutine derivau(nsd,numLadosElem,conecLadaisElem,listaDosElemsPorFace,du,inicio,fim)
    !
    !     calcula as derivadas na malha original
    !
    use mMalha,         only: numelReserv, numLadosReserv
    !
    implicit none
    !
    integer :: nsd,numLadosElem
    integer :: conecLadaisElem(numLadosElem,numelReserv)
    integer :: listaDosElemsPorFace(numLadosElem, numLadosReserv)
    real(8) :: du(nsd, numelReserv)
    !
    integer :: nel,n(numLadosElem)
    real(8) :: u1,u2,u3,u4,u5,u6,ux,uy,uz
    integer  :: inicio,fim
    !
    !     volumes interiores
    !
    do nel=inicio,fim
        !
        !     derivada em x
        !
        n(1)= listaDosElemsPorFace(2,conecLadaisElem(4,nel))
        n(2)= listaDosElemsPorFace(3,conecLadaisElem(1,nel))
        n(3)= listaDosElemsPorFace(4,conecLadaisElem(2,nel))
        n(4)= listaDosElemsPorFace(1,conecLadaisElem(3,nel))
        if(nsd==3)n(5)= listaDosElemsPorFace(5,conecLadaisElem(5,nel))
        if(nsd==3)n(6)= listaDosElemsPorFace(6,conecLadaisElem(6,nel))
        !
        !           print*, "vizinhosLadal", nel, n(1:4)
        !
        if(n(1)==0) then
            u1=satElem(nel)
        else
            u1=satElem(n(1))
        endif
        !
        if(n(3)==0) then
            u3=satElem(nel)
        else
            u3=satElem(n(3))
        endif
        !
        ux=diffu(u1,satElem(nel),u3)
        du(1,nel)=ux
        !
        !     derivada em y
        !
        if(n(2)==0) then
            u2=satElem(nel)
        else
            u2=satElem(n(2))
        endif

        if(n(4)==0) then
            u4=satElem(nel)
        else
            u4=satElem(n(4))
        endif

        uy=diffu(u2,satElem(nel),u4)
        du(2,nel)=uy
        !
        !     derivada em z
        !
        if(nsd==3) then

            if(n(5)==0) then
                u5=satElem(nel)
            else
                u5=satElem(n(5))
            endif

            if(n(6)==0) then
                u6=satElem(nel)
            else
                u6=satElem(n(6))
            endif
            !
            uz=diffu(u5,satElem(nel),u6)
            du(3,nel)=uz

        endif
        !
    end do


    CONTAINS

    !
    ! =======================================================================
    !
    function diffu(xi,xc,xf)
    implicit none
    !
    real(8) :: diffu
    real(8) :: xi,xc,xf
    real(8) :: alfa,d1,d2,d3
    !
    !       alfa=0.0d0
    !
    !      se alfa=0, modificar a funcao dif1o
    !
    !      alfa=0.5d0
    !       alfa=1.0d0
    alfa=2.0d0
    !       alfa=1.95d0
    !      alfa=2.1d0
    !
    !      calculo das derivadas
    !
    d1= (xc-xi)*alfa
    d2= (xf-xi)*0.5d0
    d3= (xf-xc)*alfa
    !
    diffu=mm(d1,d2,d3)
    !
    end function
    !
    ! =======================================================================
    !
    function mm(a,b,c)
    !
    implicit none
    !
    real(8) :: mm
    real(8) :: a,b,c
    real(8) :: mm1
    !
    mm1=minmod2(a  ,b)
    mm =minmod2(mm1,c)
    !
    return
    !
    end function
    !
    ! =======================================================================
    !
    function minmod2(a,b)
    !
    !      If a and b have the same sign, then minmod will return whichever
    !      has smaller modulus. Otherwise, it will return 0.
    !
    implicit none
    !
    real(8) :: a,b,minmod2,z
    !
    z=0.d0
    !
    if (a*b.gt.z) then
        if (dabs(a).lt.dabs(b)) then
            minmod2 = a
        else
            minmod2 = b
        end if
    else
        minmod2 = z
    end if
    !
    return
    end function
    !
    end subroutine derivau
    !
    !=======================================================================
    !
    function udd_rt(dt,hx,hy,hz,uc,uel,unb,vel,flux,p,pc,xk, &
        &      yk,zk,xkb)
    !
    !     KT DxD na fronteira do dominio
    !
    use mMalha, only: numLadosElem, nsd
    !
    implicit none
    !
    real(8) :: dt,hx,hy,hz,uc,udd_rt
    real(8), dimension(2,numLadosElem) :: flux
    real(8), dimension(numLadosElem)   :: uel,unb,vel,p,xkb
    real(8) :: a1,a2,a3,a4,a5,a6,ff,pc
    real(8) :: pm1,pm2,pm3,pm4,pm5,pm6,xk,yk,zk
    !
    !     calcula as derivadas dos fluxos para a construcao da malha auxiliar
    !
    a1=xmax_speed(uel(1),unb(1),vel(1),xk,xkb(1),pc,p(1))
    a2=ymax_speed(uel(2),unb(2),vel(2),yk,xkb(2),pc,p(2))
    a3=xmax_speed(uel(3),unb(3),vel(3),xk,xkb(3),pc,p(3))
    a4=ymax_speed(uel(4),unb(4),vel(4),yk,xkb(4),pc,p(4))
    if(nsd==3) then
        a5=zmax_speed(uel(5),unb(5),vel(5),zk,xkb(5),pc,p(5))
        a6=zmax_speed(uel(6),unb(6),vel(6),zk,xkb(6),pc,p(6))
    endif
    !
    !      if(nnp.eq.1 .and. ns.eq.1) then
    !      eps=1.d-3
    !      if(dabs(vel(1)).gt.eps) a1=hx/dt/8.d0
    !      if(dabs(vel(2)).gt.eps) a2=hy/dt/8.d0
    !      if(dabs(vel(3)).gt.eps) a3=hx/dt/8.d0
    !      if(dabs(vel(4)).gt.eps) a4=hy/dt/8.d0
    !      if(dabs(vel(5)).gt.eps) a5=hz/dt/8.d0
    !      if(dabs(vel(6)).gt.eps) a6=hz/dt/8.d0
    !      end if

    !
    !...  Porosidades "medias"
    !
    pm1 = pc+p(1)
    pm2 = pc+p(2)
    pm3 = pc+p(3)
    pm4 = pc+p(4)
    if(nsd==3) then
        pm5 = pc+p(5)
        pm6 = pc+p(6)
    endif
    !
    !     Fluxo total
    !
    ff=((a3*p(3)*(unb(3)-uel(3)) -(flux(2,3)-flux(1,3)))/pm3   &
        &   +(a1*p(1)*(unb(1)-uel(1)) -(flux(2,1)-flux(1,1)))/pm1   &
        &                             -(flux(1,3)-flux(2,1))/pc)/hx &
        &   +((a4*p(4)*(unb(4)-uel(4))-(flux(2,4)-flux(1,4)))/pm4   &
        &   +(a2*p(2)*(unb(2)-uel(2)) -(flux(2,2)-flux(1,2)))/pm2   &
        &                             -(flux(1,4)-flux(2,2))/pc)/hy

    if(nsd==3) ff = ff +                                     &
        &   ((a6*p(6)*(unb(6)-uel(6)) -(flux(2,6)-flux(1,6)))/pm6   &
        &   +(a5*p(5)*(unb(5)-uel(5)) -(flux(2,5)-flux(1,5)))/pm5   &
        &                             -(flux(1,6)-flux(2,5))/pc)/hz
    !
    !     Evolucao no tempo (Runge-Kutta)
    !
    udd_rt=rk(uc,dt,ff)

    !
    CONTAINS
    !
    !=======================================================================
    !
    function rk(uc,dt,ff)
    !
    implicit none
    !
    real(8) :: rk,uc,dt,ff
    real(8) :: r1,r2,r3
    !
    !     Evolucao no tempo (euler)
    !
    rk=uc + dt*ff
    !
    end function rk
    !
    !=======================================================================
    !
    function rkT1(uc,dt,ff,ns,hx,hy,hz,uel,unb,vel,flux,p,pc,xk,yk,zk,xkb,nel)
    !
    implicit none
    !
    real(8) :: rkT1,uc,dt,ff
    integer :: ns,nel
    real*8  :: hx,hy,hz,pc,xk,yk,zk
    real(8), dimension(2,6) :: flux
    real(8), dimension(6)   :: uel,unb,vel,p,xkb
    !
    rkT1 = uc+dt*ff

    end function
    !
    !=======================================================================
    !
    function rkT2(uc,dt,ff,ns,hx,hy,hz,uel,unb,vel,flux,p,pc,xk,yk,zk,xkb,nel)
    !
    use mGlobaisArranjos, only : uTempoN
    implicit none
    !
    real(8) :: rkT2,uc,dt,ff
    integer :: ns,nel
    real*8  :: hx,hy,hz,pc,xk,yk,zk
    real(8), dimension(2,*) :: flux
    real(8), dimension(*)   :: uel,unb,vel,p,xkb
    !
    !     Evolucao no tempo (Runge-Kutta)
    !
    !     terceira ordem
    !
    rkT2 = uTempoN(nel)*3.d0/4.d0 + (uc+dt*ff)/4.d0

    end function
    !
    !=======================================================================
    !
    function rkT3(uc,dt,ff,ns,hx,hy,hz,uel,unb,vel,flux,p,pc,xk,yk,zk,xkb,nel)
    !
    use mGlobaisArranjos, only : uTempoN
    !
    implicit none
    !
    real(8) :: rkT3,uc,dt,ff
    integer :: ns,nel
    real*8  :: hx,hy,hz,pc,xk,yk,zk
    real(8), dimension(2,*) :: flux
    real(8), dimension(*)   :: uel,unb,vel,p,xkb

    rkT3 = uTempoN(nel)/3.d0 + (uc+dt*ff)*2.d0/3.d0

    end function

    !
    !=======================================================================
    !
    function rkT2_ordem2(uc,dt,ff,ns,hx,hy,hz,uel,unb,vel,flux,p,pc,xk,yk,zk,xkb,nel)
    !
    use mGlobaisArranjos, only : uTempoN
    implicit none
    !
    real(8) :: rkT2_ordem2,uc,dt,ff
    integer :: ns,nel
    real*8  :: hx,hy,hz,pc,xk,yk,zk
    real(8), dimension(2,*) :: flux
    real(8), dimension(*)   :: uel,unb,vel,p,xkb
    !
    !     Evolucao no tempo (Runge-Kutta)
    !
    !     terceira ordem
    !
    rkT2_ordem2 = uTempoN(nel)*1.d0/2.d0 + (uc+dt*ff)/2.d0
    !
    end function
    !
    end function udd_rt
    !
    !=======================================================================
    !
    subroutine satinit(nsd,numel,u)
    !
    use mMalha, only: xc
    use mPropGeoFisica
    !
    implicit none
    !
    integer :: nsd,numel
    real*8  :: u(numel)
    !
    integer :: i
    real*8  :: xi,xf,yi,yf,zi,zf
    !
    !     condicao inicial: agua e oleo
    !
    !.... Bloco
    !
    xi=xcbloco-xlbloco/2.d0
    xf=xcbloco+xlbloco/2.d0
    yi=ycbloco-ylbloco/2.d0
    yf=ycbloco+ylbloco/2.d0
    if(nsd==3) zi=zcbloco-zlbloco/2.d0
    if(nsd==3) zf=zcbloco+zlbloco/2.d0

    do i=1,numel

        u(i)=sinicial

        if (nsd==2) then
            if(xc(1,i).gt.xi.and.xc(1,i).lt.xf) then
                if(xc(2,i).gt.yi.and.xc(2,i).lt.yf) then
                    u(i)=sbloco
                end if
            end if
        else
            if(xc(1,i).gt.xi.and.xc(1,i).lt.xf) then
                if(xc(2,i).gt.yi.and.xc(2,i).lt.yf) then
                    if(xc(3,i).gt.zi.and.xc(3,i).lt.zf) then
                        u(i)=sbloco
                    end if
                end if
            end if
        end if
        !
    end do
    !
    end subroutine

    !
    !=======================================================================
    !
    function xmax_speed(u1,u2,vv,xk1,xk2,phi1,phi2)
    !
    implicit none
    !
    real(8) :: xmax_speed
    real(8) :: u1,u2,vv,xk1,xk2,phi1,phi2
    !
    real(8) :: du,a1,a2,aa,a,uu
    integer :: nu,n
    nu=2
    du=(u2-u1)/nu
    !
    a=0.d0
    !
    do n=1,nu+1
        uu=u1+(n-1)*du
        a1=df(uu,vv,xk1)/phi1
        a2=df(uu,vv,xk2)/phi2
        aa=max(dabs(a1),dabs(a2))
        if(aa.gt.a) a = aa
    end do
    !
    xmax_speed=a
    !
    end function
    !
    !=======================================================================
    !
    function ymax_speed(u1,u2,vv,xk1,xk2,phi1,phi2)
    !
    implicit none
    !
    real(8) :: ymax_speed
    real(8) :: u1,u2,vv,xk1,xk2,phi1,phi2
    !
    real(8) :: du,a1,a2,aa,a,uu
    integer :: nu,n
    nu=2
    du=(u2-u1)/nu
    !
    a=0.d0
    !
    do n=1,nu+1
        uu=u1+(n-1)*du
        a1=dg(uu,vv,xk1)/phi1
        a2=dg(uu,vv,xk2)/phi2
        aa=max(dabs(a1),dabs(a2))
        if(aa.gt.a) a = aa
    end do
    !
    ymax_speed=a
    !
    end function
    !
    !=======================================================================
    !
    function zmax_speed(u1,u2,vv,xk1,xk2,phi1,phi2)
    !
    implicit none
    !
    real(8) :: zmax_speed
    real(8) :: u1,u2,vv,xk1,xk2,phi1,phi2
    !
    real(8) :: du,a1,a2,aa,a,uu
    integer :: nu,n
    nu=2
    du=(u2-u1)/nu
    !
    a=0.d0
    !
    do n=1,nu+1
        uu=u1+(n-1)*du
        a1=dh(uu,vv,xk1)/phi1
        a2=dh(uu,vv,xk2)/phi2
        aa=max(dabs(a1),dabs(a2))
        if(aa.gt.a) a = aa
    end do
    !
    zmax_speed=a
    !
    end function
    !
    !=======================================================================
    !
    subroutine calc_prod(nedfl,nedg,numel,nsd,nen,ien,ieedg,vel,sat,xx,t0)
    !
    !     calcula a velocidade nos nos por elemento
    !
    use mPropGeoFisica, only: nelxReserv, nelyReserv, dimx, dimy, hy
    use mPropGeoFisica, only: xlo, xlt, nr, prod_out, np_rand_prod,ncont_prod,npprod
    use mMalha,         only: numelReserv

    implicit none
    !
    integer :: nedfl,nedg,numel,nsd,nen,infile
    integer :: np1,np2,np,l1
    integer, dimension(nen,*)     :: ien
    real(8), dimension(numel)     :: sat
    real(8), dimension(nedfl,*)   :: vel
    real(8), dimension(nsd,*)     :: xx
    real(8) :: t0,prod
    real(8) :: st,mobi,aux1,aux2,tol
    integer, dimension(nedg,*)    :: ieedg
    !
    integer :: nel
    !
    character(len=128) :: NAME,FILE_IN
    character(len=4)   :: EXT
    character(len=5)   :: C
    !
    prod = 0.0d0
    !
    !     definir tipo de numeraccao
    !
    np1 = ien(2,nelxReserv)
    np2 = ien(4,nelyReserv)
    aux1 = abs(xx(1,np1)-dimx)
    aux2 = abs(xx(2,np2)-dimy)
    tol = 1e-09
    !
    if(aux1.lt.tol) np = 1
    if(aux2.lt.tol) np = 2
    !
    !.....  loop on elements
    !
    if(np.eq.1)then
        !         do nel=nelxReserv,numel,nelxReserv
        do nel=nelxReserv,numelReserv,nelxReserv
            !
            l1   = ieedg(2,nel)
            st   = sat(nel)
            mobi = xlo(st)/xlt(st)
            prod = prod + vel(1,l1)*mobi
            !
        end do
        !
    else
        !
        do nel=numelReserv-(nelyReserv),numelReserv
            !
            l1   = ieedg(2,nel)
            st   = sat(nel)
            mobi = xlo(st)/xlt(st)
            prod = prod + vel(1,l1)*mobi
            !
        end do
        !
    end if
    !
    infile=68
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    FILE_IN=prod_out
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    !      WRITE(*,111)nr,NAME(1:LEN_TRIM(NAME))
    open(infile,FILE=NAME)
    !
    if(np_rand_prod.eq.1)then
        !
        ncont_prod = 0
        !
    end if
    !
    !       write(infile,5566)t0,prod*hy
    write(infile,5566)t0,prod*hy
    !
    ncont_prod = ncont_prod+1
    !
    if(npprod.eq.ncont_prod) close(infile)
    !
111 FORMAT('NAME OF PRODUCTION OUTPUT FILE',I5,': ',A)
5566 format(3(f15.8,2x))
113 FORMAT(I5)
    !
    return
    end subroutine
    !
    ! =======================================================================
    !
    subroutine calc_prod3D(nen,numel,nsd,conecNodaisElem,vel,sat,xx,tTransporte)
    !
    !      calcula a velocidade nos nos por elemento
    !
    use mPropGeoFisica
    !
    implicit none
    !
    integer :: nen,numel,nsd,infile
    integer :: np1,np2,np3,np
    integer, dimension(nen,*)     :: conecNodaisElem
    real(8), dimension(numel)     :: sat
    real(8), dimension(6,numel)    :: vel
    real(8), dimension(nsd,*)     :: xx
    real(8) :: tTransporte,prod,vec
    real(8) :: st,mobi,aux1,aux2,aux3,tol
    !
    integer :: nel
    !
    character(len=128) :: NAME,FILE_IN
    character(len=4)   :: EXT
    character(len=5)   :: C
    !
    prod = 0.0d0
    vec  = 0.0d0
    !
    !      definir tipo de numeraccao
    !
    np1 = conecNodaisElem(2,nelx)
    np2 = conecNodaisElem(4,nely)
    np3 = conecNodaisElem(5,nelz)
    aux1 = abs(xx(1,np1)-dimx)
    aux2 = abs(xx(2,np2)-dimy)
    aux3 = abs(xx(3,np3)-dimz)
    tol = 1e-09
    !
    if(aux1.lt.tol) np = 1
    if(aux2.lt.tol) np = 2
    if(aux3.lt.tol) np = 3 !!??
    !
    ! .....  loop on elements
    !
    np=2
    if(np.eq.1)then
        do nel=numel-(nelx*nely)+1,numel
            !
            st   = sat(nel)
            mobi = xlo(st)/xlt(st)
            prod = prod + vel(6,nel)*mobi
            vec  = vec + vel(6,nel)
            !
        end do
        !
    else
        !
        nel=numel
        st   = sat(nel)
        mobi = xlo(st)/xlt(st)
        prod = prod + vel(3,nel)*mobi + vel(4,nel)*mobi + vel(6,nel)*mobi
        vec  = vec + vel(3,nel) + vel(4,nel) + vel(6,nel)
        !
    end if
    !
    infile=68
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    FILE_IN=prod_out
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    !       WRITE(*,111)nr,NAME(1:LEN_TRIM(NAME))
    open(infile,FILE=NAME)
    !
    if(np_rand_prod.eq.1)then
        !
        ncont_prod = 0
        !
    end if
    !
    write(infile,5566)tTransporte,(prod*hy)/(vec*hy)
    !
    ncont_prod = ncont_prod+1
    !
    if(npprod.eq.ncont_prod) close(infile)
    !
    !       111  FORMAT('NAME OF PRODUCTION OUTPUT FILE',I5,': ',A)
5566 format(2(f15.8,2x))
113 FORMAT(I5)
    !
    return
    end subroutine
    !
    ! =======================================================================
    !
    function f(uu,vv,xk)
    !
    !       calcula f em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlo,xlt,gf1,rhoo,rhow

    !
    implicit none
    !
    real(8) :: f,uu,vv,xk
    !
    real(8) :: fw
    real(8) :: xlwu,xlou
    real(8) :: xemp,xg
    !
    xlwu=xlw(uu)
    xlou=xlo(uu)
    fw=xlwu/(xlwu+xlou)
    xg = gf1
    xemp=fw*xg*xlou*xk*(rhoo-rhow)
    !
    f=fw*vv-xemp
    !
    end function
    !
    ! =======================================================================
    !
    function g(uu,vv,xk)
    !
    !       calcula f em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlo,xlt,gf2,rhoo,rhow
    !
    implicit none
    real(8) :: g,uu,vv,xk
    !
    real(8) :: fw
    real(8) :: xlwu,xlou
    real(8) :: xemp,xg
    !
    xlwu=xlw(uu)
    xlou=xlo(uu)
    fw=xlwu/(xlwu+xlou)
    xg = gf2
    xemp=fw*xg*xlou*xk*(rhoo-rhow)
    g=fw*vv-xemp

    end function
    !
    ! =======================================================================
    !
    function h(uu,vv,xk)
    !
    !       calcula f em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlo,xlt,gf3,rhoo,rhow
    !
    implicit none
    real(8) :: h,uu,vv,xk
    !
    real(8) :: fw
    real(8) :: xlwu,xlou
    real(8) :: xemp,xg
    !
    xlwu=xlw(uu)
    xlou=xlo(uu)
    fw=xlwu/(xlwu+xlou)
    
    xg = gf3
    xemp=fw*xg*xlou*xk*(rhoo-rhow)
    !
    h=fw*vv-xemp
    !
    end function
    !
    ! =======================================================================
    !
    function df(uu,vv,xk)
    !
    !       calcula df/du em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlt,xlo,gf1,iflag_linear,srw,sro,xmiw,xmio,rhoo,rhow
    !
    implicit none
    real(8) :: df,uu,vv
    real(8) :: dfw,dlt,dlo,dlw,dfe
    real(8) :: uus,zero
    real(8) :: xg,xk
    !
    zero=0.d0
    !
    !  .... Buckley-Leverett
    !
    !       f=(uu**2)/(uu**2+a*(1-uu)**2)
    !        a=1.d0/2.d0
    !        a1=uu**2+a*(1.d0-uu)**2
    !        df=(2.d0*uu/a1 - uu**2*(2.d0*uu-2.d0*a*(1.d0-uu))/a1**2)*vv
    !
    !  .... Burgers, p. 464 do LeVeque
    !
    !       a=dsqrt(2.d0)/2.d0 ! 45 graus
    !
    !       f=a*uu**2/2.d0
    !
    !       df=a*uu
    !
    !..... Agua-Oleo
    !
    xg = gf1
    !
    select case(iflag_linear)
    case(1)
        !
        !..... derivada do fluxo
        !
        df=vv
        !
    case(2)
        !
        !..... derivadas das mobilidades
        !
        uus = uu-srw
        if (uus.lt.zero) uus = zero
        dlw = 2.d0*(uus)/xmiw/(1.d0-srw)**2
        !
        uus = 1.d0-sro-uu
        if (uus.lt.zero) uus = zero
        dlo =-2.d0*(uus)/xmio/(1.d0-sro)**2
        !
        dlt = dlw+dlo
        !
        !..... derivada do fluxo fracionario
        !
        !..... derivada do fluxo
        !
        dfw = dlw/xlt(uu) - xlw(uu)*dlt/xlt(uu)**2
        !
        dfe = (dfw*xlo(uu)+xlw(uu)/xlt(uu)*dlo)
        !
        df = dfw*vv-(rhoo-rhow)*xg*xk*dfe
        !
    end select
    !
    end function
    !
    !=======================================================================
    !
    function dg(uu,vv,xk)
    !
    !..... calcula dg/du em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlt,xlo,gf2,iflag_linear,srw,sro,xmiw,xmio,rhoo,rhow
    !
    !
    implicit none
    real(8) :: dg,uu,vv
    real(8) :: dfw,dlt,dlo,dlw,dfe
    real(8) :: uus,zero
    real(8) :: xg,xk
    !
    zero=0.d0
    !
    !  .... Buckley-Leverett
    !
    !      g=(uu**2)/(uu**2+a*(1-uu)**2)
    !        a=1.d0/2.d0
    !        a1=uu**2+a*(1.d0-uu)**2
    !        dg=(2.d0*uu/a1 - uu**2*(2.d0*uu-2.d0*a*(1.d0-uu))/a1**2)*vv
    !
    !  .... Burgers, p. 464 do LeVeque
    !
    !       a=dsqrt(2.d0)/2.d0 ! 45 graus
    !
    !       g=a*uu**2/2.d0
    !
    !       dg=a*uu
    !
    !  .... Agua-Oleo
    !
    xg=gf2
    !
    select case(iflag_linear)
    case(1)
        !
        !..... derivada do fluxo
        !
        dg=vv
        !
    case(2)
        !
        !..... derivadas das mobilidades
        !
        uus = uu-srw
        if (uus.lt.zero) uus = zero
        dlw = 2.d0*(uus)/xmiw/(1.d0-srw)**2
        !
        uus = 1.d0-sro-uu
        if (uus.lt.zero) uus = zero
        dlo =-2.d0*(uus)/xmio/(1.d0-sro)**2
        !
        dlt = dlw+dlo
        !
        !..... derivada do fluxo
        !
        dfw=dlw/xlt(uu) - xlw(uu)*dlt/xlt(uu)**2
        !
        dfe=(dfw*xlo(uu)+xlw(uu)/xlt(uu)*dlo)
        !
        dg=dfw*vv-(rhoo-rhow)*xg*xk*dfe
        !
    end select
    !
    end function
    !
    !=======================================================================
    !
    function dh(uu,vv,xk)
    !
    !..... calcula dh/du em um no ou volume
    !
    use mPropGeoFisica, only: xlw,xlt,xlo,gf3,iflag_linear,srw,sro,xmiw,xmio,rhoo,rhow
    !
    !
    implicit none
    real(8) :: dh,uu,vv
    real(8) :: dfw,dlt,dlo,dlw,dfe
    real(8) :: uus,zero
    real(8) :: xg,xk
    !
    zero=0.d0
    !
    !  .... Buckley-Leverett
    !
    !      g=(uu**2)/(uu**2+a*(1-uu)**2)
    !        a=1.d0/2.d0
    !        a1=uu**2+a*(1.d0-uu)**2
    !        dg=(2.d0*uu/a1 - uu**2*(2.d0*uu-2.d0*a*(1.d0-uu))/a1**2)*vv
    !
    !  .... Burgers, p. 464 do LeVeque
    !
    !       a=dsqrt(2.d0)/2.d0 ! 45 graus
    !
    !       g=a*uu**2/2.d0
    !
    !       dg=a*uu
    !
    !  .... Agua-Oleo
    !
    xg=gf3
    !
    select case(iflag_linear)
    case(1)
        !
        !..... derivada do fluxo
        !
        dh=vv
        !
    case(2)
        !
        !..... derivadas das mobilidades
        !
        uus = uu-srw
        if (uus.lt.zero) uus = zero
        dlw = 2.d0*(uus)/xmiw/(1.d0-srw)**2
        !
        uus = 1.d0-sro-uu
        if (uus.lt.zero) uus = zero
        dlo =-2.d0*(uus)/xmio/(1.d0-sro)**2
        !
        dlt = dlw+dlo
        !
        !..... derivada do fluxo
        !
        dfw=dlw/xlt(uu) - xlw(uu)*dlt/xlt(uu)**2
        !
        dfe=(dfw*xlo(uu)+xlw(uu)/xlt(uu)*dlo)
        !
        dh=dfw*vv-(rhoo-rhow)*xg*xk*dfe
        !
    end select
    !
    end function
    !
    ! =======================================================================
    !
    function dif1o(xc,xi)
    !
    !      derivada de primeira ordem
    !
    implicit none
    !
    real(8) :: dif1o,a1o,xi,xc
    !
    a1o=0.0d0
    !       a1o=1.d0
    !
    !      calculo das derivadas
    !
    dif1o= (xc-xi)*a1o
    !
    end function

    !
    !======================================================================
    !
    subroutine mass(nsd, numel,s,phi,temp)
    !
    !...  objetivo: calcula a massa total do sistema
    !
    !...  obs: atraves da fracao de volumes da fase: phi_alpha = phi*S_alpha
    !
    use mLeituraEscritaSimHidroGeoMec
    use mPropGeoFisica, only: hx, hy, hz
    !
    implicit none
    !
    integer :: numel, nsd
    real(8) :: cmass,temp
    real(8), dimension(*) :: s,phi
    !
    integer :: nel
    !
    cmass=0.d0

    do nel=1,numel
        cmass=cmass+phi(nel)*s(nel)*hx*hy
        if(nsd==3)cmass=cmass*hz
    end do
    !
    !     impressao
    !
    write(imass,"(f10.5,f25.15)") temp, cmass
    !
    end subroutine

    !=======================================================================
    !
    subroutine ktdd_rt_CORREA(nsd,nt,&
        conecLadaisElem,listaDosElemsPorFace,dt,uf,&
        v,numLadosElem,GEOTIME)
    !
    use mGlobaisEscalares, only: ns
    use mHidrodinamicaRT, only: ndofV
    use mMalha,            only: numLadosReserv, numelReserv
    use mPropGeoFisica
    !
    implicit none
    !
    integer :: nsd,nt,numLadosElem
    integer :: conecLadaisElem(numLadosElem,numelReserv)
    integer :: listaDosElemsPorFace(numLadosElem, numLadosReserv)
    real(8)   :: dt, dttotal,dtt
    real(8), dimension(numelReserv)     :: uf, phiINTERP
    real(8)   :: v(1,numLadosReserv), GEOTIME
    !
    integer :: nel
    real(8) :: umin,z,umax
    real(8), allocatable          :: du(:,:)
    real(8) :: t1
    !
    integer :: tid
    integer :: numThreads
    integer :: numPrimeiroElemento, numUltimoElemento
    integer :: inicioSol, fimSol
    integer, parameter :: meuNumMaxThreads=32
    real(8),dimension(meuNumMaxThreads)  :: uminAux,umaxAux
    !
    !.... Mensagem: KT DxD
    !
    z=0.d0
    !
    if (allocated(du) .eqv. .false.) allocate(du(nsd,numelReserv))
    du = 0.0d0
    !
    !.... identifica os volumes e nohs internos e do contorno
    !
    uf=0.d0
    call nbflux(numelReserv,uf)
    !
    tid=1
    numThreads=1
    uminAux = 1.d+9
    umaxAux =-1.d+9
    umin = 1.d+9
    umax =-1.d+9
    !
    if(tid==1) print*, "numThreads",numThreads
    !
    !.... divisao do Trabalho para a atualizacao da solucao
    !
    numPrimeiroElemento = 1
    numUltimoElemento   = numelReserv
    call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, &
        numThreads, tid-1, inicioSol, fimSol)
    !
    call timing(t1)
    !
    dttotal =0.d0
    !
    do ns=1, nt
        !
        !.... Interpolacao da porosidade
        dttotal = real(ns)*dt
        dtt     = real(ns-1)*dt
        call interpPHI(phiINTERP,phi,phi0,dtt,GEOTIME,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....    calcula as derivadas
        !
        call derivau(nsd,numLadosElem,&
            conecLadaisElem,listaDosElemsPorFace,du,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....  Calcula os fluxos e evolui no tempo
        !
        call fkt(nsd,ndofV,dt,hx,hy,hz,&
            uf,du,phiINTERP,v,numLadosElem,permkx,permky,permkz,inicioSol,fimSol)
        !
        !$OMP BARRIER
        !
        !....  Atualiza a solucao
        !
        do nel=inicioSol,fimSol
            satElem(nel)=satElemAnt(nel)
            if (satElem(nel).gt.umaxAux(tid)) umaxAux(tid) = satElem(nel)
            if (satElem(nel).lt.uminAux(tid)) uminAux(tid) = satElem(nel)
        end do
        !
        !$OMP BARRIER

        !$OMP MASTER
        !
        !....  SATURATION UPDATE WITH POWER LAW FOR TIME FRACTION
        !
        CALL SATCORRETOR(numelReserv,phi,phi0,satElem,dt,dttotal,GEOTIME)
        !
        !      write(*,'(a,i0,2x,a,3f12.5)') "passo", ns,'massa global = ',sum(satElem), maxval(satElem), minval(satElem)
        !
        !$OMP END MASTER
        !
    end do ! ns
    !
    !     calcula a derivada para a reconstrucao linear
    !
    call derivau(nsd,numLadosElem,conecLadaisElem,listaDosElemsPorFace,du,inicioSol,fimSol)
    !
    ! final da Paralelização
    !$OMP END PARALLEL
    !
    write(*,"(//,'Smax: ',f15.7,2x,'Smin: ',f15.7,//)") maxval(umaxAux), minval(uminAux)
    !
    end subroutine ktdd_rt_CORREA
    !
    !=======================================================================
    !*** NEW ***** MODIFIED SATCOR SUBROUTINE FOR GEOMECHANICS ***********
    !
    SUBROUTINE SATCORRETOR(numel,phi,phi0,s,deltat,trnstime,GEOTIME)
    !
    !     Objetivo: passo corretor para a saturacao
    !
    !----------------------------------------------------------------------
    !
    implicit none
    integer                 :: numel,NEL
    real(8), dimension(*)   :: s, phi,phi0
    real(8)                 :: trnstime, GEOTIME, deltat
    real(8)                 :: aux, aux0, aux1
    real(8)                 :: geoinv
    !...
    geoinv = 1.d0/GEOTIME
    do nel=1,numel
        !
        aux  = (phi(nel)-phi0(nel)) * geoinv
        aux0 = phi0(nel) + (trnstime-deltat) * aux
        aux1 = phi0(nel) + trnstime          * aux
        !.... CORRECAO DA SATURACAO
        !
        s(nel)=s(nel)*aux0/aux1
        !
    END DO ! nel
    !
    RETURN
    !
    END SUBROUTINE SATCORRETOR
    !
    !=======================================================================
    !*** NEW ***** MODIFIED SATCOR SUBROUTINE FOR GEOMECHANICS ***********
    !
    SUBROUTINE SATCOR(numel,phi,phi0,s,trnstime,GEOTIME)
    !
    !     Objetivo: passo corretor para a saturacao
    !
    !----------------------------------------------------------------------
    !
    implicit none
    integer                 :: numel,NEL
    real(8), dimension(*)   :: s, phi,phi0
    real(8)                 :: trnstime, GEOTIME
    real(8)                 :: aux
    !...
    aux = trnstime/GEOTIME
    !...
    do nel=1,numel
        !
        !.... CORRECAO DA SATURACAO
        !
        s(nel)=s(nel)*(phi0(nel)/phi(nel))**(trnstime/geotime)
        !
    END DO ! nel
    !
    RETURN

    END SUBROUTINE SATCOR
    !
    ! =======================================================================
    ! =======================================================================
    !
    SUBROUTINE interpPHI(phiINTERP,phi,phi0,trnstime,GEOTIME,inicio,fim)
    !
    !     Objetivo: passo corretor para a saturacao
    !
    !----------------------------------------------------------------------
    !
    implicit none
    integer                 :: inicio,fim,NEL
    real(8), dimension(*)   :: phi,phi0,phiINTERP
    real(8)                 :: trnstime, GEOTIME
    real(8)                 :: aux
    real(8)                 :: geoinv
    !...
    geoinv = 1.d0/GEOTIME
    do nel=inicio,fim
        !
        aux  = (phi(nel)-phi0(nel)) * geoinv
        phiINTERP(nel) = phi0(nel) + trnstime * aux
        !
    END DO ! nel
    !
    RETURN
    !
    END SUBROUTINE interpPHI

    end module mTransporte



