    !         programa de elementos finitos em fortran 90
    !         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
    !
    !         Eduardo Garcia e Tuane Lopes
    !         bidu@lncc.br, tuane@lncc.br
    !
    !         LNCC/MCT
    !         Petropolis, 07.2013

    !
    !**** new module *************************************************************
    !
    module mHidroDinamicaRT
    !
    implicit none
    !
    integer :: ndofV, ndofP, nlvectV, nlvectP
    real*8,  allocatable :: pressaoElem(:,:),  pressaoElemAnt(:)
    real*8,  allocatable :: velocLadal(:,:), fVeloc(:,:)
    real*8,  allocatable :: vc(:,:), ve(:,:,:)
    real*8,  allocatable :: velocNodal(:,:)
    REAL(8), ALLOCATABLE :: PRESPRODWELL(:)

    integer              :: neqV, nalhsV, nedV
    real*8,  allocatable :: alhsV(:), brhsV(:)
    integer, allocatable :: idVeloc(:,:)
    integer, allocatable :: idiagV(:)
    integer, allocatable :: lmV(:,:,:)
    integer :: numCoefPorLinhaVel, meanbwV
    character(len=7) :: optSolverV
    integer, allocatable :: AiVel(:), ApVel(:), LMstencilVel(:,:)
    INTEGER ptV(64), iparmV(64)
    REAL*8  dparmV(64)
    logical :: simetriaVel
    !
    REAL(8) :: COTAREF, PRESSAOREF
    REAL(8) :: TPRESPRODWELL,PRESMEDIAINICIAL,PRESPROD
    REAL(8) :: ZFUNDOPOCCO
    !
    INTEGER :: NPRESPRODWELL,NCONDP
    REAL(8), DIMENSION(2,400) :: PCONDP
    INTEGER, DIMENSION(400)   :: ELEM_CONDP
    !
    public  :: hidroGeomecanicaRT
    private :: montarMatrizVelocidadeRT,  limparSistEqAlgVelocidade
    private :: montarSistEqAlgVelocidade
    private :: calcularPressaoRT
    !
    contains
    !
    !**** new *******************************************************************
    !
    subroutine hidroGeomecanicaRT(satElemL,satElemL0,SIGMAT,SIGMA0,TIMEINJ)
    !
    use mGlobaisArranjos,  only : mat, c, beta
    use mGlobaisEscalares, only : tempoSolverVel,tempoMontagemVel,tempoTotalPressao
    use mGlobaisEscalares, only : dtBlocoTransp, tTransporte
    use mGlobaisEscalares, only : nvel, nnp
    use mMalha,            only : conecNodaisElem, conecLadaisElem
    use mMalha,            only : numel, numelReserv, nsd, nen
    use mMalha,            only : x, xc,numLadosReserv, numLadosElem
    use mMalha,            only : listaDosElemsPorFace
    use mPropGeoFisica,    only : phi, phi0, nelx, nely, nelz
    use mPropGeoFisica,    only : calcphi, hx, hy
    use mPropGeoFisica,    only : permkx, xkkc, xltGeo
    use mSolverPardiso,    only : solverPardisoCSR

    use mSolverGaussSkyline, only : solverGaussSkyline
    use mSolverHypre
    !
    implicit none
    !
    real*8 :: satElemL(numelReserv), satElemL0(numelReserv)
    real*8 :: SIGMAT(numelReserv), SIGMA0(numelReserv)
    !
    real*8  :: t1, t2, t3, TIMEINJ
    integer :: i, nel, arq
    character (len=80) :: nomeArqVel, nomeArqPres
    real*8 :: difP, difV, difS
    real(8) :: xlambda,xmu,bal,maxBal
    real(8), save :: tTeste
    real(8) :: elemAnt, elemDep, dx
    real(8) :: difVel, maxDifVel, minDifVel

    real*8  ::  final_res_norm, tol
    integer :: num_iterations
    !
#ifdef debug
    nomeArqVel ="velocidadeLadal_RT.txt"
    nomeArqPres="pressao_RT.txt"
    open(unit=206, file=nomeArqVel)
    open(unit=207, file=nomeArqPres)
#endif
    !
    print*, "calculando a HIDRODINAMICA"
    !
    !******************* VELOCIDADE **********************
    !
    if(.not.allocated(alhsV)) allocate(alhsV(nalhsV))
    if(.not.allocated(brhsV)) allocate(brhsV(neqV))
    !
    call limparSistEqAlgVelocidade()
    !
    call timing(t1)

    call montarSistEqAlgVelocidade(nsd,satElemL,satElemL0, &
        &     PHI,PHI0,SIGMAT,SIGMA0,TIMEINJ)
    !
    call timing(t2)
    !
#ifdef mostrarTempos
    write(*,*) ' tempo montagem do sistema da velocidade ', t2-t1
#endif
    tempoMontagemVel=tempoMontagemVel+(t2-t1)
    !
    write(*,'(a)', ADVANCE='NO') 'solucao do sistema de eq, VELOCITY, '
    call solverVel()
    call timing(t3)
    call btod(idVeloc,velocLadal,brhsV,ndofV,numLadosReserv)

#ifdef mostrarTempos
    write(*,*) ' tempo do solver velocidade', t3-t2
#endif
    tempoSolverVel=tempoSolverVel+(t3-t2)
    !
    print*, "Calculando Pressao"
    call timing(t1)
    call calcularPressaoRT (x, conecNodaisElem, conecLadaisElem,   &
        pressaoElem, pressaoElemAnt, velocLadal, &
        satElemL, satElemL0, dtBlocoTransp, SIGMAT, SIGMA0)
    !
    call calcve(numLadosElem,numelReserv,nsd,ndofV,nen,velocLadal,ve,conecLadaisElem)
    call calcvc(nsd,nen,numelReserv,vc,ve)
    call calcVelocidadeNodal(velocNodal,velocLadal,conecLadaisElem, conecNodaisElem)
    !
    call timing(t2)
    !
#ifdef mostrarTempos
    write(*,*) ' Tempo para calculo da Pressao', t3-t2
#endif
    tempoTotalPressao=tempoTotalPressao+(t2-t1)

    contains
    subroutine solverVel()
    real*8 :: solutionNorm
    integer :: i
    !
    if (optSolverV=='skyline') then
        write(*,'(2a)') ' direto ', optSolverV
        call solverGaussSkyline(alhsV,brhsV,idiagV,nalhsV, neqV, 'full')
    end if
    !
    if (optSolverV=='pardiso') then
        write(*,'(2a)') ' direto ', optSolverV
        call solverPardisoCSR(alhsV, brhsV, ApVel, AiVel, ptV, iparmV, dparmV,  &
            &        neqV, nalhsV, simetriaVel, 'vel', 'full')
    end if
    !
    if(optSolverV=='hypre') then

        write(*,'(2a)') ' iterativo ', optSolverV

        call fecharMatriz_HYPRE         (A_HYPRE_V, parcsr_A_V)
        call fecharVetor_HYPRE          (b_HYPRE_V, par_b_V   )
        call fecharVetor_HYPRE          (u_HYPRE_V, par_u_V   )

        if(.not.allocated(initialGuess_V)) then
            allocate(initialGuess_V(neqV)); initialGuess_V=0.0
        endif

        solver_id_V  = 1
        precond_id_V = 1
        tol = 1.0e-12
        call resolverSistemaAlgHYPRE (A_HYPRE_V, parcsr_A_V, b_HYPRE_V, par_b_V, u_HYPRE_V, par_u_V, &
            solver_V, solver_id_V, precond_id_V, tol,    &
            num_iterations, final_res_norm, initialGuess_V, brhsV, rows_V, neqV, myid, mpi_comm)

        call extrairValoresVetor_HYPRE(u_HYPRE_V, 1, neqV, rows_V,BRHSV)
        initialGuess_V=brhsV

        call destruirMatriz_HYPRE(A_HYPRE_V)
        call destruirVetor_HYPRE (b_HYPRE_V)
        call destruirVetor_HYPRE (u_HYPRE_V)

        call criarMatriz_HYPRE  (A_HYPRE_V, Clower_V, Cupper_V, mpi_comm )
        call criarVetor_HYPRE   (b_HYPRE_V, Clower_V, Cupper_V, mpi_comm )
        call criarVetor_HYPRE   (u_HYPRE_V, Clower_V, Cupper_V, mpi_comm )

    endif!
    write(*,*) "Valores nos extremos do vetor solucao Velocidade RT,  "
    write(*,'(5e16.8)') brhsV(1     :5)
    write(*,'(5e16.8)') brhsV(neqV-4: neqV)
    solutionNorm = 0.0
    do i = 1, neqV
        solutionNorm = solutionNorm + brhsV(i)**2
    end do
    write(*,'(a,1e16.8)') "Euclid norms of computed vel solution: ", sqrt(solutionNorm)

    end subroutine solverVel

    end subroutine hidroGeomecanicaRT
    !
    !**** new **********************************************************************
    !
    subroutine limparSistEqAlgVelocidade()
    !
    velocLadal  = 0.0d00
    if(allocated(brhsV))  brhsV = 0.0d00
    if(allocated(alhsV))  alhsV = 0.0d00
    !
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine montarSistEqAlgVelocidade(nsd,satElemL,satElemL0, &
        &           PHI,PHI0,SIGMAT,SIGMA0,TIMEINJ)
    !
    use mGlobaisEscalares,    only: dtBlocoTransp, nnp, nvel
    use mMalha,               only: x, conecNodaisElem
    use mMalha,               only: conecLadaisElem, numLadosElem
    use mMalha,               only: numLadosReserv, numelReserv, numel
    use mMalha,               only: listaDosElemsPorFace
    !
    implicit none
    !
    integer, intent(in) :: nsd
    integer :: i
    real*8  :: satElemL(numelReserv), satElemL0(numelReserv)
    REAL(8) :: PHI(numelReserv), PHI0(numelReserv)
    REAL(8) :: SIGMAT(numelReserv), SIGMA0(numelReserv)
    REAL(8) :: TIMEINJ
    !
    if (nlvectV.gt.0) CALL LOAD(idVeloc,fVeloc,brhsV,ndofV, &
        &                       numLadosReserv,nlvectV)

    if (nlvectV.gt.0) call ftodTIME(idVeloc,velocLadal,fVeloc, &
        &                       ndofV,numLadosReserv,nlvectV,TIMEINJ)
    !
    call montarMatrizVelocidadeRT(x, conecNodaisElem, &
        &     conecLadaisElem, alhsV, brhsV, idiagV, lmV, &
        &     velocLadal, pressaoElemAnt, satElemL, satElemL0, &
        &     dtBlocoTransp, PHI, PHI0, SIGMAT, SIGMA0)

    end subroutine montarSistEqAlgVelocidade
    !
    !**** new **********************************************************************
    !
    subroutine montarMatrizVelocidadeRT (x,conecNodaisElem, &
        &           conecLadaisElem, alhs, brhs,idiagV,lmV,velocLadal, &
        &           pressao, satElemL, satElemL0, dtBlocoTransp, &
        &           PHI, PHI0, SIGMAT, SIGMA0 )

    use mSolverGaussSkyline,only:addrhs, addlhs
    use mGlobaisEscalares, only: nrowsh, npint, nnp
    use mGlobaisEscalares, only: ligarBlocosHetBeta, iflag_beta
    use mGlobaisEscalares, only: S3DIM
    use mGlobaisArranjos,  only: grav, c, mat, beta
    use mfuncoesDeForma,   only: shlq, shlqrt,shgq, shgqrt
    use mfuncoesDeForma,   only: shlq3d, shg3d, shlqrt3d, shgqrt3d
    use mMalha,            only: local, nsd
    use mMalha,            only: numnp, numel, numelReserv
    use mMalha,            only: numLadosElem, numLadosReserv, nen
    use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv
    use mPropGeoFisica,    only: xkkc,xkkcGeo,xlt,xltGeo, xlo, xlw
    use mPropGeoFisica,    only: perm, permkx, permky, permkz
    use mPropGeoFisica,    only: gf1, gf2, gf3, PWELL, rhow,rhoo
    use mSolverPardiso,    only: addlhsCRS
    !..4COMPRESSIBILITY:
    use mPropGeoFisica,    only: YOUNG, POISVECT, BULK,GRAINBLK
    use mPropGeoFisica,    ONLY: BULKWATER, BULKOIL
    !        USE mMCMC,             ONLY: ELEM_CONDP
    use mPropGeoFisica,    only: hx,hy,hz
    use mSolverHypre
    !
    !.... program to calculate stifness matrix and force array for a
    !     singular problem element in one dimension
    !     form and assemble into the global left-hand-side matrix
    !                  and right-hand side vector
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8,  intent(in) :: x(nsd,numnp)
    integer, intent(in) :: conecNodaisElem(nen,numel)
    integer, intent(in) :: conecLadaisElem(numLadosELem,numelReserv)
    real(8), intent(inout) :: alhs(nalhsV), brhs(neqV)
    integer, intent(in) :: idiagV(neqV)
    integer, intent(in) :: lmV(ndofV,numLadosElem,numelReserv)
    real*8,  intent(inout) :: velocLadal(ndofV,numLadosReserv)
    !old      real*8,  intent(in)    :: pressao(numelReserv), sw(numelReserv)
    !Next 2 Lines Modified iterative Hidrodinamics and Transport
    real*8,  intent(in) :: pressao(numelReserv)
    real*8,  intent(in) :: satElemL(numelReserv)
    real*8,  intent(in) :: satElemL0(numelReserv)
    real*8,  intent(in) :: dtBlocoTransp
    !4COMPRESSIBILITY:
    real*8,  intent(in) :: PHI(numelReserv), PHI0(numelReserv)
    real*8,  intent(in) :: SIGMAT(numelReserv), SIGMA0(numelReserv)
    !
    real*8 :: xl(nsd,nen), dl(ndofV,nen)
    real*8 :: fluxl(ndofV,nen), vnl(nsd,numLadosElem)
    real*8 :: shgrt(nrowsh,numLadosElem,npint)
    real*8 :: shlrt(nrowsh,numLadosElem,npint)
    real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
    real*8 :: det(npint), w(npint)
    !
    real*8  :: elfrt(numLadosElem), elert(numLadosElem,numLadosElem)
    integer :: nee, neesq
    !
    integer :: nel, m
    integer :: l, i, j, k
    real*8  :: pi
    real*8  :: pix, piy, piz, sx, sy, sz, cx, cy, cz
    !
    real(8) :: xk,yk,zk,delta,theta
    real(8) :: h,h2,xx,yy,zz,du,dux,duy,duz,xindi,xindj,ff
    real(8) :: c1,phii1,phii2,phii3,phij1,phij2,phij3
    real(8) :: a11,a22,a33,divphij,divphii
    real(8) :: swint,pp
    real(8) :: dt
    !
    integer :: tid, omp_get_thread_num,omp_get_num_threads
    integer :: numPrimeiroElemento, numUltimoElemento
    integer :: numThreads, inicioSol, fimSol
    !
    logical :: diag,zerol, quad,lsym
    real(8) :: xxlt,xxlo,xxlw,lambdab
    !...
    !... BEGIN CATALOG FOR ITERATIVE FORMULATION FOR COMPRESSIBILTY
    !
    REAL(8) :: BULKROCK, BETACOM, BIOTCOEF, BETAT, BETAT0, BETA4P
    !
    REAL(8) :: ALAM, AMU2, XBULK, POISSON, BETAGEO, COEFSIGM
    REAL(8) :: DIFSIGMA, AT
    INTEGER :: NELWELL
    REAL(8) :: XNELWELL,PRESAUX
    !
    !... END CATALOG FOR FOR ITERATIVE FORMULATION
    !
    nee = numLadosElem*ndofV; neesq = nee*nee
    shlrt=0.d0
    tid=1
    numThreads=1
    !
    diag = .false.
    quad = .true.
    pi   = 4.0d0*datan(1.0d0)
    dt   = dtBlocoTransp
    at = dfloat(nnp)*dt
    !
    if(nsd==2) then
        call shlq  (shl,w,npint,nen)
        call shlqrt(numLadosElem,npint,w,shlrt)
    else
        call shlq3d  (shl,w,npint,nen)
        call shlqrt3d(numLadosElem,npint,w,shlrt)
    endif
    !
    ! $OMP PARALLEL FIRSTPRIVATE(tid) &
    ! $OMP PRIVATE (numPrimeiroElemento, numUltimoElemento, inicioSol,fimSol) &
    ! $OMP PRIVATE (xxlw,xxlo,xxlt,lambdab,xk,yk,h,h2,hx,hy,delta,theta, xl,fluxl,shg,shgrt,a11,a22,a33) &
    ! $OMP PRIVATE (divphij,divphii ,phii1, phii2, phij1, phij2, AMU2, ALAM, XBULK,det,c1,xindi,xindj,vnl,nel,i,j,l,pp,m) &
    ! $OMP REDUCTION(+:elert,elfrt) !,brhs,alhs)
    !
    ! #ifdef withOMP
    !         tid=tid+omp_get_thread_num()
    !         numThreads=omp_get_num_threads()
    ! #endif
    !
    !       if(tid==1) print*, "Em Velocidade, numThreads=",numThreads
    !
    numPrimeiroElemento = 1
    numUltimoElemento   = numelReserv
    call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, &
        &     numThreads, tid-1, inicioSol, fimSol)
    !
    !rt   calcula as normais externas locais
    !
    vnl      = 0.d0
    vnl(2,1) =-1.d0
    vnl(1,2) = 1.d0
    vnl(2,3) = 1.d0
    vnl(1,4) =-1.d0
    if (nsd==3) then
        vnl(3,5)=-1.d0
        vnl(3,6)= 1.d0
    endif
    !!
    do 500 nel=inicioSol,fimSol
        !
        !....    clear stiffness matrix and force array
        elfrt=0.0d00; elert=0.0d00
        !
        !....    localize coordinates and Dirichlet b.c.
        !
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
        call local(conecLadaisElem(1,nel),velocLadal,fluxl, &
            &        numLadosElem,nedV,nedV)
        !
        m = mat(nel)
        quad = .true.
        !....    if (nen.eq.4.and.conecLadaisElem(3,nel).eq.conecLadaisElem(4,nel)) quad = .false.
        !
        !....    chama a shg. eh necessario pra calcular a fonte e o determinante
        !
        if (nsd==2) then
            call shgq (xl,det,shl,shg,npint,nel,quad,nen)
        else
            call shg3d(xl,det,shl,shg,npint,nel,nen)
        endif
        !
        !....  length of the element
        h2 = 0
        h2 = h2 + (xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,4))**2
        if (nsd==3) h2 = h2 + (xl(3,1)-xl(3,5))**2
        h  = dsqrt(h2)/2.d00
        h2 = h*h
        !
        if (nsd==2) then
            call shgqrt(numLadosElem,npint,hx,hy,shlrt,shgrt)
        else
            call shgqrt3d(numLadosElem,npint,hx,hy,hz,shlrt,shgrt)
        endif
        !
        !..... form stiffness matrix
        !
        pp       = pressao(nel)
        !
        DIFSIGMA = (SIGMAT(NEL)-SIGMA0(NEL))/S3DIM
        !
        !...  COMPRESSIBILIDADE EQUIVALENTE ESCOAMENTO BIFASICO
        !
        theta=1.0d0
        !
        BULKROCK = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
        BETAGEO  = 1.0D0/BULKROCK
        BIOTCOEF = 1.0D0 - BULKROCK/GRAINBLK(1)
        BETAT0   = BETASTAR(PHI(NEL),BIOTCOEF,BETAGEO,satElemL0(NEL))
        BETAT    = BETASTAR(PHI(NEL),BIOTCOEF,BETAGEO,satElemL(NEL))
        BETA4P   = BETAT0/BETAT
        COEFSIGM = BIOTCOEF*BETAGEO/BETAT


        !         write(*,*) nel,difsigma,sigmat(nel),sigma0(nel)!betageo,betat0,betat,biotcoef
        !
        do l=1,npint
            c1 = w(l)*det(l)
            !
            !.... calcula x, y, a pressao e suas derivadas no ponto de integracao
            !
            xx  = 0.d0
            yy  = 0.d0
            if (nsd==3) zz = 0.d0
            du  = 0.d0
            dux = 0.d0
            duy = 0.d0
            if (nsd==3) duz = 0.0D0
            !
            !.... swint=0.d0
            !
            do i=1,nen
                xx = xx + shl(nrowsh,i,l)*xl(1,i)
                yy = yy + shl(nrowsh,i,l)*xl(2,i)
                if (nsd==3) zz = zz + shl(nrowsh,i,l)*xl(3,i)
                du  = du  + shl(nrowsh,i,l)*dl(1,i)
                dux = dux + shg(1,i,l)*dl(1,i)
                duy = duy + shg(2,i,l)*dl(1,i)
                if (nsd==3) duz = duz + shg(3,i,l)*dl(1,i)
                !swint=swint+shg(3,i,l)*se(i,nel)
            end do
            !
            xxlw = xlw(satElemL(nel))
            xxlo = xlo(satElemL(nel))
            xxlt = xltGeo(satElemL(nel))
            !
            lambdab = (xxlw*rhow+xxlo*rhoo)/xxlt
            !
            !.... CALCULANDO A PERMEABILIDADE FUNCCAO DA POROSIDADE EULERIANA
            !
            XK  = XKKCGEO(PHI0(NEL),PERMKX(NEL))*XXLT
            A11 = 1.0D0/XK
            !
            YK  = XKKCGEO(PHI0(NEL),PERMKY(NEL))*XXLT
            A22 = 1.0D0/YK
            !
            IF (NSD==3) THEN
                ZK  = XKKCGEO(PHI0(NEL),PERMKZ(NEL))*XXLT
                A33 = 1.0D0/ZK
            ENDIF
            !.. OLD:::  parametros do penalty:  delta=dt/beta
            !
            DELTA = DT/BETAT
            !
            !.... vetor de carga - RHS - f  =  gf0
            !
            gf1 = grav(1)
            gf2 = grav(2)
            if (nsd==3) gf3 = grav(3)
            !
            pix = pi*xx
            piy = pi*yy
            if (nsd==3) piz = pi*zz
            !
            sx = dsin(pix)
            sy = dsin(piy)
            if (nsd==3) sz = dsin(piz)
            cx = dcos(pix)
            cy = dcos(piy)
            if (nsd==3) cz = dcos(piz)
            !
            !.... fonte NO BALANCCO DE MASSA DO FLUIDO: ff
            !
            ff = 0.d0 !gf1*sx*sy+gf2*cx*cy
            !
            do j=1,numLadosElem
                xindj = vnl(1,j)+vnl(2,j)
                if (nsd==3) xindj = xindj +vnl(3,j)
                phij1   = shgrt(1,j,l)*c1*xindj
                phij2   = shgrt(2,j,l)*c1*xindj
                if (nsd==3) phij3 = shgrt(3,j,l)*c1*xindj
                divphij = shgrt(nrowsh,j,l)*c1*xindj
                !
                !.... termo de fonte
                !
                elfrt(j) = elfrt(j) & ! falta a integral no bordo
                &                  + BETA4P*pp*divphij & ! efeito da pressao anterior
                &                  + delta*ff*divphij  &  ! fonte de massa
                &                  + lambdab*(gf1*phij1+gf2*phij2) ! fonte gravidade
                if (nsd==3) elfrt(j) = elfrt(j)+lambdab*gf3*phij3
                !
                elfrt(j) = elfrt(j)- COEFSIGM*DIFSIGMA*DIVPHIJ
                !
                !.... loop nos lados: indice i
                !
                do i=1,numLadosElem
                    xindi = vnl(1,i)+vnl(2,i)
                    if (nsd==3) xindi = xindi+vnl(3,i)
                    phii1 = shgrt(1,i,l)*xindi
                    phii2 = shgrt(2,i,l)*xindi
                    if (nsd==3) phii3 = shgrt(3,i,l)*xindi
                    divphii = shgrt(nrowsh,i,l)*xindi
                    !
                    !.... ...  matriz de rigidez
                    !
                    elert(i,j) = elert(i,j) + a11*phii1*phij1  &
                        + a22*phii2*phij2  &
                        + delta*theta*divphii*divphij
                    if (nsd==3) elert(i,j) = elert(i,j)+a33*phii3*phij3
                end do ! numLadosElem,i
                !
            end do ! numLadosElem,j
            !
        end do ! nint,l
        !
        !.... computation of Dirichlet b.c. contribution
        !.... CALCULA A PRESSAO NO POCCO DE PRODUCAO A PARTIR DO VALOR DA PRESSAO
        !     NO FUNDO DO POCCO E LEVANDO EM CONSIDERACAO A GRAVIDADE
        IF (NSD.EQ.2) THEN
            DO NELWELL=1,NCONDP
                IF(NEL.EQ.ELEM_CONDP(NELWELL))THEN
                    CALL DIRICHLET_POCO_GRAVIDADE(PRESPROD,NEL,HY,NCONDP, &
                        GF2,NELWELL,satElemL,pressao)
                    XINDJ = VNL(1,2)+VNL(2,2)
                    ELFRT(2)=ELFRT(2)-(PRESPRODWELL(NELWELL))*XINDJ*(HY)
                    EXIT
                END IF
            END DO
            !             IF (MOD(NEL,nelxReserv).EQ.0) THEN
            !                XINDJ = VNL(1,2)+VNL(2,2)
            !                K = NEL/nelxReserv
            !!                 write(2028,4500) NEL/nelxReserv, PWELL(NEL/nelxReserv)
            !                ELFRT(2)=ELFRT(2)-deplet(at,PWELL(K))*XINDJ*HY
            !             END IF
        ENDIF
        !
        call ztest(fluxl,nee,zerol)

        if(.not.zerol) then
            call kdbc(elert,elfrt,fluxl,nee)
        end if
        !
        lsym=.true.

        if (optSolverV=='skyline')   then
            call addlhs(alhs,elert,lmV(1,1,nel),idiagV,nee,diag,lsym)
        endif
        if (optSolverV=='pardiso')   then
            call addlhsCRS(alhs,elert,lmV(1,1,nel),ApVel, AiVel,nee)
        endif
        if (optSolverV=='hypre')   then
            call addnslHYPRE(A_HYPRE_V, elert, LMV(1,1,nel), nee, lsym)
        endif

        call addrhs(brhs,elfrt,lmV(1,1,nel),nee)
        !
500 continue

    ! $OMP END PARALLEL

    if (optSolverV=='hypre')   then
        do i = 1, neqV
            rows_V(i) = i-1
        end do
        call atribuirValoresVetor_HYPRE(b_HYPRE_V, 1, neqV, rows_V, BRHS)
    endif

    RETURN
4500 FORMAT(I8,X,40(1PE15.8,2X))
    !
    end subroutine
    !*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
    !
    FUNCTION deplet(tempo, pinit)
    !
    use mGlobaisEscalares, only: tt
    !
    !... COMPUTE LINEAR FACTOR TO DEPLET PRESSURE
    !
    IMPLICIT NONE
    !
    REAL(8) :: tempo, deplet, depletlin, pinit
    !
    !      IF (TEMPO.LT.0.5D0*TT) THEN
    !            DEPLETLIN = pinit !  - 0.5D0*TEMPO/TT
    !         ELSE
    !       DEPLETLIN = pinit*(1.0D0-0.183333D0*(tempo/tt))
    !      ENDIF
    !
    !      DEPLETLIN = pinit - 7584.233D0*TEMPO
    DEPLET =  pinit  !0.0d0 !1.0D0 !   DEPLETLIN !
    !
    END FUNCTION
    !
    !*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
    !
    FUNCTION deplet_old(tempo)
    !
    use mGlobaisEscalares, only: tt
    !
    !... COMPUTE LINEAR FACTOR TO DEPLET PRESSURE
    !
    IMPLICIT NONE
    !
    REAL(8) :: tempo, deplet_old, depletlin
    !
    IF (TEMPO.LT.0.5D0*TT) THEN
        DEPLETLIN = 1.0D0 !  - 0.5D0*TEMPO/TT
    ELSE
        DEPLETLIN = 0.01D0
    ENDIF
    !
    !      IF (TEMPO.LT.0.5D0*tt) THEN
    !            depletlin = 1.0d0 - 0.5d0*tempo/tt
    !         ELSE
    !            depletlin = 0.75D0
    !      ENDIF
    !
    deplet_old =  DEPLETLIN ! 0.0d0  !
    !
    END FUNCTION
    !
    !*** NEW ***** FUNCTION TO COMPUTE RECIPROCAL OF BIOT MODULUS ***********
    !
    FUNCTION BIOTMOD(POROSITY,ALPHA,SAT)
    !
    !... COMPUTE RECIPROCAL OF BIOT MODULUS,THUS 1/M: REFERENCE CHENG.PDF
    !
    use mPropGeoFisica,    ONLY: BULKWATER, BULKOIL, GRAINBLK
    use mPropGeoFisica,    ONLY: RHOW, RHOO
    !
    IMPLICIT NONE
    !
    REAL(8) :: ALPHA, SATDENSW, SATDENSO, DENSITEQ
    REAL(8) :: BIOTMOD, POROSITY, SAT, BULKFLUID
    REAL(8) :: BIOTAGUA, BIOTOLEO, ONEOVERN
    !
    !new
    ONEOVERN = (ALPHA-POROSITY)/GRAINBLK(1)
    !new
    BIOTAGUA = POROSITY/BULKWATER + ONEOVERN
    !new
    BIOTOLEO = POROSITY/BULKOIL   + ONEOVERN
    !
    !new
    BIOTMOD = SAT*BIOTAGUA + (1.0D0-SAT)*BIOTOLEO
    !
    !old      SATDENSW  = SAT*RHOW
    !old      SATDENSO  = (1.0D0-SAT)*RHOO
    !old      DENSITEQ  = SATDENSW + SATDENSO
    !
    !old      BULKFLUID = (SATDENSW/BULKWATER + SATDENSO/BULKOIL)/DENSITEQ
    !
    !old      BIOTMOD   = POROSITY*BULKFLUID + (ALPHA-POROSITY)/BULKSOLID
    !
    END FUNCTION
    !
    !*** NEW ***** FUNCTION TO COMPUTE BETA_STAR COMPRESSIBILITY **********
    !
    FUNCTION BETASTAR(POROSITY,ALPHA,BETA,SAT)
    !
    !... COMPUTE BETA_STAR COMPRESSIBILITY REFERENCE TUTORIAL.PDF
    !
    IMPLICIT NONE
    !
    !      REAL(8) :: BIOTMOD, ALPHA, BETA, BETASTAR, POROSITY, SAT
    REAL(8) :: ALPHA, BETA
    REAL(8) :: BETASTAR, POROSITY, SAT
    !
    BETASTAR = BIOTMOD(POROSITY,ALPHA,SAT) + BETA*ALPHA**2
    !
    END FUNCTION
    !
    !**** new *************************************************************
    !
    !old      subroutine calcularPressaoRT (x, conecNodaisElem, conecLadaisElem,  &
    !old                          pressaoElem, pressaoElemAnt, velocLadal, satElem, dtBlocoTransp,         &
    !old                          SIGMAT,SIGMA0)
    !
    !Next Line Modified for Iterative Hidrodinamic and Transport
    !
    subroutine calcularPressaoRT (x, conecNodaisElem, &
        &           conecLadaisElem, pressaoElem, pressaoElemAnt,&
        &           velocLadal, satElemL, satElemL0, dtBlocoTransp,&
        &           SIGMAT,SIGMA0)
    !
    !.... program to calculate stifness matrix and force array for a
    !     singular problem element in one dimension
    !     form and assemble into the global left-hand-side matrix
    !                  and right-hand side vector
    !
    use mFuncoesDeForma,   only: shlqrt, shgqrt, shlqrt3d, shgqrt3d
    use mMalha,            only: local, nsd
    use mMalha,            only: numnp, numel, numelReserv
    use mMalha,            only: nen, numLadosElem, numLadosReserv
    use mGlobaisArranjos,  only: c, mat, beta
    use mGlobaisEscalares, only: nrowsh, iflag_beta
    use mGlobaisEscalares, only: npint, ligarBlocosHetBeta
    use mGlobaisEscalares, only: S3DIM
    use mPropGeoFisica,    only: YOUNG, POISVECT, GRAINBLK, phi, BULK
    use mPropGeoFisica,    only: hx,hy,hz
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8,  intent(in) :: x(nsd,numnp)
    integer, intent(in) :: conecNodaisElem(nen,numel)
    integer, intent(in) :: conecLadaisElem(numLadosElem,numelReserv)
    real*8,  intent(inout) :: pressaoElem(ndofP,numelReserv)
    real*8,  intent(inout) :: pressaoElemAnt(numelReserv)
    real*8,  intent(in) :: velocLadal(ndofV,numLadosReserv)
    real*8,  intent(in) :: satElemL(numelReserv)
    real*8,  intent(in) :: satElemL0(numelReserv)
    real*8,  intent(in) :: dtBlocoTransp
    real*8,  intent(in) :: SIGMAT(numelReserv), SIGMA0(numelReserv)
    !
    real*8 :: xl(nsd,nen)
    real*8, dimension(nrowsh,numLadosElem,npint) :: shgrt, shlrt
    real*8  :: w(npint)
    real(8) :: dt
    !
    integer :: nel,m,i,l
    real(8) :: vnl(nsd,numLadosElem)
    real(8) :: xind,div,fe
    real(8) :: pi,delta,gf1,gf2,gf3
    real(8) :: sx,sy,sz,cx,cy,cz,pix,piy,piz,xg,yg,zg,ff
    real(8) :: dif, aux
    !...
    !... BEGIN CATALOG FOR ITERATIVE FORMULATION
    !
    REAL(8) :: XBULK, POISSON, XTRACCO, BULKROCK, BETACOM, BIOTCOEF
    !
    REAL(8) :: BETAT,BETAT0,BETA4P,COEFSIGM,BETAGEO,DIFSIGMA
    !
    !... END CATALOG FOR FOR ITERATIVE FORMULATION
    !
    dt=dtBlocoTransp
    !
    if(nsd==2) then
        call shlqrt(numLadosElem,npint,w,shlrt)
    else
        call shlqrt3d(numLadosElem,npint,w,shlrt)
    endif
    !
    ! vetores normais
    vnl      = 0.d0
    vnl(2,1) =-1.d0
    vnl(1,2) = 1.d0
    vnl(2,3) = 1.d0
    vnl(1,4) =-1.d0
    if (nsd==3) then
        vnl(3,5)=-1.d0
        vnl(3,6)= 1.d0
    endif
    !
    if (nsd==2) then
        call shgqrt(nen,npint,hx,hy,shlrt,shgrt)
    else
        call shgqrt3d(numLadosElem,npint,hx,hy,hz,shlrt,shgrt)
    endif
    !
    do nel=1,numelReserv
        !
        !     Numeracao local das faces
        !
        !                 3
        !             ________
        !            /  6    /|
        !           /_______/ |
        !           |       |2|
        !         4 |   1   | /
        !           |_______|/
        !
        !               5
        !
        !
        !...  fonte
        !
        FF = 0.0D0 !gf1*sx*sy+gf2*cx*cy
        !
        DIFSIGMA=(SIGMAT(NEL)-SIGMA0(NEL))/S3DIM
        !
        !.... COMPRESSIBILIDADE EQUIVALENTE ESCOAMENTO BIFASICO
        !
        BULKROCK = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
        BETAGEO  = 1.0D0/BULKROCK
        BIOTCOEF = 1.0D0 - BULKROCK/GRAINBLK(1)
        BETAT0   = BETASTAR(PHI(NEL),BIOTCOEF,BETAGEO,satElemL0(NEL))
        BETAT    = BETASTAR(PHI(NEL),BIOTCOEF,BETAGEO,satElemL(NEL))
        BETA4P   = BETAT0/BETAT
        COEFSIGM = BIOTCOEF*BETAGEO/BETAT
        !
        !.. :::  parametros do penalty   OLD_LINE:   delta=dt/beta
        !
        DELTA = DT/BETAT
        !
        div = 0.d0
        do i=1,numLadosElem
            xind = vnl(1,i)+vnl(2,i)
            if (nsd==3) xind = xind + vnl(3,i)
            fe  = velocLadal(1,conecLadaisElem(i,nel))*xind
            div = div + shgrt(nrowsh,i,1)*fe
        end do ! numLadosElem
        !
        !.... calculo da pressao sobre o elemento
        !
        ff = ff-div
        !
        XTRACCO=BETA4P*pressaoElemAnt(NEL) - COEFSIGM*DIFSIGMA
        pressaoElem(1,nel)=RK(XTRACCO,DELTA,FF)
        !
    end do ! nel
    !
    return
    !
    end  subroutine calcularPressaoRT
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
    !     Evolucao no tempo (Runge-Kutta)
    !
    !     primeira ordem
    !
    rk=uc + dt*ff
    !
    return
    !
    end function rk
    !
    !=======================================================================
    !
    subroutine calcvc(nsd,nen,numel,v,ve)
    !
    !     calcula a velocidade no centro dos elementos
    !
    implicit none
    !
    integer  :: nen,numel,nel,no,nsd
    real(8)  :: vx,vy,vz
    real(8), dimension (nsd,*) :: v
    real(8), dimension (nsd,nen,*) :: ve
    !
    do nel=1,numel
        !
        vx=0.d0
        vy=0.d0
        vz=0.d0
        !
        do no=1,nen
            vx=vx+ve(1,no,nel)
            vy=vy+ve(2,no,nel)
            if(nsd==3) vz=vz+ve(3,no,nel)
        end do
        !
        v(1,nel)=vx/nen
        v(2,nel)=vy/nen
        if(nsd==3) v(3,nel)=vz/nen
        !
    end do
    !
    end subroutine
    !
    !============================================================================
    !
    subroutine calcVelocidadeNodal(velocNodal, velocLadal, conecLadaisElem, conecNodaisElem)
    !
    !     calcula a velocidade nos nos por elemento

    use mMalha, only: nen, nsd, numel, numelReserv, numLadosElem, numnpReserv, numLadosReserv
    !
    implicit none
    !
    real*8 :: velocNodal(nsd,numnpReserv), velocLadal(1, numLadosReserv)
    integer :: conecLadaisElem(numLadosElem, numelReserv)
    integer :: conecNodaisElem(nen, numel)
    !
    integer :: nel,no1,no2,no3,no4,no5,no6,no7,no8,j,i
    integer :: lado1, lado2, lado3, lado4, lado5, lado6
    !
    !.....  loop on elements
    !
    do nel=1,numelReserv
        !
        lado1=conecLadaisElem(1,nel)
        lado2=conecLadaisElem(2,nel)
        lado3=conecLadaisElem(3,nel)
        lado4=conecLadaisElem(4,nel)
        if(nsd==3)then
            lado5=conecLadaisElem(5,nel)
            lado6=conecLadaisElem(6,nel)
        endif

        no1=conecNodaisElem(1,nel)
        no2=conecNodaisElem(2,nel)
        no3=conecNodaisElem(3,nel)
        no4=conecNodaisElem(4,nel)
        if(nsd==3) then
            no5=conecNodaisElem(5,nel)
            no6=conecNodaisElem(6,nel)
            no7=conecNodaisElem(7,nel)
            no8=conecNodaisElem(8,nel)
        endif

        !
        if(nsd==2) then
            velocNodal(1,no1) = velocLadal(1,lado4)
            velocNodal(2,no1) = velocLadal(1,lado1)

            velocNodal(1,no2) = velocLadal(1,lado2)
            velocNodal(2,no2) = velocLadal(1,lado1)

            velocNodal(1,no3) = velocLadal(1,lado2)
            velocNodal(2,no3) = velocLadal(1,lado3)

            velocNodal(1,no4) = velocLadal(1,lado4)
            velocNodal(2,no4) = velocLadal(1,lado3)
        else
            velocNodal(1,no1) = velocLadal(1,lado4)
            velocNodal(2,no1) = velocLadal(1,lado1)
            velocNodal(3,no1) = velocLadal(1,lado5)

            velocNodal(1,no2) = velocLadal(1,lado2)
            velocNodal(2,no2) = velocLadal(1,lado1)
            velocNodal(3,no2) = velocLadal(1,lado5)

            velocNodal(1,no3) = velocLadal(1,lado2)
            velocNodal(2,no3) = velocLadal(1,lado3)
            velocNodal(3,no3) = velocLadal(1,lado5)

            velocNodal(1,no4) = velocLadal(1,lado4)
            velocNodal(2,no4) = velocLadal(1,lado3)
            velocNodal(3,no4) = velocLadal(1,lado5)

            velocNodal(1,no5) = velocLadal(1,lado4)
            velocNodal(2,no5) = velocLadal(1,lado1)
            velocNodal(3,no5) = velocLadal(1,lado6)

            velocNodal(1,no6) = velocLadal(1,lado2)
            velocNodal(2,no6) = velocLadal(1,lado1)
            velocNodal(3,no6) = velocLadal(1,lado6)

            velocNodal(1,no7) = velocLadal(1,lado2)
            velocNodal(2,no7) = velocLadal(1,lado3)
            velocNodal(3,no7) = velocLadal(1,lado6)

            velocNodal(1,no8) = velocLadal(1,lado4)
            velocNodal(2,no8) = velocLadal(1,lado3)
            velocNodal(3,no8) = velocLadal(1,lado6)

        endif
        !
    end do

    !       do i=1, numnpReserv
    !       write(233,*) i, velocNodal(1:2,i)
    !       end do

    !
    return
    end subroutine
    !
    !============================================================================
    !
    subroutine calcve(nedg,numel,nsd,nedfl,nen,fluxo,vel,ieedg)
    !
    !     calcula a velocidade nos nos por elemento

    use mGlobaisEscalares, only: nnp
    !
    implicit none
    !
    integer :: numel,nsd,nen,nedfl,nedg
    real(8), dimension(nedfl,*)   :: fluxo
    real(8), dimension(nsd,nen,*) :: vel
    integer, dimension(nedg,*) :: ieedg
    !
    integer :: nel,n1,n2,n3,n4,n5,n6,j
    !
    !.....  loop on elements
    !
    do nel=1,numel
        !
        n1=ieedg(1,nel)
        n2=ieedg(2,nel)
        n3=ieedg(3,nel)
        n4=ieedg(4,nel)
        if(nsd==3) then
            n5=ieedg(5,nel)
            n6=ieedg(6,nel)
        endif
        !
        if(nsd==2) then
            vel(1,1,nel) = fluxo(1,n4)
            vel(2,1,nel) = fluxo(1,n1)
            vel(1,2,nel) = fluxo(1,n2)
            vel(2,2,nel) = fluxo(1,n1)
            vel(1,3,nel) = fluxo(1,n2)
            vel(2,3,nel) = fluxo(1,n3)
            vel(1,4,nel) = fluxo(1,n4)
            vel(2,4,nel) = fluxo(1,n3)
        endif
        !
        if(nsd==3) then
            vel(1,1,nel) = fluxo(1,n4)
            vel(2,1,nel) = fluxo(1,n1)
            vel(3,1,nel) = fluxo(1,n5)

            vel(1,2,nel) = fluxo(1,n2)
            vel(2,2,nel) = fluxo(1,n1)
            vel(3,2,nel) = fluxo(1,n5)

            vel(1,3,nel) = fluxo(1,n2)
            vel(2,3,nel) = fluxo(1,n3)
            vel(3,3,nel) = fluxo(1,n5)

            vel(1,4,nel) = fluxo(1,n4)
            vel(2,4,nel) = fluxo(1,n3)
            vel(3,4,nel) = fluxo(1,n5)

            vel(1,5,nel) = fluxo(1,n4)
            vel(2,5,nel) = fluxo(1,n1)
            vel(3,5,nel) = fluxo(1,n6)

            vel(1,6,nel) = fluxo(1,n2)
            vel(2,6,nel) = fluxo(1,n1)
            vel(3,6,nel) = fluxo(1,n6)

            vel(1,7,nel) = fluxo(1,n2)
            vel(2,7,nel) = fluxo(1,n3)
            vel(3,7,nel) = fluxo(1,n6)

            vel(1,8,nel) = fluxo(1,n4)
            vel(2,8,nel) = fluxo(1,n3)
            vel(3,8,nel) = fluxo(1,n6)

        endif
        !
    end do

1001 format('nnp=',i2,1x,'elemnt 602 lado=',i1,1x,'vel_x=',1PE15.8,2X,'vel_y=',1PE15.8)

    !
    return
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE DIRICHLET_POCO_GRAVIDADE_APROX(BWP,NEL,HYY,NP,GRAVI,NELWELL,satElemL)
    !
    !        use mMCMC,          only : ELEM_CONDP
    use mPropGeoFisica, only : RHOW, RHOO
    use mMalha,         only : numelReserv
    !
    IMPLICIT NONE
    !
    REAL(8)               :: HYY,GRAVI,BWP,RHO
    INTEGER               :: NP,NELWELL,NEL,I
    real*8,  intent(in)   :: satElemL(numelReserv)
    !
    RHO=0.0
    DO I=1,NP
        RHO = RHO + RHOW*satElemL(ELEM_CONDP(I)) +  &
            RHOO*(1.0-satElemL(ELEM_CONDP(I)))
    END DO
    RHO = RHO/REAL(NP)
    PRESPRODWELL(NELWELL) = BWP+(REAL(NELWELL-1))*HYY*GRAVI*RHO
    RETURN
    END SUBROUTINE DIRICHLET_POCO_GRAVIDADE_APROX
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    SUBROUTINE DIRICHLET_POCO_GRAVIDADE(BHP,NEL,HYY,NP,GRAVI,NELWELL,satElemL,pressao)
    !
    !        use mMCMC,             only : ELEM_CONDP
    use mPropGeoFisica,    only : RHOW, RHOO, BULKWATER, BULKOIL
    use mMalha,            only : xc,numelReserv
    !
    IMPLICIT NONE
    !
    REAL(8) :: HYY,GRAVI,BHP,CW,CO,A,B,DELTAZ,RHO
    REAL(8) :: RHOW_BHP,RHOO_BHP,SAT_MEAN,AUX,RHOM
    INTEGER :: NELWELL,NEL,NP,I
    real*8,  intent(in)   :: satElemL(numelReserv)
    real*8,  intent(in)   :: pressao(numelReserv)
    !
    AUX = xc(2,nel) -HYY/2.0d0
    SAT_MEAN = 0.0
    DO I=1,NP
        !           SAT_MEAN = 1.0
        SAT_MEAN = SAT_MEAN +satElemL(ELEM_CONDP(I))
    END DO
    SAT_MEAN = SAT_MEAN/REAL(NP)
    CW=1.0D0/BULKWATER
    CO=1.0D0/BULKOIL
    RHOW_BHP = RHOW*(1.0+CW*(BHP-PRESSAOREF))
    RHOO_BHP = RHOO*(1.0+CO*(BHP-PRESSAOREF))
    A = CW*RHOW_BHP*SAT_MEAN + CO*RHOO_BHP*(1.0-SAT_MEAN)
    RHO = A*GRAVI
    RHOM = RHOW_BHP*SAT_MEAN + RHOO_BHP*(1.0-SAT_MEAN)
    B = (-A*PRESSAOREF + RHOM)*GRAVI
    !        PRESPRODWELL(NELWELL) = (RHOM/A)*(EXP(RHO*(AUX-COTAREF))-1.0)+PRESSAOREF
    PRESPRODWELL(NELWELL) = (RHOM/A)*(EXP(RHO*(AUX-ZFUNDOPOCCO))-1.0)+BHP
    !        write(*,*)xc(2,nel),PRESPRODWELL(NELWELL),pressao(NEL)
    !
    RETURN
    !
    END SUBROUTINE DIRICHLET_POCO_GRAVIDADE
    !
    !=======================================================================
    !
    SUBROUTINE leituraCoordenadasPoco(N,PCOND)

    use mInputReader, only: findKeyword,file_lines

    !
    implicit none
    !
    integer  :: N
    REAL(8), DIMENSION(2,100) :: PCOND
    !
    integer                        :: idata,i,j,ISTAT
    real(8)                        :: TOL

    integer*4 ::  keyword_line
    character(len=50) keyword_name
    !
    keyword_name = "coordenadas_poco"
    keyword_line = findKeyword(keyword_name)
    if (keyword_line.eq.0) then
        return
    endif

    TOL=1.e-6
    !
    read (file_lines(keyword_line:),  *)  N
    keyword_line = keyword_line + 1
    write(*,101)N

    !
    do i=1,N
        read (file_lines(keyword_line:),  *) (pcond(j,i),j=1,2)
        keyword_line = keyword_line + 1
        write(*,201)i,pcond(1,i),pcond(2,i)
    enddo
    !
    do i=1,N
        do j=1,2
            pcond(j,i)=pcond(j,i)+TOL
        enddo
    enddo

    !
100 format(I5)
101 format('Numero de Pontos de Controle:',i5,/)
200 format(3F12.7)
201 format('Coordenada do Ponto',i5,' (',F10.4,',',F10.4,')')

    !
    end SUBROUTINE

    !
    !===========================================================================
    !
    subroutine lerParametrosHidrodinamica_DS
    !
    use mGlobaisEscalares
    use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
    !
    IMPLICIT NONE
    !
    INTEGER I, ISimulator
    !
    CHARACTER(LEN=128) :: NAMEFILE, TEXT, FLAG1
    character(len=50)  :: keyword_name
    integer :: ierr
    !
    keyword_name='pressao_de_referencia'
    call readRealKeywordValue(keyword_name, PRESSAOREF, PRESSAOREF, ierr)
    IF(ierr.ne.1) THEN
        WRITE(*,*)'PRESSAO DE REFERENCIA:',PRESSAOREF
    ELSE
        WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DE REFERENCIA'
        STOP 200
    END IF
    WRITE(*,*)'########################################'
    !
    keyword_name='cota_de_referencia'
    call readRealKeywordValue(keyword_name, COTAREF, COTAREF, ierr)
    IF(ierr.ne.1) THEN
        WRITE(*,*)'COTA DE REFERENCIA:',COTAREF
    ELSE
        WRITE(*,*) 'ERRO NA LEITURA DA COTA DE REFERENCIA'
        STOP 205
    END IF
    !
    WRITE(*,*)'########################################'
    !
    keyword_name='pressao_no_fundo_do_poco_de_producao'
    call readRealKeywordValue(keyword_name, PRESPROD, PRESPROD, ierr)
    IF(ierr.ne.1) THEN
        WRITE(*,*)'PRESSAO NO POCO DE PRODUCAO:',PRESPROD
    ELSE
        WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DO POCCO'
        STOP 210
    END IF
    WRITE(*,*)'########################################'
    !
    keyword_name='cota_do_fundo_do_poco'
    call readRealKeywordValue(keyword_name, ZFUNDOPOCCO, ZFUNDOPOCCO, ierr)
    IF(ierr.ne.1) THEN
        WRITE(*,*)'COTA DO FUNDO DO POCCO:',ZFUNDOPOCCO
    ELSE
        WRITE(*,*) 'ERRO NA LEITURA DA COTA DO FUNDO DO POCCO'
        STOP 215
    END IF
    WRITE(*,*)'########################################'
    !
    keyword_name='tempo_para_a_descompressao_do_poco'
    call readRealKeywordValue(keyword_name, TPRESPRODWELL, TPRESPRODWELL, ierr)
    IF(ierr.ne.1) THEN
        WRITE(*,*)'TEMPO PARA DESCOMPRESSAO DO POCCO:',TPRESPRODWELL
    ELSE
        WRITE(*,*) 'ERRO NA LEITURA DO TEMPO PARA A DESCOMPRESSAO'
        STOP 220
    END IF
    WRITE(*,*)'########################################'
    !
    RETURN
    !
    END SUBROUTINE

    END MODULE mHidrodinamicaRT
    
    
    
    module mHidrodinamicaGalerkin
    
    implicit none
    public montarEstruturasDadosPressaoSkyline, montarSistemaEquacoesPressao, solveGalerkinPressao
    
    !module variables
    integer :: hgNumPassosTempo
    real*8 :: hgTempoTotal

    !various integers
    integer :: hgNeq
    integer :: hgNdof, hgNlvect
    
    !global arrays
    real*8, allocatable :: hgF(:,:)
    
    !flow variables
    real*8 :: hgInitialPressure
    real*8, allocatable :: hgPrevPressure(:,:), hgPressure(:,:)
    
    !data processing arrays
    integer, allocatable :: hgId(:,:), hgLm(:,:,:)
    
    !skyline data
    integer :: hgNalhs
    integer, allocatable :: hgIdiag(:)
    real*8, allocatable :: hgAlhs(:), hgBrhs(:)
    
    
    contains
    
    !subroutines
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine montarEstruturasDadosPressaoSkyline(conecNodaisElem, nen, nel)
    !function imports
    use mMalha, only:formlm
    
    !variable input
    integer :: nen, nel, conecNodaisElem(nen, nel)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    !generation of lm array
    allocate(hgLm(hgNdof, nen, nel))
    hgLm = 0.0
    call formlm(hgId, conecNodaisElem, hgLm, hgNdof, hgNdof ,nen, nel)
    
    !compute column heights in global left-hand-side matrix
    allocate(hgIdiag(hgNeq))
    hgIdiag = 0.0
    call colht(hgIdiag, hgLm, hgNdof, nen, nel, hgNeq)
    
    !compute diagonal addresses of left-hand-side matrix
    call diag(hgIdiag, hgNeq, hgNalhs)
    
    !initialize alhs
    allocate(hgAlhs(hgNalhs))
    hgAlhs = 0.0
    
    !initialize brhs
    allocate(hgBrhs(hgNeq))
    
    end subroutine montarEstruturasDadosPressaoSkyline
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine montarSistemaEquacoesPressao(curTimeStep, conecNodaisElem, nen, nel, nnp, nsd, x)
    
    !variable input
    integer :: curTimeStep, conecNodaisElem(nen,nel), nen, nel, nnp, nsd
    real*8 :: x(nsd,nnp)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    hgBrhs = 0.0
    if (curTimeStep == 1) then
        if (hgNlvect >= 1) then
            call load(hgId, hgF, hgBrhs, hgNdof, nnp, hgNlvect)
            call dirichletConditions(hgId, hgPressure, hgF, hgNdof, nnp, hgNlvect)
        endif
    endif
    
    call calcCoeficientesSistemaEquacoesPressao(curTimeStep, conecNodaisElem, nen, nel, nnp, nsd, x)
    
    end subroutine montarSistemaEquacoesPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine calcCoeficientesSistemaEquacoesPressao(curTimeStep, conecNodaisElem, nen, nel, nnp, nsd, x)
    !variable imports
    use mGlobaisEscalares, only: nrowsh, npint, betaCompressibility
    use mGlobaisArranjos,  only: mat, c, grav
    
    !function imports
    use mFuncoesDeForma, only: shlt, shlq, shlq3d
    use mFuncoesDeForma, only: shgq, shg3d
    use mSolverGaussSkyline, only: addlhs, addrhs
    use mMalha, only: local
    
    !variable input
    integer :: curTimeStep, conecNodaisElem(nen,nel), nen, nel, nnp, nsd
    real*8 :: x(nsd, nnp)
    
    !variables
    real*8 :: pi !pi constant
    
    integer :: nee ! number of element equations
    real*8 :: elementK(nen, hgNdof * nen), elementF(hgNdof * nen) !element K matrix, element f matrix
    integer :: m !material index
    
    real*8 :: shG(nrowsh,nen,npint), shL(nrowsh,nen,npint) !global shape function with derivative, local shape function with derivative
    real*8 :: det(npint), w(npint) ! jacobian determinant, gauss integration weight
    
    real*8 :: xL(nsd, nnp), pressureL(1, nen), prevPressureL(1, nen) !local position, local pressure, previous local pressure
    real*8 :: kX, kY !x component of permeability, y component of permeability
    
    integer :: ni, nj !node position considering multiple degrees of freedom
    real*8 :: pss
    real*8 :: uup
    real*8 :: djx, djy, djn, dix, diy, din
    
    integer :: curElement, l, i, j ! iterators
    real*8 :: temp1, gf1, gf2, gf3 !temporary variables
    
    logical :: quad, zerodl !is diagonal, is degenerated triange, is zero
    
    real*8 :: deltaT
    
    !------------------------------------------------------------------------------------------------------------------------------------
    if (nsd > 2) then
        stop 9
    end if
    
    ! set some constants
    nee = hgNdof * nen
    pi = 4.d00*datan(1.d00)
    deltaT = hgTempoTotal / hgNumPassosTempo
    
    w = 0.0
    shL = 0.0

    ! construct shape functions and calculates quadrature weights
    if(nen==3) call shlt(shL,w,npint,nen)
    if(nen==4) call shlq(shL,w,npint,nen)
    if(nen==8) call shlq3d(shL,w,npint,nen)
    
    do curElement = 1, nel ! foreach element
        !clear stiffness matrix and force array
        elementK=0.0
        elementF=0.0
        
        call local(conecNodaisElem(1,curElement), x, xL, nen, nsd, nsd) !localize coordinates
        call local(conecNodaisElem(1,curElement), hgPressure, pressureL, nen, hgNdof, hgNdof) !localize dirichlet b.c.
        if (hgNumPassosTempo > 1) then
            call local(conecNodaisElem(1,curElement), hgPrevPressure, prevPressureL, nen, hgNdof, hgNdof) !localize dirichlet b.c.
        end if
        
        ! set material index
        m = mat(curElement)
        
        ! check if element has any coalesced nodes
        quad = .true.
        if (nen.eq.4.and.conecNodaisElem(3,curElement).eq.conecNodaisElem(4,curElement)) quad = .false.
        
        ! calculates global derivatives of shape functions and jacobian determinants
        if(nen==3) call shgq  (xl,det,shL,shG,npint,nel,quad,nen)
        if(nen==4) call shgq  (xl,det,shL,shG,npint,nel,quad,nen)
        if(nen==8) call shg3d (xl,det,shL,shG,npint,nel,nen)
        
        ! retrieve permeabilities from material
        kX = c(1,m)
        kY = c(2,m)
        
        do l = 1, npint ! foreach integration point
            temp1 = w(l)*det(l)
            
            gf1 = grav(1)
            gf2 = grav(2)
            
            pss = 0.0
            if (hgNumPassosTempo > 1) then
                uup = 0.0
                do j = 1, nen ! foreach node on element
                    uup = uup + shg(3, j, l) * prevPressureL(1, j)
                end do
            end if
            
            do j = 1, nen !foreach node on element
                nj = hgNdof * j
                djx = shg(1,j,l) * temp1
                djy = shg(2,j,l) * temp1
                
                djn = shg(nrowsh,j,l) * temp1
                
                ! source terms
                if(hgNumPassosTempo > 1) elementF(nj) = elementF(nj) + (betaCompressibility * djn * uup) 
            end do
        end do

        do i = 1, nen
            do j = 1, nen !for each pair of nodes in the element
                nj = hgNdof * j
                djx = shg(1, j, l) * temp1
                djy = shg(2, j, l) * temp1
                djn = shg(nrowsh, j, l) * temp1
                
                ni = hgNdof * i
                dix = shg(1, i, l)
                diy = shg(2, i, l)
                din = shg(3, i, l)

                if(hgNumPassosTempo > 1) then
                    elementK(nj, ni) = elementK(nj, ni) + deltaT * Kx * dix * djx
                    elementK(nj, ni) = elementK(nj, ni) + deltaT * Ky * diy * djy
                    elementK(nj, ni) = elementK(nj, ni) + betaCompressibility * din * djn
                else
                    elementK(nj, ni) = elementK(nj, ni) + Kx * dix * djx
                    elementK(nj, ni) = elementK(nj, ni) + Ky * diy * djy
                endif
            end do
        end do
        
        ! computate dirichlet bc contribution to the system
        call zTest(pressureL, nee, zerodl)
        if(.not.zerodl) then
            call kdbc(elementK, elementF, pressureL, nee)
        endif
        
        ! assemble element stifness matrix and force array into the global matrixes
        if(curTimeStep == 1) call addlhs (hgAlhs, elementK, hgLm, hgIdiag, nee, .false., .true.)
        call addrhs(hgBrhs, elementF, hgLm, nee)
    end do 
    
    end subroutine calcCoeficientesSistemaEquacoesPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine solveGalerkinPressao(curTimeStep, nnp)
    !function imports
    use mSolverGaussSkyline, only: factor, back
    !variables input
    integer :: curTimeStep, nnp
    
    !------------------------------------------------------------------------------------------------------------------------------------
    !solve by LU decomposition
    if (curTimeStep == 1) call factor(hgAlhs, hgIdiag, hgNeq)
    call back(hgAlhs, hgBrhs, hgIdiag, hgNeq)
    
    !store the result
    call btod(hgId, hgPressure, hgBrhs, hgNdof, nnp)
    
    end subroutine solveGalerkinPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    end module mHidrodinamicaGalerkin