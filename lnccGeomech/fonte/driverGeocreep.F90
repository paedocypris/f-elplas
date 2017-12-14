    !
    !         programa de elementos finitos em fortran 90
    !         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
    !
    !         Eduardo Garcia e Tuane Lopes
    !         bidu@lncc.br, tuane@lncc.br
    !
    !         LNCC/MCT
    !         Petropolis, 07.2013
    !
    !     ************************************************************
    !     *                                                          *
    !     *                                                          *
    !     *         A LINEAR STATIC FINITE ELEMENT PROGRAM FOR       *
    !     *                                                          *
    !     *                 GALERKIN METHOD                          *
    !     *                                                          *
    !     *                                                          *
    !     ************************************************************
    !
    program reservoirSimulator
    !
    use mLeituraEscrita,   only: fecharArquivosBase
    use mLeituraEscritaSimHidroGeoMec,   only: fecharArquivosSimHidroGeoMec
    !
    implicit none

    !
    !.... initialization phase
    !
    print*, ""
    print*, "PREPROCESSAMENTO"
    call preprocessador_DS()
    !
    !.... solution phase
    !
    print*, ""
    print*, "Iniciando o PROCESSAMENTO..."

    !processa o escoamento
    call processamentoOneWayPlast()
    call processamentoTwoWayPlast(1)
    call processamentoTwoWayPlast(2)
    call processamentoTwoWayElast()

    !
    call fecharArquivosBase()
    call fecharArquivosSimHidroGeoMec()

    !
    end program reservoirSimulator
    !
    !**** NEW ********************************************************************
    !
    subroutine preprocessador_DS()
    !
    use mLeituraEscritaSimHidroGeoMec, only : CYLINDER, SOLIDONLY
    use mGlobaisArranjos,  only: listaSolverDisponivel
    use mGlobaisEscalares
    use mGeomecanica, only: ndofD, nlvectD
    !
    use mGeomecanica,      only: idDesloc, idiagD, neqD, nalhsD
    !
    use mMalha,            only: nen, nsd
    use mMalha,            only: numLadosElem, numLadosReserv
    use mMalha,            only: numnp, x, numnpReserv

    use mMalha,            only: IrregMesh
    use mMalha,            only: I4SeaLoad
    use mMalha,            only: FNCPROCESS, FNCINIT, FNCSOLID
    !
    use mLeituraEscrita,                 only: iecho, echo, iflag_tipoPrint
    use mLeituraEscritaSimHidroGeoMec,   only: leituraGeoformations_DS, leituraGeoformations3D_DS, ifdata
    use mLeituraEscritaSimHidroGeoMec,   only: inittime,lerDataIn_DS,lerNumericParam_DS
    use mLeituraEscritaSimHidroGeoMec,   only: lerCreepParam_DS, lerGeoMechParam_DS
    use mLeituraEscritaSimHidroGeoMec,   only: lerGeoMechParam3D
    use mLeituraEscritaSimHidroGeoMec,   only: SETUPDX, abrirArquivosResultados, lerSimulatorParam_DS, lerRandfilesIn_DS
    use mLeituraEscritaSimHidroGeoMec,   only: readSetupPhaseDS
    !
    use mPropGeoFisica,    only: nr
    use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv
    use mPropGeoFisica,    only: XTERLOAD, GEOMECLAW
    !
    use mGeomecanica,      only: fDesloc, InSeaLoad, optSolverD
    use mGeomecanica, only: initStress
    use mInputReader,      only: readInputFileDS, LEIturaGERacaoCOORdenadasDS, readIntegerKeywordValue
    use mInputReader,      only: leituraCodigosCondContornoDS,leituraValoresCondContornoDS

    !
    implicit none

    character(len=50) :: keyword_name
    integer :: ierr
    !
    tempoTotalVelocidade  = 0.d00
    tempoTotalPressao     = 0.d00
    tempoTotalTransporte  = 0.d00
    tempoTotalGeomecanica = 0.d00
    tempoMontagemGeo=0.d00; tempoSolverGeo=0.d00
    tempoMontagemVel=0.d00; tempoSolverVel=0.d00
    !
    !.... Tipo Malha
    !
    novaMalha = .false.
    !
    !.... input phase
    !
    INITS3    = .TRUE.
    CYLINDER  = .FALSE.
    SALTCREEP = .FALSE.
    SOLIDONLY = .FALSE.
    LDRAINED  = .FALSE.
    optSolverD='skyline'

    open(unit=ifdata, file= 'inputDS.dat'  )
    call readInputFileDS(ifdata)

    call readSetupPhaseDS(nlvectD)

    call identificaSolversDisponiveis(listaSolverDisponivel)
    call verificarSolver(optSolverD, listaSolverDisponivel)

    !
    !.... STRESS DIMENSION PHISICS: S3DIM
    !
    S3DIM = FNCPROCESS(TypeProcess)
    !
    INITS3 = FNCINIT(TypeProcess)
    !
    SOLIDONLY = FNCSOLID(TypeProcess)
    !
    IF (SOLIDONLY) INITS3 = .FALSE.
    !
    IF (INITS3) LDRAINED = .FALSE.
    !
    !.... Logical novaMalha
    IF (IrregMesh.eq.1) novaMalha=.true.
    !
    !
    if (nsd==2) then
        numLadosReserv=(nelxReserv+1)*(nelyReserv+1)*2
        numLadosReserv=numLadosReserv - ((nelxReserv+1)+(nelyReserv+1))
    endif
    if (nsd==3) then
        numLadosReserv=(nelxReserv*nelyReserv*nelzReserv)
        numLadosReserv=numLadosReserv+(nelxReserv*nelyReserv)
        numLadosReserv=numLadosReserv+((nelxReserv+1)*nelyReserv*nelzReserv)
        numLadosReserv=numLadosReserv+ (nelyReserv+1)*nelxReserv*nelzReserv
    endif
    !
    print*, "numLadosReserv=", numLadosReserv
    !
    !.... inicializa os parametros da simulacao
    !
    call lerDataIn_DS
    CALL lerGeoMechParam_DS
    !
    call lerNumericParam_DS
    call lerCreepParam_DS
    call lerSimulatorParam_DS
    !
    call abrirArquivosResultados
    call lerRandfilesIn_DS

    if(iflag_tipoPrint.eq.3)  CALL SETUPDX(NUMDX,  NITGEO, NITHIDRO, SALTCREEP)
    !
    IF ((SOLIDONLY).AND.(SALTCREEP)) THEN
        CYLINDER = .TRUE.
        GEOMECLAW(1) = 'CREEP'
    ENDIF
    !
    !.... initialization phase
    !
    numLadosElem = 2*NSD
    nen          = 2**NSD
    ndofD        = NSD
    !
    nr    = 1
    !
    !.... input coordinate data well
    !
    !print*,"leituraCoordenadasPoco"
    !call leituraCoordenadasPoco(NCONDP,PCONDP)
    !
    call alocarMemoria()
    !
    !.... input coordinate data
    !
    print*, "LEITuraGERAcaoCOORdenadas_DS"
    call leituraGeracaoCoordenadasDS(x, nsd, numnp, iprtin, iecho)

    print*, "leituraGeoformations_DS"
    if(nsd==2) call leituraGeoformations_DS  (x, nsd, numnp)
    if(nsd==3) call leituraGeoformations3D_DS(x, nsd, numnp)
    !
    !.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
    !
    keyword_name = "codigos_cond_contorno_desloc"
    call leituraCodigosCondContornoDS(keyword_name,idDesloc,ndofD,numnp,neqD,iecho,iprtin)
    allocate(idiagD(neqD));  idiagD=0
    !
    !.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
    !
    keyword_name = "valores_cond_contorno_desloc"
    if(nlvectD.gt.0) call leituraValoresCondContornoDS(keyword_name, fDesloc, ndofD, numnp, 0, nlvectD, iprtin, iecho)
    
    keyword_name = "initStress"
    call readIntegerKeywordValue(keyword_name, initStress, 1, ierr)
    !
    !.... NEXT MULTIPLY SEALOAD ON Y DIRECTION 2-D MODEL NEUMANN CONDITIONS
    !
    IF (I4SeaLoad.EQ.1) CALL InSeaLoad(FDESLOC,NDOFD,NSD,NUMNP, &
        &                         NLVECTD,XTERLOAD)

    !
    !.... input element data
    !
    call TOPologiaMALhaSistEQUAcoesDS(NALHSD, NEQD)
    !
    !.... inicializa os tempos de impressao
    !
    call inittime
    !
    !    inicializa a pressão para o galerkin
    !
    call initGalerkinPressure(numnpReserv)
    
    !inicializa variáveis de drucker prager
    call initDPProperties
    !
1000 format(20a4)
3000 format(//&
        ' input data print code . . . . . . . . . . .(iprtin) = ',i10//5x,&
        '    eq. 0, print nodal and element input data          ',   /5x,&
        '    eq. 1, do not print nodal and element input data   ',   /5x, &
        ' number of space dimensions  . . . . . . . .(nsd   ) = ',i10)
4000 format(5x,&
        ' number of nodal points  . . . . . . . . .  (numnp ) = ',i10//5x,&
        ' number of Lados         . . . . . . . . .(numLados) = ',i10//5x,&
        ' number of nodal degrees-of-freedom  . . . (ndofP ) = ',i10//5x,&
        ' number of nodal degrees-of-freedom  . . . (ndofV ) = ',i10//5x,&
        ' number of nodal degrees-of-freedom  . . . (ndofD ) = ',i10//5x,&
        ' number of load vectors  . . . . . . . . . (nlvectP) = ',i10//5x,&
        ' number of load vectors  . . . . . . . . . (nlvectV) = ',i10//5x,&
        ' number of load vectors  . . . . . . . . . (nlvectD) = ',i10//5x)
5000 FORMAT(5X,  &
        &' MESH DATA FOR RESERVOIR AND OVERBUDEN DOMAINS:     '//5X,     &
        &'  ELEMENTS IN X-DIRECTION GLOBAL DOMAIN. . (NELX ) = ',I10//5X, &
        &'  ELEMENTS IN Y-DIRECTION GLOBAL DOMAIN. . (NELY ) = ',I10//5X, &
        &'  ELEMENTS IN Z-DIRECTION GLOBAL DOMAIN. . (NELZ ) = ',I10//5X, &
        &'  NODAL POINTS FOR GLOBAL DOMAIN . . . . . (NUMNP) = ',I10//5X, &
        &'  ELEMENTS IN X-DIRECTION RESERVOIR. . . . (NELX1) = ',I10//5X, &
        &'  ELEMENTS IN Y-DIRECTION RESERVOIR. . . . (NELY1) = ',I10//5X, &
        &'  ELEMENTS IN Z-DIRECTION RESERVOIR. . . . (NELZ1) = ',I10//5X, &
        &'  NODAL POINTS FOR RESERVOIR . . . . . . . (NUMNP1)= ',I10//5X)
    !
6000 FORMAT(5X,  &
        &' MESH STRUCTURE AND READ INPUT DATA OPTIONS:           '//5X, &
        &' Reservoir Fine Mesh and External Gross  (IrregMesh) = ',I2/5X,&
        &'    eq. 0, Not Read External Files                     ',  /5x,&
        &'    eq. 1, Read Files with Coordinates and Conectivies ',  //5x,&
        &' Read Geomechanical Dirichlet Conditions (Dirichlet) = ',I2/5X,&
        &'    eq. 0, Not Read External Files                     ',  /5x,&
        &'    eq. 1, Read External File                          ',  //5x,&
        &' Read Geomechanical Neumann   Conditions (Neumann  ) = ',I2/5X,&
        &'    eq. 0, Not Read External Files                     ',  /5x,&
        &'    eq. 1, Read External File                          ',  //5x,&
        &' Multiply Sea Load on Neumann Conditions (I4SeaLoad) = ',I2/5X,&
        &'    eq. 0, Not Multiply by External Load               ',   /5x,&
        &'    eq. 1, Multiply by External Load                   ',  //)
9001 FORMAT("Solvers escolhidos para a solucao do sistema de equacoes: Hidro:",A7, ", Geo:", A7)
    !
    end subroutine preprocessador_DS
    !
    !**** new *******************************************************************
    !
    subroutine verificarSolver(optSolver_, listaSolverDisponivel_)
    !
    character(len=7), intent(in) :: optSolver_
    logical, intent(in) :: listaSolverDisponivel_(*)
    !
    if(optSolver_=='pardiso')then
        if (listaSolverDisponivel_(2).eqv..false.) then
            print*, "O Solver escolhido ...,  ", optSolver_,",  não está disponível"
            stop 100
        endif
    endif

    if(optSolver_=='hypre') then
        if(listaSolverDisponivel_(3).eqv..false.) then
            print*, "O Solver escolhido ...,  ", optSolver_, ",  não está disponível"
            stop 100
        endif
    endif

    end subroutine verificarSolver
    !
    !**** new *******************************************************************
    !
    subroutine identificaSolversDisponiveis(listaSolverDisponivel_)

    logical, intent(inout) :: listaSolverDisponivel_(*)

    listaSolverDisponivel_(1)=.true. !skyline
    listaSolverDisponivel_(2)=.false. !pardiso
    listaSolverDisponivel_(3)=.false. !hypre

    print*, "Solvers disponiveis:"
    print*, "                      SKYLINE"
    if(listaSolverDisponivel_(2).eqv..true.) print*, "                      PARDISO"
    if(listaSolverDisponivel_(3).eqv..true.) print*, "                      HYPRE"

    end subroutine
    !
    !**** new ******************************************************************
    !
    SUBROUTINE UPDTINIT(DIS,DIS0,NDOFD,NUMNP,LFLAG)
    !
    IMPLICIT NONE
    !
    LOGICAL :: LFLAG
    INTEGER :: I,NODE, NDOFD, NUMNP
    REAL(8), DIMENSION(NDOFD,NUMNP) :: DIS, DIS0
    !
    IF (.NOT.LFLAG) RETURN

    DO 200 NODE=1, NUMNP
        DO 100 I=1,NDOFD
            DIS(I,NODE) = DIS(I,NODE)-DIS0(I,NODE)
100     CONTINUE
200 CONTINUE
    !
    RETURN
4500 FORMAT(I8,X,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW ******************************************************************
    !
    FUNCTION TIMELOAD(REALTIME)
    !
    use mGlobaisEscalares, only: YEARINJ
    !
    REAL*8  :: X, TIMELOAD, REALTIME, TOL
    !
    TOL = 1e-6
    IF (REALTIME-YEARINJ.LE.TOL) THEN
        X = 0.0D0
    ELSE
        X = 1.0D0
    ENDIF
    !
    TIMELOAD = X
    !
    END FUNCTION
    !
    !**** NEW ******************************************************************
    !
    FUNCTION SIGNORM(A,B,N)
    !
    !.... FUNCTION TO COMPUTE SUPREMUM NORM OF ARRAY DIFFERENCE
    !
    INTEGER :: I, N
    REAL(8) :: SIGNORM, XMAXIMO, XDIFF
    REAL(8), DIMENSION(N) :: A, B
    !
    XMAXIMO = 1.0D0
    !
    DO 100 I=1,N
        XDIFF = DABS(A(I)-B(I))
        XMAXIMO = DMAX1(XDIFF,XMAXIMO)
100 CONTINUE
    !
    SIGNORM = XMAXIMO
    !
    END FUNCTION
    !
    !**** NEW ******************************************************************
    !
    FUNCTION VELNORM(A,B,N)
    !
    !.... FUNCTION TO COMPUTE EUCLIDEAN NORM OF VELOCITIES DIFFERENCE
    !
    INTEGER :: I, N
    REAL(8) :: VELNORM, XMAXIMO, XDIFF
    REAL(8), DIMENSION(1,N) :: A, B
    !
    XMAXIMO = 0.0D0
    !
    DO 100 I=1,N
        XDIFF   = DABS(A(1,I)-B(1,I))
        XMAXIMO = DMAX1(XDIFF,XMAXIMO)
100 CONTINUE
    !
    VELNORM = XMAXIMO
    !
4000 FORMAT(i5,2x,2(1PE15.8,2x))
    END FUNCTION
    !
    !**** NEW **********************************************************************
    !
    subroutine alocarMemoria()
    use mGlobaisArranjos,  only: uTempoN, mat, grav, beta

    !
    use mMalha,            only: nsd, numel, numelReserv, numnp
    use mMalha,            only: nen
    use mMalha,            only: x, xc
    !
    use mMalha,            only: listaDosElemsPorNo
    use mMalha,            only: conecNodaisElem
    use mGeomecanica,      only: EINELAS
    use mGeomecanica,      only: nrowB, nintD, nrowB2
    use mGeomecanica,      only: idDesloc, lmD, ndofD, nlVectD, fDesloc
    use mGeomecanica,      only: idDis, dis, dis0, vdP, dtrl
    use mGeomecanica,      only: divU, divU0, strss, strss0, hmTTG, eCreep
    use mGeomecanica,      only: geoPrsr, avStrs, avCrep, strs3D
    use mGeomecanica,      only: sigmaT, sigma0
    use mPropGeoFisica,    only: GEOFORM, MASCN, MASCN0

    !
    implicit none
    !
    !malha
    !
    allocate(x(nsd,numnp))
    x  = 0.0d0
    allocate(xc(nsd,numel))
    xc = 0.0d0
    allocate(conecNodaisElem(nen,numel))
    conecNodaisElem = 0
    allocate(listaDosElemsPorNo(nen,numnp))
    listaDosElemsPorNo = 0
    !
    !
    !material
    allocate(mat(numel)); mat=0.d0
    !
    !gravidade
    allocate(grav(3));   grav=0.d0
    !
    !tempo
    allocate(uTempoN(numel))
    !
    !.... GEOMECANICS 2-D MODEL ARRAYS
    !
    NROWB = 2*NSD
    NINTD = 2**NSD
    NROWB2 = NROWB*NROWB
    ALLOCATE(idDesloc(ndofD,NUMNP));    idDesloc  = 0
    ALLOCATE(LMD(ndofD,NEN,NUMEL));     LMD       = 0
    !
    IF (nlvectD.ne.0) THEN
        ALLOCATE(fDesloc(ndofD,numnp));   fDesloc = 0.0d0
    ENDIF
    ALLOCATE(IDDIS(NDOFD, NUMNP));        IDDIS   = 0
    ALLOCATE(DIS(NDOFD, NUMNP));          DIS     = 0.0D0
    ALLOCATE(DIS0(NDOFD, NUMNP));         DIS0    = 0.0D0
    ALLOCATE(VDP(NDOFD,NUMNP));           VDP     = 0.0D0
    ALLOCATE(DTRL(NDOFD,NUMNP));          DTRL    = 0.0D0
    ALLOCATE(DIVU  (NUMEL));              DIVU    = 0.0D0
    ALLOCATE(DIVU0 (NUMEL));              DIVU0   = 0.0D0
    !
    ALLOCATE(STRSS (NUMEL,NINTD,NROWB));  STRSS   = 0.0D0
    ALLOCATE(HMTTG (nrowb2,NINTD,numel)); HMTTG   = 0.0D0
    ALLOCATE(ECREEP(NUMEL,NINTD,NROWB));  ECREEP  = 0.0D0
    !
    ALLOCATE(GEOPRSR(NUMEL));             GEOPRSR = 0.0D0
    ALLOCATE(STRSS0(NROWB,NUMEL));        STRSS0  = 0.0D0
    ALLOCATE(AVSTRS(NROWB,NUMEL));        AVSTRS  = 0.0D0
    ALLOCATE(AVCREP(NUMEL,NROWB));        AVCREP  = 0.0D0
    ALLOCATE(STRS3D(NROWB,NUMEL));        STRS3D  = 0.0D0
    !
    ALLOCATE(SIGMAT(numelReserv));        SIGMAT  = 0.0D0
    ALLOCATE(SIGMA0(numelReserv));        SIGMA0  = 0.0D0
    ALLOCATE(GEOFORM(NUMEL));             GEOFORM = 'NONE'
    !
    !... ALLOCATE MEMORY 4 MASS CONTENT
    !
    allocate(MASCN (numelReserv));        MASCN   = 0.0D0
    allocate(MASCN0(numelReserv));        MASCN0  = 0.0D0
    !
    !... allocate memory for plasticity
    allocate(EINELAS(NUMEL,NINTD,NROWB)); EINELAS = 0.0d0
    !
    !... END GEOMECANICA
    !
    allocate(beta(numel))
    !
    end subroutine

    !
    !**** NEW **** MODIFIED 4 HIERARCH MESH  ************************************
    !
    SUBroutine TOPologiaMALhaSistEQUAcoesDS(NALHSD,NEQD)
    use mGlobaisEscalares
    use mGeomecanica, only: ndofD, optSolverD
    use mGlobaisArranjos
    use mLeituraEscrita, only: nprint, prntel
    use mLeituraEscritaSimHidroGeoMec, only: GEOREGION_DS
    !
    use mGeomecanica,     only: LMD, IDIAGD, IDDESLOC
    use mGeomecanica,     only: AUXM,IOPT,NED,NED2,NEE,NEE2,NEESQ,NEESQ2,NESD,NSTR
    !
    use mMalha,          only: numel, numnp, nsd, nen
    use mMalha,          only: numnpReserv
    use mMalha,          only: x, xc, local
    use mMalha,          only: conecNodaisElem
    use mMalha,          only: listaDosElemsPorNo
    use mMalha,          only: criarListaVizinhos
    use mMalha,          only: genel
    use mMalha,          only: formlm
    !
    use mPropGeoFisica,  only: hx,hy,hz, calcdim

    use mInputReader,      only: readNodeElementsDS, readMaterialPropertiesDS, readConstantBodyForcesDS
    use mInputReader,      only: genelFacesDS, leituraGeracaoConectividadesDS
    !

    !
    implicit none
    !
    !.... program to set arrays storage
    !
    integer :: NALHSD, NEQD
    integer :: i
    real*8, allocatable :: xl(:,:)

    character(len=50) keyword_name
    integer :: ierr
    !
    !.... calculate hx, hy, hz
    !
    call calcdim(nsd,numnpReserv,x)
    WRITE(*,1001) hx, hy, hz
    !
    !.... set element parameters
    !
    allocate(npar(numParElem))
    call readNodeElementsDS
    !
    NROWSH = NSD + 1
    !
    nprint = 0
    !
    NED    = 1
    NEE    = NEN*NED
    NEESQ  = NEE*NEE
    NNP    = 0
    !..Bbar.. BEGIN
    !..Bbar..  SET ELEMENT PARAMETERS VISCOELASTIC INCOMPRESSIVEL MODEL
    !..Bbar..
    NESD   = NDOFD
    NED2   = NDOFD
    NEE2   = NEN*NED2
    NEESQ2 = NEE2*NEE2
    NSTR   = 3
    IOPT   = 1
    !
    allocate(AUXM(NUMEL,NEE2,NEE2)); AUXM = 0.0D0
    !
    !..Bbar.. END
    !
    !
    !....... set memory pointers
    !
    !     constant body forces
    !
    keyword_name = "constant_body_forces"
    call readConstantBodyForcesDS(keyword_name, ierr)
    !
    !    generation of conectivities
    !
    keyword_name = "conectividades_nodais"
    call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen, ierr)

    listaDosElemsPorNo=0
    call criarListaVizinhos(nen,numnp,numel,conecNodaisElem, &
        &     listaDosElemsPorNo  )

    !     calculo de xc (centros dos elementos)
    allocate(xl(nsd,nen)); xl=0.0
    do i=1,numel
        call local(conecNodaisElem(1,i),x,xl,nen,nsd,nsd)
        xc(1,i) = sum(xl(1,1:nen))/nen
        xc(2,i) = sum(xl(2,1:nen))/nen
        if(nsd==3)xc(3,i) = sum(xl(3,1:nen))/nen
    end do

    !
    !..Bbar.. BEGIN
    !-----------------------------------------------------------------------
    !     B-BAR FORMULATION FOR DISPLACEMENTS
    !-----------------------------------------------------------------------
    !
    !.... CLEAR IDIAG ARRAY
    !
    IDIAGD=0
    !
    !.... GENERATION OF LMD ARRAY
    !
    CALL FORMLM(idDesloc,conecNodaisElem,LMD,NDOFD,NED2,NEN,NUMEL)
    !
    !.... MODIFICATION OF IDIAGD ARRAY
    !
    CALL COLHT(IDIAGD,LMD,NED2,NEN,NUMEL,NEQD)
    !
    !.... DETERMINE ADDRESSES OF DIAGONALS IN LEFT-HAND-SIDE MATRIX
    !
    if(optSolverD=='skyline') then
        CALL DIAG(IDIAGD,NEQD,NALHSD)
    endif
    !
    !.... SETUP GEO-MECHANICAL REGIONS
    !
    CALL GEOREGION_DS(XC,NSD,NUMEL)
    !
    if(optSolverD=='pardiso') then
        write(*,*) ' não implementado '
        stop 9
    endif
    !
    return
    !
1000 format(//,&
        ' two/three-n o d e    e l e m e n t s ',//,5x,&
        ' element type number . . . . . . . . . . (ntype ) = ',i10,//5x,&
        ' number of elements  . . . . . . . . . . (numel ) = ',i10,//5x,&
        ' number of element material sets . . . . (numat ) = ',i10,//5x,&
        ' number of element nodes . . . . . . . . (nen   ) = ',i10,//5x,&
        ' number of integration points. . . . . . (npint  ) = ',i10)
1001 FORMAT(' hx= ',F12.5,2X,'hy= ',F12.5,2X,' hz= ',F12.5)
1500 FORMAT(I10)
4000 format(///,&
        ' m a t e r i a l   s e t   d a t a             ',  //5x,&
        ' number of material sets . . . . . . .(numat ) = ',i10,//,2x,&
        & 'set',4x,'Kx ',4x,'Ky',4x,'Kz')
4500 FORMAT(27I8)
5000 format(i10,5x,5f10.0)
6000 format(2x,i3,1x,5(1x,1pe11.4))
7000 format(8f10.0)
8000 format(///,&
        ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
        ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
        ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
        ' exemplo 3............................ = ',      1pe15.8,//)
9000 format('Calculo das velocidades nodais' //&
        ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
        ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
        ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
        ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
        ' memoria necessaria para a matriz do sistema  (Mbytes)  = ',e10.2)
9002 FORMAT( "Tempo de pre-processamento da Velocidade para solver externo",f12.5)
9003 FORMAT( "Tempo de pre-processamento da Geomecanic para solver externo",f12.5)
9500 FORMAT(5X,  &
        &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
        &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
        &' **                                                  **',/5X,&
        &' **                                                  **',/5X,&
        &' **  CONECTIVITIES STRUCTURE:                        **',/5X,&
        &' **  OBTAINED FROM ', A30 ,                     '    **',/5X,&
        &' **                                                  **',/5X,&
        &' **                                                  **',/5X,&
        &' **                                                  **',/5X,&
        &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
        &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X)
    !
    END  SUBROUTINE
    !
    !**************** ********* ************ *********** ********************
    !
    FUNCTION DESPRESSURIZAR(NPWELL,NNP,PREINIC,AUX,BWP)
    !
    IMPLICIT NONE
    INTEGER :: NPWELL,NNP
    REAL(8) :: AUX,PREINIC,BWP,DESPRESSURIZAR
    IF(NPWELL.GE.NNP)THEN
        BWP = PREINIC-REAL(NNP)*AUX
        write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(*,*)'PRESSAO NO FUNDO DO POCO.....:',BWP
        WRITE(*,*)'PRESSAO MINIMA INICIAL.......:',PREINIC
        write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    END IF
    DESPRESSURIZAR = BWP
    END FUNCTION DESPRESSURIZAR
    

    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine initGalerkinPressure(nnp)
    !variable imports
    use mHidrodinamicaGalerkin, only: hgNSteps, hgDts, hgNTimeSteps
    
    use mHidrodinamicaGalerkin, only: hgInitialPressure
    use mHidrodinamicaGalerkin, only: hgNdof, hgNlvect, hgNeq
    use mHidrodinamicaGalerkin, only: hgF
    use mHidrodinamicaGalerkin, only: hgId

    use mLeituraEscrita, only:iecho
    use mGlobaisEscalares, only: iprtin

    !function imports
    use mInputReader,      only: leituraCodigosCondContornoDS, leituraValoresCondContornoDS
    use mInputReader,      only: readRealKeywordValue, readIntegerKeywordValue
    use mInputReader, only: leituraTabelaTempos

    implicit none

    !variable input
    integer nnp

    !variables
    character(len=50) :: keyword_name
    real*8 :: tempoTotal
    integer :: hgIsTable
    integer :: ierr

    !------------------------------------------------------------------------------------------------------------------------------------
    keyword_name = "hgTimeIsTable"
    call readIntegerKeywordValue(keyword_name, hgIsTable, 0, ierr)
    
    if (hgIsTable == 0) then
        allocate(hgNSteps(1))
        allocate(hgDts(1))
        
        hgNTimeSteps = 1
        
        keyword_name = "numeroPassosTempoP"
        call readIntegerKeywordValue(keyword_name, hgNSteps(1), 0_4, ierr)

        keyword_name = "tempoTotal"
        call readRealKeywordValue(keyword_name, tempoTotal, 0.0d0, ierr)
        tempoTotal = tempoTotal * 2592000.0 ! converts months to seconds
        
        hgDts(1) = tempoTotal / hgNSteps(1)
    else if (hgIsTable == 1) then
        keyword_name = "tabelaTempos"
        call leituraTabelaTempos(keyword_name, hgNSteps, hgDts, hgNTimeSteps)
    end if

    keyword_name = "hgNlvect"
    call readIntegerKeywordValue(keyword_name, hgNlvect, 0_4, ierr)

    keyword_name = "pressao_inicial_fluxo"
    call readRealKeywordValue(keyword_name, hgInitialPressure, 0.0d0, ierr)


    hgNdof = 1

    allocate(hgId(hgNdof,nnp))
    hgId = 0

    if (hgNlvect.ne.0)  then
        allocate(hgF(hgNdof,nnp))
        hgF = 0.0d0
    endif

    keyword_name = "codigos_cond_contorno_fluxo"
    call leituraCodigosCondContornoDS(keyword_name,hgId,hgNdof,nnp,hgNeq,iecho,iprtin)

    keyword_name = "valores_cond_contorno_fluxo"
    call leituraValoresCondContornoDS(keyword_name,hgF, hgNdof,nnp, 1, hgNlvect, iprtin, iecho)

    end subroutine initGalerkinPressure
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine initDPProperties()
    !function imports
    
    !variables import
    use mPropGeoFisica,    only: mcFriction, mcC, dpH, dpAlpha
    
    implicit none
    !variables input
    
    !variables
    integer :: i 
    real*8, parameter :: piConst = 3.141592653589793
    real*8 :: constPlaneStrain
    real*8 :: phiRad
    
    !------------------------------------------------------------------------------------------------------------------------------------
    do i = 1,9
        phiRad = mcFriction(i) * piConst / 180
        
        constPlaneStrain = 3.0d0 / sqrt(9 + 12*(tan(phiRad)**2))
        dpAlpha(i) = constPlaneStrain * tan(phiRad)
        dpH(i) = mcC(i) / tan(phiRad)
    end do
    
    
    end subroutine initDPProperties
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine processamento(way, isPlast, plastType)
    !variable imports
    !flow
    use mHidrodinamicaGalerkin, only:hgInitialPressure
    use mHidrodinamicaGalerkin, only:hgNSteps, hgDts, hgNTimeSteps
    use mHidrodinamicaGalerkin, only:hgNdof
    use mGeomecanica, only:hmTTG, initStress
    use mMalha, only: numnp, numel, nen, nsd, conecNodaisElem
    use mMalha, only: x
    
    !mechanics
    use mGeomecanica, only: ndofD, nrowB, nintD
    
    !function imports
    use mHidrodinamicaGalerkin, only:montarEstruturasDadosPressaoSkyline
    use mHidrodinamicaGalerkin, only: incrementFlowPressureSolutionOneWay, incrementFlowPressureSolutionTwoWay
    use mGeomecanica, only: incrementMechanicPlasticSolution, montarEstruturasDadosDeslocamentoSkyline
    use mGeomecanica, only: incrementMechanicElasticSolution
    use mPropGeoFisica, only: lerPropriedadesFisicas
    
    !variables input
    implicit none
    integer :: way, isPlast, plastType
    
    !variables
    integer :: curTimeStep
    real*8 :: t, deltaT
    character(30) :: filename
    
    integer :: k, i
    
    real*8, allocatable :: u(:,:), uInit(:,:), uDif(:,:), prevUDifIt(:,:)
    real*8, allocatable :: p(:,:), prevP(:,:), prevPIt(:,:)
    real*8, allocatable :: vDarcy(:,:), vDarcyNodal(:,:)
    real*8, allocatable :: stress(:,:,:), out2waySource(:,:)
    real*8, allocatable :: stressS(:,:), prevStressS(:,:)
    real*8, allocatable :: trStrainP(:,:), prevTrStrainP(:,:)
    real*8, allocatable :: stressTotal(:,:,:)
    real*8, allocatable :: strainP(:,:,:)
    
    real*8,allocatable :: elementIsPlast(:)
    
    real*8 :: uNorm, pNorm
    logical :: converged
    
    integer :: outFileUnit
    integer :: currentTimeStepPrint
    
    real*8, external :: matrixNorm

    !------------------------------------------------------------------------------------------------------------------------------------
    !allocate memory
    allocate(u(ndofD, numnp))
    allocate(uInit(ndofD, numnp))
    allocate(uDif(ndofD, numnp))
    allocate(prevUDifIt(ndofD, numnp))
    allocate(p(hgNdof,numnp))
    allocate(prevP(hgNdof,numnp))
    allocate(prevPIt(hgNdof, numnp))
    allocate(vDarcy(nsd,numel))
    allocate(vDarcyNodal(nsd,numnp))
    allocate(stress(nrowb,nintD, numel))
    allocate(stressS(nintD, numel))
    allocate(out2waySource(nintD, numel))
    allocate(prevStressS(nintD, numel))
    allocate(trStrainP(nintD, numel))
    allocate(prevTrStrainP(nintD, numel))
    allocate(stressTotal(nrowb,nintD,numel))
    allocate(strainP(nrowb,nintD,numel))
    allocate(elementIsPlast(numel))
    
    outFileUnit = 5513
    
    !read physical properties
    call lerPropriedadesFisicas()
    
    ! mount data structure
    call montarEstruturasDadosPressaoSkyline(conecNodaisElem, nen, numel)
    call montarEstruturasDadosDeslocamentoSkyline(conecNodaisElem, nen, numel)
    
    !init pressure
    if(hgNTimeSteps > 1 .or. hgNSteps(1) > 1) then
        p = hgInitialPressure
        prevP = p
    endif
    
    !init stress and strain
    write(*,*) "Inicializa estado de tensão"
    u = 0.d0
    strainP = 0.d0
    stress = 0.d0
    elementIsPlast = 0.d0
    if (initStress == 1) then
        if (isPlast == 1) then
            call incrementMechanicPlasticSolution(conecNodaisElem, nen, numel, numnp, nsd, x, u, strainP, stress, stressS, trStrainP, stressTotal, p, elementIsPlast, 0, 1)
        else if (isPlast == 0) then
            call incrementMechanicElasticSolution(conecNodaisElem, nen, numel, numnp, nsd, x, u, stress, stressS, stressTotal, p, 0, 1)
        end if
    end if 
    uInit = u
    uDif = 0.d0
    prevTrStrainP = trStrainP

    ! print initial state
    if (isPlast == 1) then
        if (way==1) filename = "plastSolution1way"
        if (way==2) then
            if (plastType == 1) then
                filename = "plastSolution2waySilva"

                open(unit=outFileUnit, file="out/debugFiles/logSilva.txt", status='replace')
                write (outFileUnit,*) filename

            else if (plastType == 2) then
                filename =  "plastSolution2wayKim"

                open(unit=outFileUnit, file="out/debugFiles/logKim.txt", status='replace')
                write (outFileUnit,*) filename
            end if
        end if
    else if (isPlast == 0) then
        if (way == 1) filename = "elastSolution1way"
        if (way==2) then
            filename = "elastSolution2way"
            
            open(unit=outFileUnit, file="out/debugFiles/logElast.txt", status='replace')
            write (outFileUnit,*) filename
        end if
    end if
    
    call writeCurrentSolution(filename,0, p, u, stress, stressTotal, stressS, out2waySource, vDarcy, vDarcyNodal, elementIsPlast, conecNodaisElem, numnp, numel, nen, nsd, nrowb, nintD)
    
    !time loop
    t = 0
    currentTimeStepPrint = 1

    do i = 1, hgNTimeSteps
        deltaT = hgDts(i)

        do curTimeStep = 1, hgNSteps(i)
            t = t + deltaT
            write(*,*) "tempo", t, "passo de tempo", curTimeStep

            prevStressS = stressS
            prevTrStrainP = trStrainP

            ! begin split loop
            converged = .false.
            elementIsPlast = 0.d0
            do k = 1, 3000
                if (way == 1 .or. k == 1) then
                    call incrementFlowPressureSolutionOneWay(conecNodaisElem, nen, numel, numnp, nsd, x, deltaT, p, prevP, vDarcy, vDarcyNodal)
                    pNorm = 1.0
                    uNorm = 1.0
                else
                    call incrementFlowPressureSolutionTwoWay(conecNodaisElem, nen, numel, numnp, nsd, x, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, vDarcy, vDarcyNodal, hmTTG, out2waySource, plastType)
                end if

                if (isPlast == 1) then
                    call  incrementMechanicPlasticSolution(conecNodaisElem, nen, numel, numnp, nsd, x, u, strainP, stress, stressS, trStrainP, stressTotal, p, elementIsPlast, 0, 1)
                else if (isPlast == 0) then
                    call incrementMechanicElasticSolution(conecNodaisElem, nen, numel, numnp, nsd, x, u, stress, stressS, stressTotal, p, 0, 1)
                end if
                uDif = u - uInit

                if (way == 1) then
                    exit
                end if

                if (k > 1) then
                    prevPIt = prevPIt - p
                    pNorm = matrixNorm(prevPIt, hgNdof, numnp)

                    prevUDifIt = prevUDifIt - uDif
                    uNorm = matrixNorm(prevUDifIt, ndofD, numnp)

                    if (pNorm < 1.0d-6 .and. uNorm < 1.0d-6) then
                        converged = .true.
                        exit
                    end if
                end if

                write(*,*) "k=", k, "pNorm=", pNorm, "uNorm=", uNorm

                prevPIt = p
                prevUDifIt = uDif
            end do

            if ((way == 2) .and.(converged.eqv..false.)) then
                write(*,*) "fixed stress split didn't converged."
                write (outFileUnit,*) "fixed stress split didn't converged"
                exit
            end if

            !update time
            prevP = p

            !print current solution
            call writeCurrentSolution(filename,currentTimeStepPrint, p, uDif, stress, stressTotal, stressS, out2waySource, vDarcy, vDarcyNodal, elementIsPlast, conecNodaisElem, numnp, numel, nen, nsd, nrowb, nintD)

            if (way == 2) write (outFileUnit,*) "load stage = ", i,"curtimestep = ", curTimeStep, "k = ", k
            write(*,*) " "
            
            currentTimeStepPrint = currentTimeStepPrint + 1
        end do
    end do
    
    close(outFileUnit)
    
    deallocate(u)
    deallocate(uInit)
    deallocate(uDif)
    deallocate(prevUDifIt)
    deallocate(p)
    deallocate(prevP)
    deallocate(prevPIt)
    deallocate(vDarcy)
    deallocate(vDarcyNodal)
    deallocate(stress)
    deallocate(stressS)
    deallocate(stressTotal)
    deallocate(prevStressS)
    deallocate(trStrainP)
    deallocate(prevTrStrainP)
    deallocate(strainP)
    deallocate(elementIsPlast)
    
    
    end subroutine processamento
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine processamentoOneWayPlast()
    
    implicit none
    !------------------------------------------------------------------------------------------------------------------------------------
    call processamento(1, 1, 1)
    end subroutine processamentoOneWayPlast
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine processamentoTwoWayElast()
    
    implicit none
    !------------------------------------------------------------------------------------------------------------------------------------
    call processamento(2, 0, 1)
    end subroutine processamentoTwoWayElast
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine processamentoTwoWayPlast(plastType)
    
    implicit none
    integer :: plastType
    
    !------------------------------------------------------------------------------------------------------------------------------------
    call processamento(2, 1, plastType)
    
    end subroutine processamentoTwoWayPlast
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine convertStressStoPrintable(stressS, stressSPrint, npint, nel)
    
    implicit none
    !variables input
    real*8 :: stressS(npint,nel)
    real*8 :: stressSPrint(nel)
    integer :: npint, nel
    
    !variables
    integer :: i,j
    
    !------------------------------------------------------------------------------------------------------------------------------------
    do j = 1, nel
        stressSPrint(j) = 0.
        do i = 1, npint
            stressSPrint(j) = stressSPrint(j) + stressS(i,j)
        end do
        stressSPrint(j) = stressSPrint(j)/npInt
    end do
    
    end subroutine convertStressStoPrintable
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine writeCurrentSolution(filename,curTimeStep, p, u, stressEf, stressTot, stressS, out2waySource, vDarcy, vDarcyNodal, elementIsPlast, conecNodaisElem, nnp, nel, nen, nsd, nrowb, nintD)
    !function imports
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaviewOpening
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaviewIntermed_CampoEscalar
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaviewIntermed_CampoVetorial
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaviewIntermed_CampoTensorialElemento
    
    !variables import
    use mLeituraEscritaSimHidroGeoMec, only:reservHGPres
    
    implicit none
    !variables input
    character(len=*) :: filename
    integer :: curTimeStep
    real*8 :: p(1,nnp), u(nsd,nnp), stressEf(nrowb,nintD,nel), stressTot(nrowb,nintD,nel)
    real*8 :: stressS(nintD,nel), out2waySource(nintD, nel), vDarcy(nsd,nel), vDarcyNodal(nsd,nnp)
    real*8 :: elementIsPlast(nel)
    integer :: conecNodaisElem(nen,nel)
    integer :: nnp, nel, nen, nsd, nrowb, nintD
    
    !variables
    integer :: unitNumber
    real*8, allocatable :: stressSPrint(:)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    unitNumber = 13587
    
    allocate(stressSPrint(nel))
    
    call escreverArqParaviewOpening(filename, p, 1, nnp, nen, conecNodaisElem, 2, 'p', len('p'), reservHGPres, curTimeStep, unitNumber)
    call escreverArqParaviewIntermed_CampoVetorial(u, nsd, nnp, 'u', len('u'), 2, unitNumber)
    call escreverArqParaviewIntermed_CampoVetorial(vDarcy, nsd, nel, 'vDarcy', len('vDarcy'), 1, unitNumber)
    call escreverArqParaviewIntermed_CampoVetorial(vDarcyNodal, nsd, nnp, 'vDarcyNodal', len('vDarcyNodal'), 2, unitNumber)
    call escreverarqparaviewintermed_campotensorialelemento(stressef, nrowb, nintd, nel, 'stressef', len('stressef'), unitnumber)
    call escreverarqparaviewintermed_campotensorialelemento(stresstot, nrowb, nintd, nel, 'stresstot', len('stresstot'), unitnumber)
    call convertStressStoPrintable(stressS, stressSPrint, nintD, nel)
    call escreverarqparaviewintermed_campoescalar(unitnumber, stresssprint, 1, nel, 'stresss', len('stresss'), 1)
    call convertStressStoPrintable(out2waySource, stressSPrint, nintD, nel)
    call escreverarqparaviewintermed_campoescalar(unitnumber, stressSPrint, 1, nel, '2waySource', len('2waySource'), 1)
    call escreverarqparaviewintermed_campoescalar(unitnumber, elementIsPlast, 1, nel, 'elementIsPlast', len('elementIsPlast'), 1)
    close(unitNumber)
    
    deallocate(stressSPrint)
    
    end subroutine writeCurrentSolution
    !************************************************************************************************************************************
    !************************************************************************************************************************************