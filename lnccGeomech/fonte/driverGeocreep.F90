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
    use mLeituraEscrita,   only: fecharArquivosBase, abrirArquivosInformacoesMalha
    use mLeituraEscritaSimHidroGeoMec,   only: fecharArquivosSimHidroGeoMec
    use mgeomecanica,      only: STRESS_INIT
    !
    implicit none
    !
    real*8 :: t1, t2, t3, t4

    call timing(t3)

    !
    !.... initialization phase
    !
    call abrirArquivosInformacoesMalha()
    print*, ""
    print*, "PREPROCESSAMENTO"
    call timing(t1)
    call preprocessador_DS()
    call timing(t2)
    !
    !.... solution phase
    !
    print*, ""
    print*, "Iniciando o PROCESSAMENTO..."

    !processa o escoamento
    call processamentoGalerkinElastico()

    !
    call fecharArquivosBase()
    call fecharArquivosSimHidroGeoMec()

    call timing(t4)
    write(*,*) "TEMPO DE PAREDE TOTAL MEDIDO=", t4 - t3, " segundos "
    !
    end program reservoirSimulator
    !
    !**** NEW ********************************************************************
    !
    subroutine preprocessador_DS()
    !
    use mLeituraEscritaSimHidroGeoMec, only : CYLINDER, SOLIDONLY
    use mGlobaisArranjos,  only: title, listaSolverDisponivel
    use mGlobaisEscalares
    use mGeomecanica, only: ndofD, nlvectD
    use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP
    !
    use mGeomecanica,      only: idDesloc, idiagD, neqD, nalhsD
    use mHidroDinamicaRT,  only: neqV, nalhsV, idiagV, idVeloc,simetriaVel
    !
    use mMalha,            only: nen, nsd
    use mMalha,            only: numLadosElem, numLadosReserv
    use mMalha,            only: numnp, x, numelReserv, numnpReserv

    use mMalha,            only: IrregMesh, Dirichlet, Neumann
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
    use mPropGeoFisica,    only: lerPropriedadesFisicas
    use mPropGeoFisica,    only: nelx, nely, nelz
    use mPropGeoFisica,    only: nr
    use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv
    use mPropGeoFisica,    only: XTERLOAD, GEOMECLAW
    !
    use mHidrodinamicaRT,  only: fVeloc, leituraCoordenadasPoco
    use mHidroDinamicaRT,  only: NCONDP,PCONDP, lerParametrosHidrodinamica_DS, optSolverV
    use mGeomecanica,      only: fDesloc, InSeaLoad, optSolverD, simetriaGeo
    use mTransporte,       only: satElem, satinit
    use mInputReader,      only: readInputFileDS, LEIturaGERacaoCOORdenadasDS
    use mInputReader,      only: leituraCodigosCondContornoDS,leituraValoresCondContornoDS

    !
    implicit none

    character(len=50) :: keyword_name
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
    optSolverV='skyline'
    optSolverD='skyline'

    open(unit=ifdata, file= 'inputDS.dat'  )
    call readInputFileDS(ifdata)

    call readSetupPhaseDS(nlvectV, nlvectP, nlvectD, optSolverV, optSolverD)

    call identificaSolversDisponiveis(listaSolverDisponivel)
    call verificarSolver(optSolverV, listaSolverDisponivel)
    call verificarSolver(optSolverD, listaSolverDisponivel)
    write(*,9001)  optSolverV, optSolverD

    if(optSolverV=="pardiso") simetriaVel=.true.
    if(optSolverD=="pardiso") simetriaGeo=.true.
    if(optSolverV=="hypre")   simetriaVel=.false.
    if(optSolverD=="hypre")   simetriaGeo=.false.
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
    write(iecho,1000) title
    write(iecho,3000) iprtin, nsd
    write(iecho,4000) numnp, numLadosReserv, ndofP, ndofV, ndofD, &
        nlvectP, nlvectV, nlvectD
    WRITE(IECHO,5000) nelx,nely,nelz,numnp, &
        nelxReserv,nelyReserv,nelzReserv,numnpReserv
    WRITE(IECHO,6000) IrregMesh, Dirichlet, Neumann, I4SeaLoad
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
    call lerParametrosHidrodinamica_DS
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
    ndofP = 1
    ndofV = 1
    !
    numLadosElem = 2*NSD
    nen          = 2**NSD
    ndofD        = NSD
    !
    nr    = 1
    !
    !.... input coordinate data well
    !
    print*,"leituraCoordenadasPoco"
    call leituraCoordenadasPoco(NCONDP,PCONDP)
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
    !.... input boundary condition data and establish equation numbers
    !
    if (.not.SOLIDONLY) then
        keyword_name = "codigos_cond_contorno_veloc"
        call leituraCodigosCondContornoDS(keyword_name, idVeloc,ndofV,numLadosReserv,neqV,iecho,iprtin)
        !
        allocate(idiagV(neqV));  idiagV=0
    end if
    !
    !.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
    !
    keyword_name = "codigos_cond_contorno_desloc"
    call leituraCodigosCondContornoDS(keyword_name,idDesloc,ndofD,numnp,neqD,iecho,iprtin)
    allocate(idiagD(neqD));  idiagD=0
    !
    print*, "ndofV=", ndofV, "neqV=", neqV
    print*, "ndofD=", ndofD, "neqD=", neqD
    !
    !.... input nodal force and prescribed kinematic boundary-value data
    !
    if (.not.solidonly) then
        keyword_name = "valores_cond_contorno_veloc"
        if (nlvectV.gt.0) call leituraValoresCondContornoDS(keyword_name,fVeloc,ndofV,numLadosReserv,1,nlvectV,iprtin, iecho)
    end if
    !
    !.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
    !
    keyword_name = "valores_cond_contorno_desloc"
    if(nlvectD.gt.0) call leituraValoresCondContornoDS(keyword_name, fDesloc, ndofD, numnp, 0, nlvectD, iprtin, iecho)
    !
    !.... NEXT MULTIPLY SEALOAD ON Y DIRECTION 2-D MODEL NEUMANN CONDITIONS
    !
    IF (I4SeaLoad.EQ.1) CALL InSeaLoad(FDESLOC,NDOFD,NSD,NUMNP, &
        &                         NLVECTD,XTERLOAD)

    !
    !.... input element data
    !
    call TOPologiaMALhaSistEQUAcoesDS(NALHSV, NEQV, NALHSD, NEQD)
    !
    !.... inicializa os tempos de impressao
    !
    call inittime
    !
    !.... estabelece a condicao inicial para a saturacao
    !
    call satinit(nsd,numelReserv,satElem)
    !
    !    inicializa a pressão para o galerkin
    !
    call initGalerkinPressure(numnpReserv)
    !
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
    !**** new *******************************************************************
    !
    SUBROUTINE CREEP_EXAMPLE()
    !
    use mGlobaisEscalares
    use mLeituraEscrita,   only : iflag_tipoPrint
    use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXMESH
    !
    use mMalha,            only : x, conecNodaisElem
    !
    use mPropGeoFisica,    only : YOUNG, PERMKX
    use mPropGeoFisica,    only : TOLCREEP
    !
    use mGeomecanica,      only : DIS, DTRL, DIS0
    use mGeomecanica,      only : DIVU, STRSS0, GEOPRSR, POS4STRS
    use mGeomecanica,      only : NWTNITER, RESIDUAL, RESMAX
    use mGeomecanica,      only: geomechanic
    !
    implicit none
    !
    LOGICAL  :: LJUMP
    REAL(8)  :: ERRSIZE
    !
    LJUMP = .FALSE.
    !
    !.... MOUNT STIFFNESS MATRIX OF GEOMECHANIC
    !
50  CONTINUE
    !....
    NWTNITER = NWTNITER + 1
    WRITE(*,2500) NWTNITER
    !....
    IF (NWTNITER.GT.MAXITERC) LJUMP = .TRUE.
    !
    CALL GEOMECHANIC('GEOMECHANICS_CREEP')
    ERRSIZE = RESIDUAL/RESMAX
    write(*,5000) ERRSIZE
    ! !
    !.... ..TEST RESIDUAL NORM WITH TOLERANCE CRITERIA 4 CREEP
    !
    IF ((ERRSIZE.GT.TOLCREEP).AND.(.NOT.LJUMP)) GOTO 50

    !....
70  CONTINUE
    !
    !.... UPDATE DISPLACEMENT
    !
    DIS = DTRL
    !
    !.... POST-PROCESS STRESS FIELD
    !
    CALL POS4STRS(X, conecNodaisElem, STRSS0, DIVU)
    !
    !.... MOVE COMPUTED DISPLACEMENTS (DIS) TO INITIAL DISPLACEMENTS (DIS0)
    !
    DIS0 = DIS
    !
    !.... PRINT GEOMECHANICAL INITIAL CONDITIONS
    !
    if(iflag_tipoPrint.eq.3)then
        CALL PRINT_DXMESH(X,DIS0,GEOPRSR,STRSS0,conecNodaisElem, &
            &     YOUNG,PERMKX)
    end if
    !
    DIS  = 0.0D0
    DIVU = 0.0D0
    !
    WRITE(*,*) "  "
    WRITE(*,*) "*** ************* ************* ************* ***"
    WRITE(*,*) "***                                           ***"
    WRITE(*,*) "***    END OF NON-LINEAR CREEP EXAMPLE        ***"
    WRITE(*,*) "***  SEE FILE 'dxcreep01.ht01/cilindr.dat'    ***"
    WRITE(*,*) "***     FOR RADIAL STRESS OUTPUT              ***"
    WRITE(*,*) "***                                           ***"
    WRITE(*,*) "*** ************* ************* ************* ***"
    WRITE(*,*) "  "
    !
    RETURN
    !
2500 FORMAT('    NEWTON ITERATION COUNTER =',I5)
4500 FORMAT(I8,X,40(1PE15.8,2X))
5000 FORMAT('    RESIDUAL/RMAX = ',1PE15.8)
    !
    END SUBROUTINE
    !
    !**** new *******************************************************************
    !
    subroutine processamento_elast()
    !
    use mGlobaisEscalares
    use mLeituraEscrita, only: iflag_tipoPrint
    use mGeomecanica, only: ndofD
    use mHidrodinamicaRT, only: ndofV, ndofP
    !
    use mLeituraEscritaSimHidroGeoMec,   only : imprimirCondicoesIniciais
    use mLeituraEscritaSimHidroGeoMec,   only : imprimirSolucaoNoTempo
    use mLeituraEscritaSimHidroGeoMec,   only : isat,escreverArqParaviewIntermed
    use mLeituraEscritaSimHidroGeoMec,   only : iflag_sat
    use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXINFO, IFEDX
    !
    use mMalha,            only : nsd, numel,numelReserv, nen
    use mMalha,            only : numnp
    use mMalha,            only : numLadosReserv, numLadosElem
    use mMalha,            only : x
    use mMalha,            only : conecNodaisElem, conecLadaisElem
    !
    use mPropGeoFisica,    only : phi
    use mPropGeoFisica,    only : permkx, phi0
    use mPropGeoFisica,    only : lerPropriedadesFisicas
    use mPropGeoFisica,    only : PORE, PORE0, PHIEULER
    use mPropGeoFisica,    only : YOUNG, MASCN0, MASCN
    !
    use mTransporte,       only : transport, calc_prod
    use mTransporte,       only : satElem
    use mTransporte,       only : satElemL, satElemL0, satElem0
    !
    use mGeomecanica,      only : SIGMAT, SIGMA0, GEOTIME
    use mGeomecanica,      only : NROWB
    use mGeomecanica,      only : DIS, DIS0, AVSTRS, PRINT_DX
    use mGeomecanica,      only : DIVU0, DIVU, STRSS0, TMANDEL, GEOPRSR
    use mGeomecanica,      only: geomechanic
    !
    use mHidrodinamicaRT,  only : hidroGeomecanicaRT,  pressaoElem
    use mHidrodinamicaRT,  only : velocLadal, pressaoElemAnt, vc, velocNodal
    use mHidroDinamicaRT,  only : NPRESPRODWELL, PRESMEDIAINICIAL, PRESPROD
    !
    IMPLICIT NONE
    !
    !.... solution driver program
    !
    LOGICAL       :: JUMPINDXK, JUMPINDXL
    INTEGER       :: INDEXL, INDEXK
    character(21) :: labelTransp
    !
    REAL(8)       :: ERRSZSIG, ERRSZVEL, ERRSIG0, ERRVEL0
    REAL(8)       :: SIGNORM, VELNORM
    REAL*8        :: TIMEINJ, TIMELOAD, AUX
    REAL(8), external :: DESPRESSURIZAR
    REAL(8), external :: DESPRESSURIZAR_INIT
    REAL(8), DIMENSION(numelReserv)      :: SIGMAK
    REAL(8), DIMENSION(1,numLadosReserv) :: VELOCITYL
    !
    !.... imprime condicoes iniciais
    !
    call imprimirCondicoesIniciais(GEOPRSR, pressaoElem, velocLadal, velocNodal, vc, phi, &
        permkx, satElem, YOUNG, DIS, STRSS0, MASCN, ndofV, ndofP, ndofD, nrowb)
    !
    IF(iflag_tipoPrint.EQ.3)THEN
        IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN
                CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
                    &                       VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
        ENDIF
    END IF
    !
    !.... tempo de simulacao para cada bloco do transporte
    !
    dtBlocoTransp=tt/DFLOAT(nvel)
    !
    !.... INICIALIZACAO DO TEMPO DA GEOMECANICA
    !
    GEOTIME = TT/DFLOAT(NVEL)
    !
    IF (TypeProcess.EQ.'MANDEL') TMANDEL = 0.0D0
    !
    !-----------------------------------------------------------------------
    !
    !.... DIMINUICAO GRADUAL DA PRESSAO NO POCO
    !
    !-----------------------------------------------------------------------
    AUX = DESPRESSURIZAR_INIT(PRESPROD,dtBlocoTransp,nvel)
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     SIMULACAO
    !
    !-----------------------------------------------------------------------
    !
    DO NNP=1,NVEL

        ! DESPRESURIZACAO LENTA DO POCO DE PRODUCAO
        PRESPROD = DESPRESSURIZAR(NPRESPRODWELL,NNP,PRESMEDIAINICIAL,AUX,PRESPROD)
        !
        write(*,1000) tTransporte,nnp,nvel
        !...  ..SETUP SATURATION FOR ITERATIVE HIDRODINAMICS AND TRANSPORT
        satElem0       = satElem
        satElemL       = satElem
        !...  ..UPDATE FOR MACROTIME EVOLUTION
        pressaoElemAnt = pressaoElem(1,:)
        SIGMA0         = SIGMAT
        PORE0          = PORE
        DIVU0          = DIVU
        MASCN0         = MASCN
        !
        IF (TypeProcess.EQ.'MANDEL') THEN
            TMANDEL = TMANDEL + GEOTIME
            CALL GEOMECHANIC('MANDEL_DATA_EXAMPL')
        ENDIF
        !
        !...  ..COMPUTE INJECTION GEOMECHANICAL TIME
        TIMEINJ   = TIMELOAD(GEOTIME*DFLOAT(NNP))
        !...  ..BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND TRANSPORT (INDEX L)
        VELOCITYL = velocLadal
        INDEXL    = 0
        JUMPINDXL = .FALSE.
100     CONTINUE
        !INDEXL =  INDEXL+1
        DO  INDEXL=1,NITHIDRO
            !IF (INDEXL.EQ.NITHIDRO) JUMPINDXL = .TRUE.
            !...  .. UPDATE ARRAY SATURATIONS FOR ITERATIVE HIDRO-TRANSPORT
            satElemL0 = satElemL
            satElemL  = satElem
            satElem   = satElem0
            !
            IF (INDEXL.EQ.1) satElemL = satElemL0
            !...  .. BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND GEOMECHANIC (INDEX K)
            SIGMAK  = SIGMAT
            INDEXK  = 0
            do INDEXK = 1, NITGEO
                IF (INDEXK.EQ.NITGEO) JUMPINDXK = .TRUE.
                IF (INDEXK.EQ.1     ) SIGMAT= SIGMA0
                WRITE(*,2000) INDEXL, INDEXK, NITGEO
                PHI  = PORE
                PHI0 = PHIEULER
                call hidroGeomecanicaRT(satElemL,satElemL0,SIGMAT,SIGMA0,TIMEINJ)
                !...  .. ...
                CALL GEOMECHANIC('RIGHT_SIDE_2_SOLVE')
                CALL UPDTINIT(DIS,DIS0,NDOFD,NUMNP,INITS3)
                CALL GEOMECHANIC('ELAST_STRSS_SIGMAT')
                CALL GEOMECHANIC('RESETS_FORCE_VECTR')
                !.... ..
                ERRSZSIG = SIGNORM(SIGMAT,SIGMAK,numelReserv)
                IF (INDEXK.EQ.1) THEN
                    ERRSIG0 = ERRSZSIG
                ELSE
                    ERRSIG0  = DMAX1(ERRSIG0,ERRSZSIG)
                ENDIF
                ERRSZSIG = ERRSZSIG/ERRSIG0
                WRITE(*,3000) 'SIGMA',INDEXK, ERRSZSIG
                IF (ERRSZSIG.LT.TOLSIGMA) exit
            end do   !.... do INDEXK = 1, NITGEO
            CALL GEOMECHANIC('UPDAT_MASS_CONTENT')
            !          IF (NSD.EQ.2) THEN
            call transport(velocLadal,GEOTIME)
            ERRSZVEL = VELNORM(velocLadal,VELOCITYL,numLadosReserv)
            IF (INDEXL.EQ.1) THEN
                ERRVEL0 = ERRSZVEL
            ELSE
                ERRVEL0  = DMAX1(ERRVEL0,ERRSZVEL)
            ENDIF
            ERRSZVEL = ERRSZVEL/ERRVEL0
            WRITE(*,3000) 'VELOC',indexl,ERRSZVEL
            VELOCITYL = velocLadal
            !...
            IF ((ERRSZVEL.LT.TOLVELOC)) exit ! .AND.(.NOT.JUMPINDXL)) GOTO 100
        end do  ! END LOOP FOR INDEX L: ITERATION HIDRODINAMICS AND TRANSPORT
        !...
        tTransporte=tTransporte+dtBlocoTransp
        !
        call calc_prod(ndofV,numLadosElem,numelReserv,&
            &     nsd,nen, conecNodaisElem, conecLadaisElem,velocLadal, &
            &     satElem,x,tTransporte)
        !
        if (iflag_sat==1) then
            if (iflag_tipoPrint==1) then
                call gerarLabel(labelTransp,tTransporte)
                call escreverArqParaviewIntermed(isat, satElem, ndofV, &
                    numel, trim(labelTransp), len(trim(labelTransp)))
            endif
        endif
        !
        !.... imprime solucao intermediaria no tempo
        !
        call imprimirSolucaoNoTempo(satElem,DIS,PORE,YOUNG,GEOPRSR,pressaoElem, velocLadal, velocNodal, &
            &       VC, AVSTRS, MASCN, tTransporte, NDOFV,NDOFP, NDOFD, nrowb)
        !
        IF(iflag_tipoPrint.EQ.3)THEN
            IF (NUMDX.GT.0) THEN
                IF (MOD(NNP,NUMDX).EQ.0) THEN
                    CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
                        &                       VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
                ENDIF
            ENDIF
        END IF
        !
    END DO ! nnp=1,nvel
    !
    !.....CLOSE SERIES OF DATA AT NODAL POINTS
    !
    IF(iflag_tipoPrint.EQ.3)THEN
        IF(NUMDX.GT.0) CALL PRINT_DXINFO('CLOSE_FEDX_FILE',IFEDX, &
            &                                     NUMNP,NUMNP)
    END IF

    tempoTotalVelocidade=tempoMontagemVel+tempoSolverVel
    tempoTotalGeomecanica=tempoMontagemGeo+tempoSolverGeo
    !
    write(*,*) " "
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da velocidade     = ", tempoTotalVelocidade
    write(*,*) "      Montagem velocidade  =", tempoMontagemVel
    write(*,*) "      Solver velocidade    =", tempoSolverVel
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da pressao        = ", tempoTotalPressao
    write(*,*) "**********************************************"
    write(*,*) "Tempo total do transporte  = ", tempoTotalTransporte
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da geomecanica    = ", tempoTotalGeomecanica
    write(*,*) "      Montagem geomecanica =", tempoMontagemGeo
    write(*,*) "      Solver geomecanica   =", tempoSolverGeo
    write(*,*) "**********************************************"
    write(*,*) "SOMATORIO DOS TEMPOS =", tempoTotalVelocidade + &
        tempoTotalPressao    + &
        tempoTotalTransporte + &
        tempoTotalGeomecanica,  " segundos "

    return

1000 FORMAT(/,'###########################################',/, &
        '###########################################',/, &
        'tempo inicial: ',1PE15.8,5x,'Passo: ',i5, 2x,'de ',i5,/, &
        '###########################################',/)
1010 format('###########',/,'Fim da realizacao:', &
        & i5,/,'###########',/)
2000 FORMAT('ITERATIVE L-INDEX = ',I2,2X,'ITERATIVE K-INDEX =',I2,' OF ',I2)
2500 FORMAT('ITERATIVE WAY_S COUNTER = ',I2,' OF ',I2)
3000 FORMAT(/'PROPORTIONAL ERROR OF ',A5,' ON ITERATION ',I2,X,'IS ',1PE15.8/)
    !
    END SUBROUTINE
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
    !**** new *******************************************************************
    !
    subroutine processamento_creep()
    use mLeituraEscrita, only: iflag_tipoPrint

    use mGlobaisEscalares
    use mGeomecanica, only: ndofD
    use mHidrodinamicaRT, only: ndofV, ndofP
    !
    use mLeituraEscritaSimHidroGeoMec,   only : imprimirCondicoesIniciais,imprimirSolucaoNoTempo
    use mLeituraEscritaSimHidroGeoMec,   only : escreverArqParaviewIntermed
    use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXINFO, IFEDX
    use mLeituraEscritaSimHidroGeoMec,   only : isat, iflag_sat
    !
    use mMalha,            only : nsd, numel, numelReserv, nen
    use mMalha,            only : numnp
    use mMalha,            only : numLadosReserv, numLadosElem
    use mMalha,            only : x
    use mMalha,            only : conecNodaisElem, conecLadaisElem
    !
    use mPropGeoFisica,    only : phi, phi0
    use mPropGeoFisica,    only : permkx
    use mPropGeoFisica,    only : lerPropriedadesFisicas
    use mPropGeoFisica,    only : DTCREEP, TOLCREEP
    use mPropGeoFisica,    only : PORE, PORE0, PHIEULER
    use mPropGeoFisica,    only : YOUNG, MASCN0, MASCN
    !
    use mTransporte,       only : transport, satElem, calc_prod
    use mTransporte,       only : satElemL, satElemL0, satElem0
    use mGeomecanica,      only : SIGMAT, SIGMA0, GEOTIME, RESMAX
    use mGeomecanica,      only : NWTNITER, RESIDUAL
    use mGeomecanica,      only : NROWB
    use mGeomecanica,      only : GEOSETUP, DIS, DIS0, AVSTRS, PRINT_DX
    use mGeomecanica,      only : DIVU0, DIVU, DTRL, STRSS0, GEOPRSR
    use mGeomecanica,      only: geomechanic
    !
    use mHidrodinamicaRT,  only : hidroGeomecanicaRT, vc, velocLadal, velocNodal
    use mHidrodinamicaRT,  only : pressaoElem, pressaoElemAnt
    use mHidroDinamicaRT,  only : NPRESPRODWELL,PRESMEDIAINICIAL,PRESPROD
    !
    implicit none
    !
    !.... solution driver program
    !
    LOGICAL       :: JUMPINDXK, JUMPINDXL, LJUMP
    INTEGER       :: INDEXL, INDEXK
    character(21) :: labelTransp
    !
    REAL(8)       :: ERRSZSIG, ERRSZVEL, ERRSIG0, ERRVEL0
    REAL(8)       :: ERRSIZE
    REAL(8)       :: SIGNORM, VELNORM
    REAL*8        :: TIMEINJ, TIMELOAD, AUX
    REAL(8), external :: DESPRESSURIZAR
    REAL(8), external :: DESPRESSURIZAR_INIT
    REAL(8), DIMENSION(numelReserv)      :: SIGMAK
    REAL(8), DIMENSION(1,numLadosReserv) :: VELOCITYL
    !
    LJUMP = .FALSE.
    !

    !.... imprime condicoes iniciais
    !
    call imprimirCondicoesIniciais(GEOPRSR, pressaoElem, velocLadal, velocNodal, vc, phi, &
        permkx, satElem, YOUNG, DIS, STRSS0, MASCN, ndofV, ndofP, ndofD, nrowb)
    !
    IF(iflag_tipoPrint.EQ.3)THEN
        IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN
                CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
                    &                       VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
        ENDIF
    END IF
    !
    !.... tempo de simulacao para cada bloco do transporte
    !
    dtBlocoTransp=tt/dfloat(nvel)
    !
    !.... INICIALIZACAO DAS VARIAVEIS DA GEOMECANICA
    !
    GEOTIME = TT/DFLOAT(NVEL)
    DTCREEP = DFLOAT(NCREEP)*GEOTIME
    !
    !-----------------------------------------------------------------------
    !
    !.... DIMINUICAO GRADUAL DA PRESSAO NO POCO
    !
    AUX = DESPRESSURIZAR_INIT(PRESPROD,dtBlocoTransp,nvel)
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     SIMULACAO
    !
    !-----------------------------------------------------------------------
    !
    DO NNP=1,NVEL
        !
        ! DESPRESURIZACAO LENTA DO POCO DE PRODUCAO
        PRESPROD = DESPRESSURIZAR(NPRESPRODWELL,NNP, &
            &                              PRESMEDIAINICIAL,AUX,PRESPROD)
        write(*,1000) tTransporte,nnp,nvel
        !...  ..SETUP SATURATION FOR ITERATIVE HIDRODINAMICS AND TRANSPORT
        satElem0       = satElem
        satElemL       = satElem
        !...  ..UPDATE FOR MACROTIME EVOLUTION
        pressaoElemAnt = pressaoElem(1,:)
        SIGMA0         = SIGMAT
        PORE0          = PORE
        DIVU0          = DIVU
        MASCN0         = MASCN
        !...  ..COMPUTE INJECTION GEOMECHANICAL TIME
        TIMEINJ        = TIMELOAD(GEOTIME*DFLOAT(NNP))
        !...  ..BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND TRANSPORT (INDEX L)
        VELOCITYL = velocLadal
        INDEXL    = 0
        JUMPINDXL = .FALSE.
100     CONTINUE
        INDEXL =  INDEXL+1
        !            DO 500 INDEXL=1,NITHIDRO
        IF (INDEXL.EQ.NITHIDRO) JUMPINDXL = .TRUE.
        !...  .. UPDATE ARRAY SATURATIONS FOR ITERATIVE HIDRO-TRANSPORT
        satElemL0 = satElemL
        satElemL  = satElem
        satElem   = satElem0
        !
        IF (INDEXL.EQ.1) satElemL = satElemL0
        !...  .. BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND GEOMECHANIC (INDEX K)
        SIGMAK  = SIGMAT
        INDEXK  = 0
        JUMPINDXK = .FALSE.
300     CONTINUE
        do INDEXK = 1, NITGEO
            !INDEXK = INDEXK+1
            !               DO 400 INDEXK=1,NITGEO
            IF (INDEXK.EQ.NITGEO) JUMPINDXK = .TRUE.
            IF (INDEXK.EQ.1) SIGMAT = SIGMA0
            RESMAX   = 1.0D0
            NWTNITER = 0
            WRITE(*,2000) INDEXL, INDEXK, NITGEO
            PHI  = PORE
            PHI0 = PHIEULER
            call hidroGeomecanicaRT(satElemL,satElemL0,SIGMAT, &
                &              SIGMA0,TIMEINJ)
            !...  .. ...
350         CONTINUE
            !...  .. ...
            NWTNITER = NWTNITER + 1
            WRITE(*,2500) INDEXK,NITGEO,NWTNITER
            !...  .. ...
            IF (NWTNITER.GT.MAXITERC) LJUMP = .TRUE.
            !
            CALL GEOMECHANIC('GEOMECHANICS_CREEP')
            ERRSIZE = residual/resmax
            write(*,5000) ERRSIZE
            !
            !.... .. ... NEXT TEST RESIDUAL NORM WITH TOLERANCE CRITERIA 4 CREEP
            !
            IF ((ERRSIZE.GT.TOLCREEP).AND.(.NOT.LJUMP)) GOTO 350
            !.... .. ...
            !.... .. ...  UPDATE NON-LINEAR DISPLACEMENT
            !.... .. ...
            DIS = DTRL
            !
            !.... .. ...  CORRECTION WITH INITIAL DISPLACEMENTS
            !
            CALL UPDTINIT(DIS,DIS0,NDOFD,NUMNP,INITS3)
            !.... .. ...
            call GEOMECHANIC('CREEP_STRSS_SIGMAT')
            !.... .. ...
            ERRSZSIG = SIGNORM(SIGMAT,SIGMAK,numelReserv)
            IF (INDEXK.EQ.1) THEN
                ERRSIG0 = ERRSZSIG
            ELSE
                ERRSIG0  = DMAX1(ERRSIG0,ERRSZSIG)
            ENDIF
            ERRSZSIG = ERRSZSIG/ERRSIG0
            WRITE(*,3000) 'SIGMA',INDEXK, ERRSZSIG
            SIGMAK = SIGMAT
            !ERRSZSIGPreviews=ERRSZSIG
            IF (ERRSZSIG.LT.TOLSIGMA) exit
            !IF ((ERRSZSIG.GT.TOLSIGMA).OR.(.NOT.JUMPINDXK)) GOTO 300
        end do
        !IF ((ERRSZSIG.GT.TOLSIGMA).AND.(.NOT.JUMPINDXK)) GOTO 300
        ! 400       CONTINUE ! END LOOP FOR INDEX_K: ITERATION HIDRODINAMICS AND GEOMECHANIC
        !.... .. ...
        call GEOMECHANIC('UPDAT_MASS_CONTENT')

        !...
        call transport(velocLadal,GEOTIME)

        ERRSZVEL = VELNORM(velocLadal,VELOCITYL,numLadosReserv)
        IF (INDEXL.EQ.1) THEN
            ERRVEL0 = ERRSZVEL
        ELSE
            ERRVEL0  = DMAX1(ERRVEL0,ERRSZVEL)
        ENDIF
        ERRSZVEL = ERRSZVEL/ERRVEL0
        WRITE(*,3000) 'VELOC',indexl,ERRSZVEL
        VELOCITYL = velocLadal
        !...
        IF ((ERRSZVEL.GT.TOLVELOC).AND.(.NOT.JUMPINDXL)) GOTO 100
        !...
        !500   CONTINUE   ! END LOOP FOR INDEX L: ITERATION HIDRODINAMICS AND TRANSPORT
        !...
        tTransporte=tTransporte+dtBlocoTransp
        !
        call calc_prod(ndofV,numLadosElem,numelReserv, &
            &      nsd,nen, conecNodaisElem, conecLadaisElem,velocLadal, &
            &      satElem, x,tTransporte)
        !
        if (iflag_sat==1) then
            if (iflag_tipoPrint==1) then
                call gerarLabel(labelTransp,tTransporte)
                call escreverArqParaviewIntermed(isat, satElem, ndofV, &
                    numel, trim(labelTransp), len(trim(labelTransp)))
            endif
        endif

        !
        !.... imprime solucao intermediaria no tempo
        !
        call imprimirSolucaoNoTempo(satElem,DIS,PORE,YOUNG,GEOPRSR,pressaoElem, velocLadal, velocNodal, &
            &       VC, AVSTRS, MASCN, tTransporte, NDOFV,NDOFP, NDOFD, nrowb)

        !
        IF(iflag_tipoPrint.EQ.3)THEN
            IF (NUMDX.GT.0) THEN
                IF (MOD(NNP,NUMDX).EQ.0) THEN
                    CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
                        &                       VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
                ENDIF
            ENDIF
        END IF
        !
    end do ! nnp=1,nvel
    !
    !.....CLOSE SERIES OF DATA AT NODAL POINTS
    !
    IF(iflag_tipoPrint.EQ.3)THEN
        IF(NUMDX.GT.0) CALL PRINT_DXINFO('CLOSE_FEDX_FILE',IFEDX, &
            &                                    NUMNP,NUMNP)
    END IF


    tempoTotalVelocidade=tempoMontagemVel+tempoSolverVel
    tempoTotalGeomecanica=tempoMontagemGeo+tempoSolverGeo
    !
    write(*,*) " "
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da velocidade=", tempoTotalVelocidade
    write(*,*) "      Montagem velocidade=", tempoMontagemVel
    write(*,*) "      Solver velocidade=", tempoSolverVel
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da pressao   =", tempoTotalPressao
    write(*,*) "**********************************************"
    write(*,*) "Tempo total do transporte=", tempoTotalTransporte
    write(*,*) "**********************************************"
    write(*,*) "Tempo total da geomecanica=", tempoTotalGeomecanica
    write(*,*) "      Montagem geomecanica=", tempoMontagemGeo
    write(*,*) "      Solver geomecanica=", tempoSolverGeo
    write(*,*) "**********************************************"
    write(*,*) "TEMPO TOTAL DE EXECUCAO=", &
        tempoTotalVelocidade+tempoTotalPressao+tempoTotalTransporte+tempoTotalGeomecanica

    !
1000 FORMAT(/,'###########################################',/, &
        '###########################################',/, &
        'tempo inicial: ',1PE15.8,5x,'Passo: ',i5, 2x,'de ',i5,/, &
        '###########################################',/)
1010 format('###########',/,'Fim da realizacao:', &
        & i5,/,'###########',/)
2000 FORMAT('ITERATIVE L-INDEX = ',I2,2X,'ITERATIVE K-INDEX =',I2,' OF ',I2)
2500 FORMAT('ITERATIVE WAY_S COUNTER = ',I2,' OF ',I2,2X,&
        &        '    NEWTON ITERATION COUNTER =',I5)
3000 FORMAT(/'PROPORTIONAL ERROR OF ',A5,' ON ITERATION ',I2,X,'IS ',1PE15.8/)
5000 FORMAT('     RESIDUAL/RMAX = ',1PE15.8//)


    return
    end subroutine processamento_creep
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
    use mHidrodinamicaRT,  only: ndofV, nlvectV, ndofP

    !
    use mMalha,            only: nsd, numel, numelReserv, numnp, numnpReserv
    use mMalha,            only: numLadosReserv, nen, numLadosElem
    use mMalha,            only: x, xc
    !
    use mMalha,            only: listaDosElemsPorNo
    use mMalha,            only: listaDosElemsPorFace
    use mMalha,            only: conecNodaisElem, conecLadaisElem
    use mHidroDinamicaRT,  only: idVeloc, lmV
    use mTransporte,       only: satElemAnt, satElem
    use mTransporte,       only: satElemL, satElemL0, satElem0
    use mGeomecanica,      only: EINELAS
    use mGeomecanica,      only: nrowB, nintD, nrowB2
    use mGeomecanica,      only: idDesloc, lmD, ndofD, nlVectD, fDesloc
    use mGeomecanica,      only: idDis, dis, dis0, vdP, dtrl
    use mGeomecanica,      only: divU, divU0, strss, strss0, hmTTG, eCreep
    use mGeomecanica,      only: geoPrsr, avStrs, avCrep, strs3D
    use mGeomecanica,      only: sigmaT, sigma0
    use mPropGeoFisica,    only: GEOFORM, MASCN, MASCN0
    use mHidrodinamicaRT,  only: pressaoElem, pressaoElemAnt,PRESPRODWELL,NCONDP
    use mHidrodinamicaRT,  only: vc, ve, velocLadal, fVeloc, velocNodal
    
    use mGeomecanica, only: disInc

    !
    implicit none
    !
    !campos
    !
    allocate(pressaoElem(ndofP,numelReserv))
    pressaoElem = 0.0D0
    allocate(pressaoElemAnt(numelReserv))
    pressaoElemAnt = 0.0D0
    allocate(velocLadal (ndofV,numLadosReserv))
    velocLadal = 0.0D0
    ALLOCATE(PRESPRODWELL(NCONDP))
    PRESPRODWELL = 0.0D0
    allocate(ve (nsd,nen,numelReserv))
    ve = 0.0D0
    allocate(vc (nsd,numelReserv))
    vc = 0.0D0
    allocate(velocNodal(nsd,numnpReserv))
    velocNodal = 0.0D0
    allocate(satElem(numelReserv))
    satElem = 0.0D0
    allocate(satElemAnt(numelReserv))
    satElemAnt  = 0.0D0
    !
    ! NEXT 3 LINES NEW SATURATION ARRAY FOR ITERATIVE HIDRODINAMICS AND TRANSPORT
    !
    ALLOCATE(satElem0(numelReserv))
    satElem0    = 0.0D0
    ALLOCATE(satElemL(numelReserv))
    satElemL    = 0.0D0
    ALLOCATE(satElemL0(numelReserv))
    satElemL0   = 0.0D0
    !
    !malha
    !
    allocate(x(nsd,numnp))
    x  = 0.0d0
    allocate(xc(nsd,numel))
    xc = 0.0d0
    allocate(conecNodaisElem(nen,numel))
    conecNodaisElem = 0
    allocate(conecLadaisElem(numLadosElem,numelReserv))
    conecLadaisElem = 0
    allocate(listaDosElemsPorNo(nen,numnp))
    listaDosElemsPorNo = 0
    allocate(listaDosElemsPorFace(numLadosElem,numLadosReserv))
    listaDosElemsPorFace = 0
    !
    !contorno
    !
    allocate(idVeloc(ndofV,numLadosReserv))
    idVeloc = 0
    allocate(lmV(ndofV,numLadosElem,numelReserv))
    lmV = 0
    !
    if (nlvectV.ne.0)  then
        allocate(fVeloc(ndofV,numLadosReserv))
        fVeloc = 0.0d0
    endif
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
    ALLOCATE(HMTTG (NUMEL,NINTD,NROWB2)); HMTTG   = 0.0D0
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
    allocate(disInc(ndofD, numnp)); disInc = 0.d0
    !
    !... END GEOMECANICA
    !
    allocate(beta(numel))
    !
    end subroutine

    !
    !**** NEW **** MODIFIED 4 HIERARCH MESH  ************************************
    !
    SUBroutine TOPologiaMALhaSistEQUAcoesDS(NALHSV,NEQV,NALHSD,NEQD)
    use mGlobaisEscalares
    use mGeomecanica, only: ndofD, optSolverD
    use mHidrodinamicaRT, only: ndofV, optSolverV
    use mGlobaisArranjos
    use mLeituraEscrita, only: iecho, nprint, prntel
    use mLeituraEscritaSimHidroGeoMec, only: GEOREGION_DS
    !
    use mHidroDinamicaRT, only: lmV, idVeloc, idiagV, nedV
    use mHidroDinamicaRT, only: meanbwV
    use mGeomecanica,     only: LMD, IDIAGD, IDDESLOC
    use mGeomecanica,     only: AUXM,IOPT,NED,NED2,NEE,NEE2,NEESQ,NEESQ2,NESD,NSTR
    !
    use mMalha,          only: numel, numnp, nsd, nen
    use mMalha,          only: numLadosElem,numLadosReserv
    use mMalha,          only: numelReserv,numnpReserv
    use mMalha,          only: x, xc, local
    use mMalha,          only: conecNodaisElem, conecLadaisElem
    use mMalha,          only: listaDosElemsPorNo
    use mMalha,          only: listaDosElemsPorFace
    use mMalha,          only: criarListaVizinhos
    use mMalha,          only: genel
    use mMalha,          only: formlm
    !
    use mPropGeoFisica,  only: hx,hy,hz, calcdim
    use mPropGeoFisica,  only: nelxReserv, nelyReserv

    use mInputReader,      only: readNodeElementsDS, readMaterialPropertiesDS, readConstantBodyForcesDS
    use mInputReader,      only: genelFacesDS, leituraGeracaoConectividadesDS
    !

    !
    implicit none
    !
    !.... program to set arrays storage
    !
    integer :: NALHSV, NEQV, NALHSD, NEQD
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
    nedV   = ndofV
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
    allocate(c(6,numat)); c=0.d0
    !
    write(iecho,1000) ntype,numel,numat,nen,npint
    !
    !      read material properties
    !
    keyword_name = "material_properties"
    call readMaterialPropertiesDS(keyword_name, iecho, ierr)
    !
    !     constant body forces
    !
    keyword_name = "constant_body_forces"
    call readConstantBodyForcesDS(keyword_name, iecho, ierr)
    !
    !    generation of conectivities
    !
    keyword_name = "conectividades_nodais"
    call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen, ierr)
    !
    keyword_name = "conectividades_ladais"
    if (nsd==2) then
        call leituraGeracaoConectividadesDS(keyword_name, conecLadaisElem, mat, numLadosELem, ierr)
    else
        call genelFacesDS(keyword_name, conecLadaisElem, numLadosElem, nelxReserv, nelyReserv, ierr)
    endif

    listaDosElemsPorNo=0
    call criarListaVizinhos(nen,numnp,numel,conecNodaisElem, &
        &     listaDosElemsPorNo  )
    listaDosElemsPorFace=0
    call criarListaVizinhos(numladosElem,numLadosReserv, &
        &     numelReserv,conecLadaisElem,listaDosElemsPorFace)
    !
    !.... generation of velocity lm array
    !
    if (NEQV > 0) then
        call formlm(idVeloc,conecLadaisElem,lmV,nedV,nedV,&
            &     numLadosElem,numelReserv)
        !
        !.... modification of idiag array
        !
        call colht(idiagV,lmV,nedV,numLadosElem,numelReserv,neqV)
        !
        if(optSolverV=='skyline') then
            call diag(idiagV,neqV,nalhsV)
        endif
        !
        meanbwV = nalhsV/neqV
        write(iecho,9000) neqV, nalhsV, meanbwV, 8.0*nalhsV/1000/1000
    end if

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

    if(optSolverV=='pardiso') then
        write(*,*) ' não implementado '
        stop 9
    endif
    !
    if(optSolverD=='pardiso') then
        write(*,*) ' não implementado '
        stop 9
    endif
    !
    print*, "NALHSV=", NALHSV, "NALHSD=", NALHSD

    if(optSolverV=="hypre") then
        write(*,*) ' não implementado '
        stop 9
    endif

    if(optSolverD=="hypre") then
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
    !
    !===========================================================================
    !
    FUNCTION DESPRESSURIZAR_INIT(BWP,DT,nvel)
    !
    use mHidroDinamicaRT,  only : pressaoElem
    use mHidroDinamicaRT,  ONLY : TPRESPRODWELL,NPRESPRODWELL,PRESMEDIAINICIAL
    use mHidroDinamicaRT,  ONLY : NCONDP,ELEM_CONDP
    !
    IMPLICIT NONE
    !
    REAL(8) :: BWP,DT,AUX
    REAL(8) :: DESPRESSURIZAR_INIT
    INTEGER :: NEL,NVEL
    !
    PRESMEDIAINICIAL = -1E30
    DO NEL=1,NCONDP
        IF(pressaoElem(1,ELEM_CONDP(NEL)).GE.PRESMEDIAINICIAL)THEN
            PRESMEDIAINICIAL=pressaoElem(1,ELEM_CONDP(NEL))
        end IF
    END DO
    WRITE(*,*)'#############################################################'
    WRITE(*,*)'PRESSAO INICIAL NO FUNDO DO POCO DE PRODUCAO: ###############'
    WRITE(*,*)PRESMEDIAINICIAL
    WRITE(*,*)'#############################################################'
    !
    NPRESPRODWELL = NINT(TPRESPRODWELL/DT)
    IF(TPRESPRODWELL.GT.REAL(NVEL)*DT)THEN
        AUX = 0.0
        NPRESPRODWELL=1
    ELSE
        IF(NPRESPRODWELL.EQ.0)THEN
            AUX = 0.0
        ELSE
            AUX = (PRESMEDIAINICIAL-BWP)/REAL(NPRESPRODWELL)
        END IF
    END IF
    !
    WRITE(*,*)'#############################################################'
    WRITE(*,*)'TEMPO TOTAL PARA CHEGAR A PRESSAO MINIMA NO POCO DE PRODUCAO:'
    WRITE(*,*)TPRESPRODWELL
    WRITE(*,*)'#############################################################'
    !
    DESPRESSURIZAR_INIT = AUX
    !
    END FUNCTION DESPRESSURIZAR_INIT

    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine initGalerkinPressure(nnp)
    !variable imports
    use mHidrodinamicaGalerkin, only: hgNumPassosTempo, hgTempoTotal
    use mHidrodinamicaGalerkin, only: hgInitialPressure, hgPrevPressure, hgPressure
    use mHidrodinamicaGalerkin, only: hgNdof, hgNlvect, hgNeq
    use mHidrodinamicaGalerkin, only: hgF
    use mHidrodinamicaGalerkin, only: hgId

    use mLeituraEscrita, only:iecho
    use mGlobaisEscalares, only: iprtin, betaCompressibility

    !function imports
    use mInputReader,      only: leituraCodigosCondContornoDS, leituraValoresCondContornoDS
    use mInputReader,      only: readRealKeywordValue, readIntegerKeywordValue

    implicit none

    !variable input
    integer nnp

    !variables
    character(len=50) :: keyword_name
    integer :: ierr

    !------------------------------------------------------------------------------------------------------------------------------------
    keyword_name = "numeroPassosTempoP"
    call readIntegerKeywordValue(keyword_name, hgNumPassosTempo, 0_4, ierr)

    keyword_name = "tempoTotal"
    call readRealKeywordValue(keyword_name, hgTempoTotal, 0.0d0, ierr)
    hgTempoTotal = hgTempoTotal * 2592000.0 ! converts months to seconds

    keyword_name = "hgCompressibilidade"
    call readRealKeywordValue(keyword_name, betaCompressibility, 0.0d0, ierr)

    keyword_name = "hgNlvect"
    call readIntegerKeywordValue(keyword_name, hgNlvect, 0_4, ierr)

    keyword_name = "pressao_inicial_fluxo"
    call readRealKeywordValue(keyword_name, hgInitialPressure, 0.0d0, ierr)


    hgNdof = 1
    allocate(hgPrevPressure(hgNdof, nnp))
    hgPrevPressure = 0.0d0
    allocate(hgPressure(hgNdof, nnp))
    hgPressure=0.0d0

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
    subroutine processamentoGalerkinElastico()
    !variable imports
    !flow
    use mHidrodinamicaGalerkin, only:hgNumPassosTempo, hgTempoTotal, hgInitialPressure
    use mHidrodinamicaGalerkin, only:hgPrevPressure, hgPressure
    use mHidrodinamicaGalerkin, only:hgNdof
    use mMalha, only: numnp, numel, nen, nsd, conecNodaisElem
    use mMalha, only: x
    use mLeituraEscritaSimHidroGeoMec, only:ihgPres, reservHGPres
    !mechanics
    use mGeomecanica, only: dis0
    use mGeomecanica, only: disInc
    use mGeomecanica, only: ndofD
    use mLeituraEscritaSimHidroGeoMec, only:iDis, reservDesloc

    !function imports
    use mHidrodinamicaGalerkin, only:montarEstruturasDadosPressaoSkyline
    use mHidrodinamicaGalerkin, only:montarSistemaEquacoesPressao
    use mHidrodinamicaGalerkin, only:solveGalerkinPressao
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaview, escreverArqParaviewIntermed
    use mLeituraEscritaSimHidroGeoMec, only:escreverArqParaviewVector, escreverArqParaviewIntermed_CampoVetorial
    use mGeomecanica, only:stress_init
    use mGeomecanica, only:incrementMechanicSolution

    !variables
    implicit none
    integer :: curTimeStep, i, j
    real*8 :: t, deltaT
    character(21) :: labelTempo, num

    !------------------------------------------------------------------------------------------------------------------------------------
    if(hgNumPassosTempo > 1) then
        ! inicializar a solução no passo anterior
        hgPrevPressure = hgInitialPressure
    endif

    ! print initial state
    call escreverArqParaview(ihgPres, hgPrevPressure, hgNdof, numnp, nen, conecNodaisElem, 2, 't=0.0', len('t=0.0'), reservHGPres)
    
    ! init stress and print stress
    call stress_init()
    call escreverArqParaviewVector('dis', dis0, nDofD, numnp, nen, conecNodaisElem, 2, 'total', len('t=0.0'), reservDesloc, iDis)
    
    ! solve initialization incrementally
    call incrementMechanicSolution(conecNodaisElem, nen, numel, numnp, nsd, x, disInc)
    call escreverArqParaviewIntermed_CampoVetorial('dis',disInc, ndofD, numnp, 'incrFinal', len('incrFinal'), 2, reservDesloc, iDis)

    ! mount data structure
    call montarEstruturasDadosPressaoSkyline(conecNodaisElem, nen, numel)

    !time loop
    t = 0
    deltaT = hgTempoTotal / hgNumPassosTempo
    do curTimeStep = 1, hgNumPassosTempo
        t = t + deltaT
        write(*,*) "tempo", t, "passo de tempo", curTimeStep

        ! mount equation system
        call montarSistemaEquacoesPressao(curTimeStep, conecNodaisElem, nen, numel, numnp, nsd, x)

        !solve
        call solveGalerkinPressao(curTimeStep, numnp)

        !print current solution
        write(num,'(f12.5)') t / 2592000
        labelTempo="t="//ADJUSTL(num)
        call escreverArqParaviewIntermed(ihgPres, hgPressure, hgNdof, numnp, trim(labelTempo), len(trim(labelTempo)))

        ! copy solution to previous array
        do j = 1, numnp
            do i = 1, hgNdof
                hgPrevPressure(i,j) = hgPressure(i,j)
            end do
        end do
    end do

    end subroutine processamentoGalerkinElastico
    !************************************************************************************************************************************
    !************************************************************************************************************************************