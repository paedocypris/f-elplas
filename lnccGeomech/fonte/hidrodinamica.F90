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
    module mHidrodinamicaGalerkin

    implicit none
    public montarEstruturasDadosPressaoSkyline, montarSistemaEquacoesPressao, &
        & solveGalerkinPressao, incrementFlowPressureSolution

    !module variables
    real*8 :: hgTempoTotal
    integer, allocatable :: hgNSteps(:)
    real*8, allocatable :: hgDts(:)
    integer :: hgNTimeSteps
    real*8 :: hgTimeConvertion

    !various integers
    integer :: hgNeq, hgNdof
    integer :: hgNlvect, hgNltftn, hgNptslf

    !global arrays
    real*8, allocatable :: hgF(:,:,:), hgG(:,:,:)

    !flow variables
    real*8 :: hgInitialPressure

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
    if(.not.allocated(hgLm)) allocate(hgLm(hgNdof, nen, nel))
    hgLm = 0.0
    call formlm(hgId, conecNodaisElem, hgLm, hgNdof, hgNdof ,nen, nel)

    !compute column heights in global left-hand-side matrix
    if(.not.allocated(hgIdiag)) allocate(hgIdiag(hgNeq))
    hgIdiag = 0.0
    call colht(hgIdiag, hgLm, hgNdof, nen, nel, hgNeq)

    !compute diagonal addresses of left-hand-side matrix
    call diag(hgIdiag, hgNeq, hgNalhs)

    !initialize alhs
    if(.not.allocated(hgAlhs)) allocate(hgAlhs(hgNalhs))

    !initialize brhs
    if(.not.allocated(hgBrhs)) allocate(hgBrhs(hgNeq))

    end subroutine montarEstruturasDadosPressaoSkyline
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine montarSistemaEquacoesPressao(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, tangentMatrix, out2waySource, way, plastType)

    !variable input
    integer ::conecNodaisElem(nen,nel), nen, nel, nnp, nsd
    real*8 :: x(nsd,nnp), curT, deltaT
    real*8 :: p(hgNdof,nnp), prevP(hgNdof,nnp)
    real*8 :: stressS(:, :), prevStressS(:, :), trStrainP(:,:), prevTrStrainP(:,:)
    real*8 :: tangentMatrix(:,:,:)
    real*8 :: out2waySource(:,:)
    integer :: way, plastType
    
    !variables
    real*8 :: g1(hgNltftn)

    !------------------------------------------------------------------------------------------------------------------------------------
    hgAlhs = 0.0
    hgBrhs = 0.0
    
    if (hgNlvect >= 1) then
        if (hgNltftn >= 1) then
            call lfac(hgG, curT, g1, hgNltftn, hgNptslf)
        else
            g1 = 1.d0
        endif
        
        call loadtime(hgId, hgF, hgBrhs, g1, hgNdof, nnp, hgNlvect,deltaT)
        call ftod(hgId, p, hgF, g1, hgNdof, nnp, hgNlvect)
    endif

    call calcCoeficientesSistemaEquacoesPressao(conecNodaisElem, nen, nel, nnp, nsd, x, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, tangentMatrix, out2waySource, way, plastType)

    end subroutine montarSistemaEquacoesPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine calcCoeficientesSistemaEquacoesPressao(conecNodaisElem, nen, nel, nnp, nsd, x, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, tangentMatrix, out2waySource, way, plastType)
    !variable imports
    use mGlobaisEscalares, only: nrowsh, npint

    !function imports
    use mFuncoesDeForma, only: shlt, shlq, shlq3d
    use mFuncoesDeForma, only: shgq, shg3d
    use mSolverGaussSkyline, only: addlhs, addrhs
    use mMalha, only: local
    use mPropGeoFisica, only: calcMatrixBulk, calcBulkFromMatrix, calcBiotCoefficient
    use mPropGeoFisica, only: totalCompressibility
    use mPropGeoFisica, only: calcKozenyCarmanPerm, calcKozenyCarmanPermInitialPoro

    !variable input
    integer :: conecNodaisElem(nen,nel), nen, nel, nnp, nsd
    real*8 :: x(nsd, nnp)
    real*8 :: deltaT
    real*8 :: p(hgNdof, nnp), prevP(hgNdof, nnp)
    real*8 :: stressS(npint, nel), prevStressS(npint, nel)
    real*8 :: out2waySource(npint, nel)
    real*8 :: trStrainP(npint, nel), prevTrStrainP(npint, nel)
    real*8 :: tangentMatrix(16,npint,nel)
    integer :: way, plastType

    !variables
    integer :: nee ! number of element equations
    real*8 :: elementK(nen, hgNdof * nen), elementF(hgNdof * nen) !element K matrix, element f matrix

    real*8 :: shG(nrowsh,nen,npint), shL(nrowsh,nen,npint) !global shape function with derivative, local shape function with derivative
    real*8 :: det(npint), w(npint) ! jacobian determinant, gauss integration weight

    real*8 :: xL(nsd, nen), pressureL(1, nen), prevPressureL(1, nen) !local position, local pressure, previous local pressure
    real*8 :: kX, kY !x component of permeability, y component of permeability

    integer :: ni, nj !node position considering multiple degrees of freedom
    real*8 :: prevPLInt
    real*8 :: djx, djy, djn, dix, diy, din

    integer :: curElement, l, i, j ! iterators
    real*8 :: cl !temporary variables

    logical :: quad, zerodl !is diagonal, is degenerated triange, is zero
    
    real*8 :: betaTotalCompressibility
    real*8 :: matrixBulk
    real*8 :: psi
    real*8, allocatable :: curPerm(:)
    integer :: tensDim

    !------------------------------------------------------------------------------------------------------------------------------------
    if (nsd > 2) then
        stop 9
    end if

    ! set some constants
    nee = hgNdof * nen
    
    if (nsd == 2) then
        tensDim = 3
    else
        tensDim = 6
    end if
    if (.not.allocated(curPerm)) allocate(curPerm(tensDim))

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
        
        call local(conecNodaisElem(1,curElement), x, xL, nen, nsd, nsd)
        call local(conecNodaisElem(1,curElement), p, pressureL, nen, hgNdof, hgNdof)
        call local(conecNodaisElem(1,curElement), prevP, prevPressureL, nen, hgNdof, hgNdof)

        ! check if element has any coalesced nodes
        quad = .true.
        if (nen.eq.4.and.conecNodaisElem(3,curElement).eq.conecNodaisElem(4,curElement)) quad = .false.

        ! calculates global derivatives of shape functions and jacobian determinants
        if(nen==3) call shgq  (xl,det,shL,shG,npint,curElement,quad,nen)
        if(nen==4) call shgq  (xl,det,shL,shG,npint,curElement,quad,nen)
        if(nen==8) call shg3d (xl,det,shL,shG,npint,curElement,nen)
        
        ! retrieve permeabilities from material
        curPerm = calcKozenyCarmanPermInitialPoro(curElement, tensDim)
        !curPerm = calcKozenyCarmanPerm(curElement,tensDim)
        kX = curPerm(1)
        kY = curPerm(2)

        do l = 1, npint ! foreach integration point
            cl = w(l)*det(l)
            
            !get mechanical parameters
            if (way == 1 .or. plastType == 1) then
                matrixBulk = calcMatrixBulk(curElement)
            else if (plastType == 2) then
                matrixBulk = calcBulkFromMatrix(tangentMatrix(1:16,l,curElement))
            end if
            betaTotalCompressibility = totalCompressibility(curElement, matrixBulk)
            psi = calcBiotCoefficient(curElement)/matrixBulk

            do i = 1, nen
                do j = 1, nen !for each pair of nodes in the element
                    nj = hgNdof * j
                    djx = shg(1, j, l)
                    djy = shg(2, j, l)
                    djn = shg(nrowsh, j, l)

                    ni = hgNdof * i
                    dix = shg(1, i, l)
                    diy = shg(2, i, l)
                    din = shg(nrowsh, i, l)

                    elementK(nj, ni) = elementK(nj, ni) + deltaT * Kx * dix * djx * cl !kn A
                    elementK(nj, ni) = elementK(nj, ni) + deltaT * Ky * diy * djy * cl !kn A
                    elementK(nj, ni) = elementK(nj, ni) + betaTotalCompressibility * din * djn * cl !B
                end do
            end do
            
            prevPLInt = dot_product(shg(nrowsh, 1:nen, l),prevPressureL(1,1:nen))
            
            do j = 1, nen !foreach node on element
                nj = hgNdof * j
                djn = shg(nrowsh,j,l)

                ! source terms
                elementF(nj) = elementF(nj) + (betaTotalCompressibility * djn * prevPLInt) * cl !B \csi^(n-1)
                if (way == 2) then
                    out2waySource(l, curElement) = (-psi) * (stressS(l,curElement) - prevStressS(l,curElement)) * djn * cl
                    
                    elementF(nj) = elementF(nj) + out2waySource(l, curElement)
                    if (plastType == 1) then
                        elementF(nj) = elementF(nj) - (trStrainP(l, curElement) - prevTrStrainP(l, curElement)) * djn * cl
                    end if
                end if
            end do
        end do

        ! computate dirichlet bc contribution to the system
        call zTest(pressureL, nee, zerodl)
        if(.not.zerodl) then
            call kdbc2(elementK,elementF, pressureL, nee, hgLm(1,1,curElement), nel)
        endif

        ! assemble element stifness matrix and force array into the global matrixes
        call addlhs (hgAlhs, elementK, hgLm(1,1,curElement), hgIdiag, nee, .false., .true.)
        call addrhs(hgBrhs, elementF, hgLm(1,1,curElement), nee)
    end do

    end subroutine calcCoeficientesSistemaEquacoesPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine solveGalerkinPressao(nnp, p)
    !function imports
    use mSolverGaussSkyline, only: factor, back
    !variables input
    integer :: nnp
    real*8 :: p(hgNdof, nnp)

    !------------------------------------------------------------------------------------------------------------------------------------
    !solve by LU decomposition
    call factor(hgAlhs, hgIdiag, hgNeq)
    call back(hgAlhs, hgBrhs, hgIdiag, hgNeq)

    !store the result
    call btod(hgId, p, hgBrhs, hgNdof, nnp)

    end subroutine solveGalerkinPressao
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine incrementFlowPressureSolution(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, vDarcy, vDarcyNodal, tangentMatrix, out2waySource, way, plastType)
    implicit none
    ! variables input
    integer :: conecNodaisElem(nen, nel), nen, nel, nnp, nsd
    real*8 :: x(nsd,nnp), curT, deltaT
    real*8 :: p(hgNdof,nnp), prevP(hgNdof,nnp)
    real*8 :: stressS(:,:), prevStressS(:,:)
    real*8 :: trStrainP(:,:), prevTrStrainP(:,:)
    real*8 :: vDarcy(:,:), vDarcyNodal(:,:)
    real*8 :: out2waySource(:,:)
    real*8 :: tangentMatrix(:,:,:)
    integer :: way, plastType
    
    !------------------------------------------------------------------------------------------------------------------------------------
    ! mount equation system
    call montarSistemaEquacoesPressao(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, tangentMatrix, out2waySource, way, plastType)
    
    !solve
    call solveGalerkinPressao(nnp, p)
    
    !process darcy velocity
    call calcDarcyVelocity(conecNodaisElem, nen, nel, nnp, nsd, p, vDarcy, vDarcyNodal)
    
    end subroutine incrementFlowPressureSolution
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine incrementFlowPressureSolutionOneWay(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, vDarcy, vDarcyNodal)
    ! variables input
    integer :: conecNodaisElem(nen, nel), nen, nel, nnp, nsd
    real*8 :: x(nsd,nnp), curT, deltaT
    real*8 :: p(hgNdof,nnp), prevP(hgNdof,nnp)
    real*8 :: vDarcy(:,:), vDarcyNodal(:,:)
    
    !variables
    real*8 :: stressS(1,1), prevStressS(1,1)
    real*8 :: trStrainP(1,1), prevTrStrainP(1,1)
    real*8 :: out2waySource(1,1)
    real*8 :: tangentMatrix(1,1,1)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    call incrementFlowPressureSolution(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, vDarcy, vDarcyNodal, tangentMatrix, out2waySource, 1, 1)
    
    end subroutine incrementFlowPressureSolutionOneWay
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine incrementFlowPressureSolutionTwoWay(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, vDarcy, vDarcyNodal, tangentMatrix, out2waySource, plastType)
    ! variables input
    integer :: conecNodaisElem(nen, nel), nen, nel, nnp, nsd
    real*8 :: x(nsd,nnp), curT, deltaT
    real*8 :: p(hgNdof,nnp), prevP(hgNdof,nnp)
    real*8 :: stressS(:,:), prevStressS(:,:)
    real*8 :: trStrainP(:,:), prevTrStrainP(:,:)
    real*8 :: vDarcy(:,:), vDarcyNodal(:,:)
    real*8 :: out2waySource(:,:)
    real*8 :: tangentMatrix(:,:,:)
    integer :: plastType
    
    !------------------------------------------------------------------------------------------------------------------------------------
    call incrementFlowPressureSolution(conecNodaisElem, nen, nel, nnp, nsd, x, curT, deltaT, p, prevP, stressS, prevStressS, trStrainP, prevTrStrainP, vDarcy, vDarcyNodal, tangentMatrix, out2waySource, 2, plastType)
    
    end subroutine incrementFlowPressureSolutionTwoWay
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine calcDarcyVelocity(conecNodaisElem, nen, nel, nnp, nsd, p, vDarcy, vDarcyNodal)
    !function imports
    use mFuncoesDeForma, only: shlt, shlq, shlq3d
    use mFuncoesDeForma, only: shgq, shg3d
    use mFuncoesDeForma, only: shlqen
    use mPropGeoFisica, only: calcKozenyCarmanPerm, calcKozenyCarmanPermInitialPoro
    use mMalha, only: local
    
    !variables import
    use mGlobaisEscalares, only: nrowsh, npint
    use mMalha, only: x
    
    implicit none
    !variables input
    integer :: conecNodaisElem(nen, nel), nen, nel, nnp, nsd
    real*8 :: p(hgNdof,nnp), vDarcy(nsd,nel), vDarcyNodal(nsd,nnp)
    
    !variables
    real*8 :: det(nen), w(npint)
    real*8 :: shL(nrowsh,nen,npint), shG(nrowsh,nen,npint)
    real*8 :: shLen(nrowsh,nen,nen), shGen(nrowsh,nen,nen)
    
    real*8 :: xL(nsd,nen), pressureL(hgNdof,nen)
    real*8 :: velL(nsd)
    
    real*8 :: vDarcyNodalAccum(nsd,nnp)
    integer :: globalNodeIndex, numNodalContributions(nsd,nnp)
    
    real*8 :: curPerm(3)
    
    logical :: quad
    real*8 :: gradP(nsd)
    integer :: curElement, l, dir, i
    real*8 :: area, c1
    
    !------------------------------------------------------------------------------------------------------------------------------------
    vDarcy = 0.
    vDarcyNodalAccum = 0.
    numNodalContributions = 0
    
    ! construct shape functions and calculates quadrature weights
    if(nen==3) call shlt(shL,w,npint,nen)
    if(nen==4) then
        call shlq(shL,w,npint,nen)
        call shlqen(shLen, nen)
    end if
    if(nen==8) call shlq3d(shL,w,npint,nen)
    
    !loop over elements
    do curElement = 1, nel
        !localize variables
        call local(conecNodaisElem(1,curElement), x, xL, nen, nsd, nsd)
        call local(conecNodaisElem(1,curElement), p, pressureL, nen, hgNdof, hgNdof)
        
        ! check if element has any coalesced nodes
        quad = .true.
        if (nen.eq.4.and.conecNodaisElem(3,curElement).eq.conecNodaisElem(4,curElement)) quad = .false.

        ! calculates global derivatives of shape functions and jacobian determinants
        if(nen==3) call shgq  (xl,det,shL,shG,npint,nel,quad,nen)
        if(nen==4) then 
            call shgq(xl,det,shL,shG,npint,nel,quad,nen)
            call shgq(xL,det,shLen,shGen,nen,nel,quad,nen)
        end if
        if(nen==8) call shg3d (xl,det,shL,shG,npint,nel,nen)

        ! set material index
        curPerm = calcKozenyCarmanPermInitialPoro(curElement,3) !Porosities are global variable, care
        !curPerm = calcKozenyCarmanPerm(curElement,3) !Porosities are global variable, care
        
        do dir=1,nsd
            area = 0.d0
            do l=1,npint
                c1 = det(l)*w(l)
                area = area + c1
                
                gradP(dir) = dot_product(shG(dir,1:nen,l), pressureL(1,1:nen))
                velL(dir) = -curPerm(dir) * gradP(dir)

                vDarcy(dir,curElement) = vDarcy(dir,curElement) + velL(dir)*c1
            end do
            vDarcy(dir,curElement) = vDarcy(dir,curElement)/area
            
            do i = 1, nen
                globalNodeIndex = conecNodaisElem(i,curElement)
                
                gradP(dir) = dot_product(shGen(dir,1:nen,i), pressureL(1,1:nen))
                velL(dir) = -curPerm(dir) * gradP(dir)
                
                vDarcyNodalAccum(dir,globalNodeIndex) = vDarcyNodalAccum(dir,globalNodeIndex) + velL(dir)
                numNodalContributions(dir,globalNodeIndex) = numNodalContributions(dir,globalNodeIndex) + 1
            end do
        end do
    end do
    
    do i = 1, nnp
        do dir = 1,nsd
            vDarcyNodal(dir,i) = vDarcyNodalAccum(dir,i)/numNodalContributions(dir,i)
        end do
    end do

    end subroutine calcDarcyVelocity
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    end module mHidrodinamicaGalerkin