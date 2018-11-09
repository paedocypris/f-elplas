    !     *********************
    !     *                                                                 *
    !     *            * * *    G E O M E C H A N I C      * * *            *
    !     *                                                                 *
    !     *                                                                 *
    !     *           A  FORTRAN 77  PROGRAM FOR SIMULATION                 *
    !     *                                                                 *
    !     *       GEOMECHANICAL BEHAVIOUR WITH HETEROGENEOUS YOUNG          *
    !     *                                                                 *
    !     *       MODULUS OF ELASTIC SATURATED POROUS MEDIA                 *
    !     *                                                                 *
    !     *       FRAMEWORK: FINITE ELEMENT AND B--BAR METHODS              *
    !     *                                                                 *
    !     *       BASED ON LINEAR CODE OF MARCIO MURAD - LNCC               *
    !     *                                                                 *
    !     *    --------------------------------------------------------     *
    !     *                                                                 *
    !     *       NON-LINEAR AND B-BAR FEM IMPLEMENTATION                   *
    !     *       FORMULATION JESUS ALEXEI - UFF-LNCC                       *
    !     *                                                                 *
    !     *                      MODIFY ON 30/11/2015                       *
    !     *                                                                 *
    !**** new module **********
    !
    module mGeomecanica
    !
    use mGlobaisEscalares, only: novaMalha

    implicit none

    integer              :: ndofD, nlvectD, nltftnD, nptslfD
    integer              :: neqD, nalhsD
    real*8,  allocatable :: alhsD(:), clhsD(:), brhsD(:)
    integer, allocatable :: idDesloc(:,:)
    integer, allocatable :: idiagD(:)
    integer, allocatable :: lmD(:,:,:)
    integer              :: numCoefPorLinhaGeo
    character(len=7)     :: optSolverD
    integer, allocatable :: AiGeo(:), ApGeo(:), LMstencilGeo(:,:)
    INTEGER              :: ptD(64), iparmD(64)
    REAL*8               :: dparmD(64)
    logical :: simetriaGeo
    !
    ! ESCALARES
    REAL(8) :: KGYNG, RHOYNG ! YOUNG's GEOMETRIC MEAN, STRENGHT
    INTEGER :: NDAUX, NDT
    REAL*8  :: TMANDEL, GEOTIME
    INTEGER :: NROWB, NINTD, NESD, NED, NED2, NSTR
    INTEGER :: NEESQ2, NEESQ, NEE, NEE2, IOPT, NLOOPS
    !
    INTEGER, ALLOCATABLE :: IDDIS(:,:)
    REAL(8), ALLOCATABLE :: GEOPRSR(:)
    REAL(8), ALLOCATABLE :: SIGMAT(:),SIGMA0(:),DIVU(:),DIVU0(:)
    REAL(8), ALLOCATABLE :: DIS(:,:), DIS0(:,:), VDP(:,:)
    real(8), allocatable :: STRSS(:,:,:), STRSSP(:,:,:), EINELAS(:,:,:)
    ! CREEP
    REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: DTRL, AVCREP, AVSTRS
    REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: STRSS0, STRS3D
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HMTTG, ECREEP, AUXM
    INTEGER :: NROWB2, NWTNITER
    REAL*8  :: RESMAX, RESIDUAL
    !
    !plasticity
    real*8,  allocatable :: gmF(:,:,:), gmG(:,:,:)
    !
    CHARACTER(len=128) :: YNG_IN
    real(8) :: kg,kgphi    ! media geometrica
    real(8) :: rho,rhophi  ! coeficiente da variancia (strenght)
    
    integer :: initStress
    !
    contains
    !
    !**** new **********************
    !
    subroutine bbarmtrx_plast(x, conecnodaiselem, disDirichlet, alhsd, brhsd, idiagd, lmd, tangentMatrix, pressure)
    !
    !.... program to calculate stifness matrix and force array for the
    !        stoke's displacement  element and
    !        assemble into the global left-hand-side matrix
    !        and right-hand side vector
    !
    !function imports
    use mSolverGaussSkyline, only: addrhs, addlhs
    use mfuncoesdeforma, only: shgq, shlq
    use mmalha, only: local
    use mPropGeoFisica, only: rhoCell, calcBiotCoefficient
    
    !variables import
    use mglobaisescalares, only: nrowsh
    use mglobaisescalares, only: ibbar
    use mMalha, only: numnp, nsd, numel, nen
    use mGlobaisArranjos, only:grav
    
    implicit none
    !variables input
    real*8 :: x(nsd,numnp)
    integer :: conecnodaiselem(nen,numel)
    real*8 :: disDirichlet(ndofD, numnp)
    real*8 :: alhsd(nalhsd), brhsd(neqd)
    integer :: idiagd(neqD)
    integer :: lmd(ned2,nen,numel)
    real*8 :: tangentMatrix(nrowb2,nintD, numel)
    real*8 :: pressure(1,numnp)
    
    !variables
    real(8) :: xl(nesd,nen), disL(ned2,nen), pL(1,nen)
    real*8 :: pLInt, densTotal
    real(8) :: elresfd(nee2), eleffmd(nee2,nee2)
    real*8 :: biotCoef
    !
    real(8), external :: rowdot, coldot
    
    logical diag,quad,zerodl,lsym
    !
    !.... local vectors and matrizes
    !
    real(8), dimension(nrowb,nesd)    :: bbarj, bbari, cepbbar
    real(8), dimension(nrowb,nrowb)   :: cep
    !
    real(8) :: shgd(nrowsh,nen,nintd), shld(nrowsh,nen,nintd)
    real(8) :: shgbr(nrowsh,nen)
    real(8) :: detd(nintd), r(nintd), wd(nintd)
    !
    real(8) :: c1
    !
    integer :: nel, i, j, l
    !
    !.... generation of local shape functions and weight values
    !
    call shlq(shld,wd,nintd,nen)
    !
    !      consistent matrix
    !
    diag       = .false.
    
    do nel=1,numel
        bbari    = 0.0d0
        bbarj    = 0.0d0
        cepbbar = 0.0d0
        cep     = 0.0d0
        eleffmd = 0.0d0
        elresfd = 0.0d0
        
        !.... localize coordinates and dirichlet b.c.
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nesd)
        call local(conecNodaisElem(1,nel),disDirichlet,disl,nen,ndofd,ned2)
        call local(conecNodaisElem(1,nel),pressure,pL,nen,1,1)
        
        quad = .true.
        if (conecnodaiselem(3,nel).eq.conecnodaiselem(4,nel)) quad = .false.
        call shgq(xl,detd,shld,shgd,nintd,nel,quad,nen)
        
        !calculate total density and biot parameters
        densTotal = rhoCell(nel)
        biotCoef = calcBiotCoefficient(nel)
        
        !
        !.... calculate mean values of shape function global derivatives
        !.... for mean-dilatational b-bar formulation
        !
        call xMeanSH(shgbr,wd,detd,r,shgd,nen,nintd,iopt,nesd,nrowsh)
        !
        !.... loop over integrationn points
        !
        do l=1,nintd
            !
            !... ... setup tangent matrix cep: order 4x4 for multiplication
            !
            call tang2qx(tangentMatrix(1:16,l,nel),cep)
            !
            c1=detd(l)*wd(l)
            
            !calculate the pressure at the integration point
            pLInt = dot_product(shgD(nrowsh,1:nen,l), pL(1,1:nen))
            
            !
            !.... ...upload b-bar matrix at node j
            !
            do j=1,nen
                call setbb(bbarj,shgd(1:nrowsh,j,l),  &
                    &           shgbr(1:nrowsh,j),r(l),nrowsh,nrowb,iopt,ibbar)
                !
                !.... .... multiply cep*bbarj ===> cepbbar
                !
                cepbbar=matmul(cep,bbarj)
                !
                !.... .... upload b-bar matrix at node i
                !
                do i=1,nen
                    call setbb(bbari,shgd(1:nrowsh,i,l),  &
                        &              shgbr(1:nrowsh,i),r(l),nrowsh,nrowb,iopt,ibbar)
                    !
                    !.... .......  mount elemnt stiffness nodal matrix:
                    !.... .......  k^e_ij= multiply bbar^t_i*(cep*bbar_j)
                    !
                    eleffmd(ned2*i-1,ned2*j-1) = eleffmd(ned2*i-1,ned2*j-1) + coldot(bbari(1:4,1),cepbbar(1:4,1),4)*c1
                    eleffmd(ned2*i-1,ned2*j) = eleffmd(ned2*i-1,ned2*j) + coldot(bbari(1:4,1),cepbbar(1:4,2),4)*c1
                    eleffmd(ned2*i,ned2*j-1) = eleffmd(ned2*i,ned2*j-1) + coldot(bbari(1:4,2),cepbbar(1:4,1),4)*c1
                    eleffmd(ned2*i,ned2*j) = eleffmd(ned2*i,ned2*j) + coldot(bbari(1:4,2),cepbbar(1:4,2),4)*c1
                end do
                
                ! calc pressure contribution
                
                elresfd(ned2*j-1) = elresfd(ned2*j-1) &
                    &                  + biotCoef * pLInt * shgD(1,j,l) * c1 & ! calc pressure contribution
                    &                  + densTotal * grav(1) * shgD(nrowsh,j,l) * c1   ! calc gravity contribution
                elresfd(ned2*j) = elresfd(ned2*j) &
                    &                  + biotCoef * pLInt * shgD(2,j,l) * c1 & ! calc pressure contribution
                    &                  + densTotal * grav(2) * shgD(nrowsh,j,l) * c1     ! calc gravity contribution
            end do
        end do
        !
        !.... computation of dirichlet b.c. contribution
        !
        call ztest(disl,nee2,zerodl)
        !
        if(.not.zerodl)call kdbcgeo(eleffmd,elresfd,disl,nee2,lmd(1,1,nel))
        !
        !.... assemble element stiffness matrix and force array into global
        !....     left-hand-side matrix and right-hand side vector
        !
        lsym=.true.

        call addlhs(alhsd,eleffmd,lmd(1,1,nel),idiagD, nee2, diag,lsym)
        call addrhs(brhsd,elresfd,lmd(1,1,nel),nee2)

    end do

    !
    return
    !
2000 format(5(1pe15.8,2x))
2222 format('elemento (nel)=',i5,2x,' gauss point (l)=',i2/5x   &
        &' c11, c21, c31, c41 =',4(1pe9.2,2x)/5x,                    &
        &' c12, c22, c32, c42 =',4(1pe9.2,2x)/5x,                    &
        &' c13, c23, c33, c43 =',4(1pe9.2,2x)/5x,                    &
        &' c14, c24, c34, c44 =',4(1pe9.2,2x)//)
3333 format('elmnt (nel)=',i5,x,'gauss pt (l)=',i2,2x,4(1pe9.2,2x))
    !
    end subroutine bbarmtrx_plast
    !
    !**** new *********************
    !
    subroutine BBARMTRX_ELAST(x, conecNodaisElem, disDirichlet, alhsd, brhsd, idiagD, lmD, pressure, isUndrained)
    !
    !.... PROGRAM TO CALCULATE STIFNESS MATRIX FOR THE ELASTIC DISPLACEMENT
    !        ELEMENT AND ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX
    !
    use mGlobaisEscalares, only: nrowsh, IBBAR

    use mPropGeoFisica,    only: GEOINDIC
    use mSolverGaussSkyline, only: addrhs, addlhs
    use mPropGeofisica, only: calcMatrixUndrainedProperties, calcMatrixProperties
    use mPropGeoFisica, only: rhoCell, calcBiotCoefficient

    use mFuncoesDeForma,   only: shgq, shlq
    use mMalha,            only: local, numnp, multab
    use mMalha,            only: nsd, numel, nen
    use mGlobaisArranjos, only:grav
    !
    implicit none
    !
    real*8,  intent(in)    :: x(nsd,numnp)
    integer, intent(in)    :: conecNodaisElem(nen,numel)
    real*8 :: disDirichlet(ndofD, numnp)
    real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqD)
    integer, intent(in)    :: idiagD(neqD)
    integer, intent(in)    :: lmD(ned2,nen,numel)
    real*8 :: pressure(1,numnp)
    integer :: isUndrained
    
    !
    real(8) :: xL(nesd,nen), disL(ned2,nen), pL(1,nen)
    real*8 :: pLInt, densTotal
    real*8 :: biotCoef
    real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
    !
    !.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
    !
    LOGICAL DIAG,QUAD,ZERODL,LSYM
    real*8, external :: rowdot, coldot
    !
    !.... LOCAL VECTORS AND MATRIZES
    !
    REAL(8) :: BBARJ(NROWB,NESD), BBARI(NROWB,NESD)
    REAL(8) :: QIXIBBAR(NROWB,NESD)
    INTEGER :: NEL, I, J, L
    REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
    REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
    REAL(8) :: SHGBR(NROWSH,NEN)
    REAL(8) :: CBBAR(NROWB, NROWB)
    REAL(8) :: C1, POISSON, YOUNGMOD
    !
    !
    !.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
    !
    CALL SHLQ(SHLD,WD,NINTD,NEN)
    !
    !.... CONSISTENT MATRIX
    !
    DIAG = .FALSE.

    DO 500 NEL=1,NUMEL
        !
        !.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY
        !
        BBARI = 0.0D0
        BBARJ = 0.0D0
        ELEFFMD=0.D0
        ELRESFD=0.D0
        !
        !.... LOCALIZE COORDINATES AND DIRICHLET B.C.
        !
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),disDirichlet,DISL,NEN,NDOFD,NED2)
        call local(conecNodaisElem(1,nel),pressure,pL,nen,1,1)
        
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
        !
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
        !
        !.... SETUP ELASTIC PARAMETERS FOR DRAINED OR UNDRAINED CONDITIONS
        !
        if (isUndrained == 1) then
            call calcMatrixUndrainedProperties(poisson, youngmod, nel)
        else if (isUndrained == 0) then
            call calcMatrixProperties(poisson, youngmod, nel)
        end if
        !
        !.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD
        !
        CALL SETUPC(CBBAR,YOUNGMOD,POISSON,NROWB,IOPT)
        !
        
        !calculate total density and biot parameters
        densTotal = rhoCell(nel)
        biotCoef = calcBiotCoefficient(nel)
        
        !
        !.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
        !.... .. FOR MEAN-DILATATIONAL B-BAR FORMULATION
        !
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
        !
        !.... .. LOOP OVER INTEGRATIONN POINTS
        !
        DO L=1,NINTD
            !
            !.... ...SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
            !
            C1=DETD(L)*WD(L)
            
            !calculate the pressure at the integration point
            pLInt = dot_product(shgD(nrowsh,1:nen,l), pL(1,1:nen))
            
            !
            !.... ... UPLOAD B-BAR MATRIX AT NODE J
            !
            DO J=1,NEN
                !
                CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L), SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
                !
                !.... .... MULTIPLY CBBAR*BBARJ ===>QIXIBBAR
                !
                CALL MULTAB(CBBAR,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
                !
                !.... .... UPLOAD B-BAR MATRIX AT NODE I
                !
                DO I=1,NEN
                    !
                    CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
                    !
                    !.... .... MOUNT ELEMNT STIFFNESS NODAL MATRIX:
                    !.... ...        K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
                    !
                    ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1) + COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4) * C1
                    !
                    ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J) + COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4) * C1
                    !
                    ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1) + COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4) * C1
                    !
                    ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J) + COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4) * C1
                end do
                
                ! calc pressure contribution
                elresfd(ned2*j-1) = elresfd(ned2*j-1) &
                    &                  + biotCoef * pLInt * shgD(1,j,l) * c1 & ! calc pressure contribution
                &                  + densTotal * grav(1) * shgD(nrowsh,j,l) * c1   ! calc gravity contribution
                elresfd(ned2*j) = elresfd(ned2*j) &
                    &                  + biotCoef * pLInt * shgD(2,j,l) * c1 & ! calc pressure contribution
                &                  + densTotal * grav(2) * shgD(nrowsh,j,l) * c1     ! calc gravity contribution
            end do
            
        end do
        !
        !      COMPUTATION OF DIRICHLET B.C. CONTRIBUTION
        !
        CALL ZTEST(DISL,NEE2,ZERODL)
        !
        IF(.NOT.ZERODL)CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL))
        !
        !.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL
        !        LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR
        !
        LSYM=.TRUE.
        
        CALL ADDLHS(ALHSD,ELEFFMD,LMD(1,1,NEL),IDIAGD,NEE2,DIAG,LSYM)
        CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2)
        !
500 CONTINUE

    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
4500 FORMAT(I8,X,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW FOR 3D GEOMECHANICAL COUPLING *************************************
    !
    SUBROUTINE DRAINCOND(EE,XNU,BSK,ELEMNT,LDRAINED)
    !
    !.... PROGRAM TO SET UP DRAINED/UNDRAINED ELASTICITY MATRIX PARAMETERS
    !
    use mPropGeoFisica, only: BULK, YOUNG, GEOFORM, GEOINDIC
    use mPropGeoFisica, only: ALFABIOT, BSKEMPTON
    !
    IMPLICIT NONE
    !
    LOGICAL :: LDRAINED
    INTEGER :: ELEMNT
    REAL*8  :: NU, BSKEMPT, ALPHARS, ONEM2NU, BEQ4
    REAL*8  :: UNDRAINNU, UNDRAINEE
    REAL*8  :: EE, XNU, BSK
    !
    NU = GEOINDIC('POISSON',GEOFORM(ELEMNT))
    !
    IF (LDRAINED) THEN
        XNU = NU
        EE  = YOUNG(ELEMNT)
        BSK = 0.0D0
        !             bulkrock  = BULK(young(ELEMNT),nu,3.0D0)
        !             write(*,1000) elemnt, bulkrock
    ELSE
        ALPHARS   = ALFABIOT(YOUNG(ELEMNT),NU,ELEMNT)
        BSKEMPT   = BSKEMPTON(YOUNG(ELEMNT),NU,ELEMNT)
        ONEM2NU   = 1.0D0-2.0D0*NU
        BEQ4      = BSKEMPT*ONEM2NU*ALPHARS
        UNDRAINNU = (3.0D0*NU+BEQ4)/(3.0D0-BEQ4)
        UNDRAINEE = YOUNG(ELEMNT)*(1.0D0-2.0D0*UNDRAINNU)/ONEM2NU
        UNDRAINEE = UNDRAINEE/(1.0D0-ALPHARS*BSKEMPT)
        !
        !elemnt,alphars, alphafs, alpharf, bulkrock
        !elemnt,undrainnu, nu, undrainee,young(elemnt)
        !

        !
        !             XNU = NU
        !             EE = YOUNG(ELEMNT)!
        !
        XNU = UNDRAINNU
        EE  = UNDRAINEE
        BSK = BSKEMPT
        !        WRITE(2021,1001) elemnt,alphars, bsk, xnu
    ENDIF
    !
    RETURN
1000 FORMAT('NEL = ',I5,1X,40(1PE15.8,2X))
1001 FORMAT(I5,1X,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !*** NEW *** FOR STOCHASTIC YOUNG MODULUS *******************************
    !
    SUBROUTINE SETUPC(CBBAR,YOUNG,POISSON,NROWB,IOPT)
    !
    !     PROGRAM TO SETUP ELASTICITY TENSOR
    !
    IMPLICIT NONE
    !
    REAL*8  :: YOUNG,POISSON
    INTEGER :: NROWB,IOPT
    REAL(8) :: CBBAR(NROWB, NROWB)
    !
    REAL(8) :: AMU2,ALAM

    !     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION
    !
    !..... SET MATERIAL PARAMETERS FOR OUT-OF-PLANE COMPONENTS
    !
    !
    CBBAR = 0.0
    AMU2 = YOUNG/(1.0D0+POISSON)
    !plane strain
    alam = poisson * young / ((1+poisson)*(1-2*poisson))
    !
    !..... COLUMN MATRIX D3
    !
    CBBAR(1,4) = ALAM
    CBBAR(2,4) = ALAM
    CBBAR(3,4) = 0.0D0
    CBBAR(4,4) = ALAM + AMU2
    !
    !..... TRANSPOSE OF D3
    !
    CBBAR(4,1) = CBBAR(1,4)
    CBBAR(4,2) = CBBAR(2,4)
    CBBAR(4,3) = CBBAR(3,4)
    !
    !..... SET MATERIAL PARAMETERS FOR IN-PLANE COMPONENTS
    !
    !..... .. SET STRESS PLANE CONDITION IOPT.EQ.0
    !
    IF (IOPT.EQ.0) ALAM = ALAM*AMU2/(ALAM + AMU2)
    !
    !..... .. MATRIX D_33: UP TRIANGULAR
    !
    CBBAR(1,1) = ALAM + AMU2
    CBBAR(1,2) = ALAM
    CBBAR(2,2) = CBBAR(1,1)
    CBBAR(1,3) = 0.0D0
    CBBAR(2,3) = 0.0D0
    CBBAR(3,3) = 0.5D0*AMU2
    !
    !..... .. MATRIX D_33: LOW TRIANGULAR
    !
    CBBAR(2,1) = CBBAR(1,2)
    CBBAR(3,1) = CBBAR(1,3)
    CBBAR(3,2) = CBBAR(2,3)
    !
100 CONTINUE
    !
    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
    !
    END SUBROUTINE
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE SETCBBM1(CBBARM1,YOUNG,POISSON,NROWB)
    !
    !.... PROGRAM TO SETUP INVERSE OF ELASTICITY TENSOR
    !....  REFERENCE PAG. 74 M.SADD ELASTICITY THEORY, APPLICATIONS BOOK
    !
    IMPLICIT NONE
    !
    INTEGER :: NROWB
    REAL(8) :: YOUNG, POISSON
    REAL(8) :: CBBARM1(NROWB, NROWB)
    !
    REAL(8) :: EEM1,XNIEEM1, XNIP1EM1
    !
    !... SET MATERIAL PARAMETERS FOR OUT-OF-PLANE COMPONENTS
    !
    CBBARM1  = 0.0D0
    EEM1     = 1.0D0/YOUNG
    XNIEEM1  = -POISSON*EEM1
    XNIP1EM1 = 2.0D0*(1.0D0 + POISSON)*EEM1
    !
    !... COLUMN MATRIX D3
    !
    CBBARM1(1,4) = XNIEEM1
    CBBARM1(2,4) = XNIEEM1
    CBBARM1(3,4) = 0.0D0
    CBBARM1(4,4) = EEM1
    !
    !..... TRANSPOSE OF D3
    !
    CBBARM1(4,1) = CBBARM1(1,4)
    CBBARM1(4,2) = CBBARM1(2,4)
    CBBARM1(4,3) = CBBARM1(3,4)
    !
    !..... SET MATERIAL PARAMETERS FOR IN-PLANE COMPONENTS
    !
    !..... .. MATRIX D_33: UP TRIANGULAR
    !
    CBBARM1(1,1) = EEM1
    CBBARM1(1,2) = XNIEEM1
    CBBARM1(2,2) = CBBARM1(1,1)
    CBBARM1(1,3) = 0.0D0
    CBBARM1(2,3) = 0.0D0
    CBBARM1(3,3) = XNIP1EM1
    !
    !..... .. MATRIX D_33: LOW TRIANGULAR
    !
    CBBARM1(2,1) = CBBARM1(1,2)
    CBBARM1(3,1) = CBBARM1(1,3)
    CBBARM1(3,2) = CBBARM1(2,3)
    !
100 CONTINUE
    !
    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
    !
    END SUBROUTINE SETCBBM1
    !**************************************************************************************
    !**************************************************************************************
    !
    !**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 760 ***************
    !
    SUBROUTINE XMEANSH(SHGBAR,W,DET,R,SHG,NEN,NINT,IOPT,NESD,NROWSH)
    !
    !.... PROGRAM TO CALCULATE MEAN VALUES OF SHAPE FUNCTION
    !        GLOBAL DERIVATIVES FOR B-BAR METHOD
    !
    !        NOTE: IF IOPT.EQ.2, DET(L)=DET(L)*R(L) UPON ENTRY
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: NEN, NROWSH, NINT, IOPT, NESD
    REAL(8) :: SHGBAR(3,NEN),W(*),DET(*),R(*),SHG(NROWSH,NEN,*)
    !
    REAL(8) :: VOLINV, TEMP1, TEMP2
    INTEGER :: L, I, J
    real*8, external :: coldot
    !
    CALL CLEAR(SHGBAR,3*NEN)
    !
    VOLINV = 1.0D0/COLDOT(W,DET,NINT)
    !
    DO 300 L=1,NINT
        TEMP1 = W(L)*DET(L)*VOLINV
        IF (IOPT.EQ.2) TEMP2 = TEMP1/R(L)
        !
        DO 200 J=1,NEN
            DO 100 I=1,NESD
                SHGBAR(I,J) = SHGBAR(I,J) + TEMP1*SHG(I,J,L)
100         CONTINUE
200     CONTINUE
        IF (IOPT.EQ.2) SHGBAR(3,J)=SHGBAR(3,J)+TEMP2*SHG(3,J,L)
300 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !
    !**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 780 ****************
    !
    SUBROUTINE SETBB(BBAR,SHG,SHGBAR,R,NROWSH,NROWB,IOPT,IBBAR)
    !
    !..... PROGRAM TO SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
    !         TWO-DIMENSIONAL CONTINUUM ELEMENTS
    !
    !         IBBAR = 0,  STANDARD B-MATRIX
    !
    !         IBBAR = 1,  MEAN-DILATATIONAL B-MATRIX
    !
    IMPLICIT NONE
    !
    !.... DEACTIVATE ABOVE CARD(S) FOR SINGLE-PRECISION OPERATION
    !
    integer :: NROWSH,NROWB,IOPT,IBBAR
    real*8 ::  SHG(NROWSH), SHGBAR(NROWSH), BBAR(NROWB,2), R
    !
    real*8 :: temp1, temp2, temp3, constb
    !
    BBAR(1,1) = SHG(1)
    BBAR(2,1) = 0.0D0
    BBAR(3,1) = SHG(2)
    BBAR(4,1) = 0.0D0
    !
    BBAR(1,2) = 0.0D0
    BBAR(2,2) = SHG(2)
    BBAR(3,2) = SHG(1)
    BBAR(4,2) = 0.0D0
    !
    !.... AXISYMMETRIC CASE
    !
    IF (IOPT.EQ.2) THEN
        BBAR(4,1) = SHG(3)/R
        BBAR(4,2) = 0.0D0
    ENDIF
    !
    IF (IBBAR.EQ.0) RETURN
    !
    CONSTB = 1.0D0/3.0D0
    !
    !.... ADD CONTRIBUTIONS TO FORM B-BAR
    !
    !.... DEFINE VALUES FOR ROW=4 CASES: PLAIN STRAIN OR AXISYMMETRIC
    !
    IF (IOPT.EQ.2) THEN
        TEMP3 = CONSTB*(SHGBAR(3)-SHG(3)/R)
        BBAR(1,1) = BBAR(1,1)+TEMP3
        BBAR(2,1) = BBAR(2,1)+TEMP3
        BBAR(4,1) = BBAR(4,1)+TEMP3
    ELSE
        BBAR(4,1) = 0.0D0
        BBAR(4,1) = 0.0D0
    ENDIF
    !
    TEMP1 = CONSTB*(SHGBAR(1)-SHG(1))
    TEMP2 = CONSTB*(SHGBAR(2)-SHG(2))
    !
    BBAR(1,1) = BBAR(1,1) + TEMP1
    BBAR(2,1) = BBAR(2,1) + TEMP1
    BBAR(4,1) = BBAR(4,1) + TEMP1
    !
    BBAR(1,2) = BBAR(1,2) + TEMP2
    BBAR(2,2) = BBAR(2,2) + TEMP2
    BBAR(4,2) = BBAR(4,2) + TEMP2
    !
    RETURN
    !
    END SUBROUTINE
    !
    !**** NEW *** FOR GEOMECHANICAL COUPLING ***************************
    !
    SUBROUTINE SETCBBAR(CBBAR,COEF,B,NROWB,NROWSHD,NESD)
    !
    !.... PROGRAM TO SETUP MULTIPLICATION OF ELASTICITY MATRIX "C" WITH BBAR-3D
    !
    IMPLICIT NONE
    !
    INTEGER :: NROWB, NROWSHD, NESD, I, J
    REAL*8, DIMENSION(NROWSHD)    :: COEF
    REAL*8, DIMENSION(NROWB,NESD) :: CBBAR, B
    !
    CBBAR = 0.0D0
    !
    CBBAR(1,1) = COEF(2)*B(1,1) + COEF(3)*(B(2,1)+B(3,1))
    CBBAR(1,2) = COEF(2)*B(1,2) + COEF(3)*(B(2,2)+B(3,2))
    CBBAR(1,3) = COEF(2)*B(1,3) + COEF(3)*(B(2,3)+B(3,3))
    !
    CBBAR(2,1) = COEF(3)*(B(1,1)+B(3,1)) + COEF(2)*B(2,1)
    CBBAR(2,2) = COEF(3)*(B(1,2)+B(3,2)) + COEF(2)*B(2,2)
    CBBAR(2,3) = COEF(3)*(B(1,3)+B(3,3)) + COEF(2)*B(2,3)
    !
    CBBAR(3,1) = COEF(3)*(B(1,1)+B(2,1)) + COEF(2)*B(3,1)
    CBBAR(3,2) = COEF(3)*(B(1,2)+B(2,2)) + COEF(2)*B(3,2)
    CBBAR(3,3) = COEF(3)*(B(1,3)+B(2,3)) + COEF(2)*B(3,3)
    !
    DO 200 I=4,NROWB
        DO 100 J=1,NESD
            CBBAR(I,J) = COEF(4)*B(I,J)
100     CONTINUE
200 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !
    !**** NEW ******************
    !
    SUBROUTINE KDBCGEO(ELEFFM,ELRESF,DL,NEE,LM)
    !
    !.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT
    !     BOUNDARY CONDITION
    !
    USE mGlobaisEscalares
    !
    IMPLICIT NONE
    !
    INTEGER :: NEE
    REAL*8  :: ELEFFM(NEE,*),ELRESF(*),DL(*)
    INTEGER :: LM(*)
    !
    INTEGER :: I,J,L
    REAL(8) :: VAL
    !
    !    THIS VERSION OF KDBC IS ONLY VALID FOR THERMOELASTIC CONS.
    !
    DO 200 J=1,NEE
        L   = LM(J)
        VAL = DL(J)
        IF(L.GT.0)      GO TO 200
        IF(VAL.EQ.0.0D0) GO TO 200
        DO 100 I=1,NEE
            ELRESF(I)=ELRESF(I)-ELEFFM(I,J)*VAL
100     CONTINUE
        !
200 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !**************************************************************************************
    !**************************************************************************************
    subroutine pos4plast(x, conecNodaisElem, u, strainP, prevStrainP, curStress, curStressS, curTrStrainP, curStressTotal, tangentMatrix, biotP, elementIsPlast, pressure, pInit, isUndrained)
    !
    !.... PROGRAM TO UPDATE STRESS FOR NON-LINEAR plasticity MODEL
    !
    !function import
    use mMalha, only: local
    use mFuncoesDeForma,   only: shgq, shlq
    use mPropGeoFisica, only: calcBiotCoefficient, calcMatrixProperties, calcMatrixUndrainedProperties
    use mPropGeoFisica, only: updatePorosity
    
    !variables import
    use mMalha,            only: nsd, numnp, numel, nen
    use mGlobaisEscalares, only: nrowsh
    use mPropGeoFisica,    only: geoform, geoindic
    use mGlobaisEscalares, only: ibbar, tolYield, tolePlas
    
    implicit none
    !variables input
    real(8), intent(in)    :: x(nsd,numnp)
    integer, intent(in)    :: conecnodaiselem(nen,numel)
    real*8 :: u(ndofD, numnp)
    real*8 :: strainP(nrowb, nintD, numel)
    real*8 :: prevStrainP(nrowb, nintD,numel)
    real*8 :: curStress(nrowB, nintD, numel)
    real*8 :: curStressTotal(nrowB, nintD, numel)
    real*8 :: curStressS(nintD, numel)
    real*8 :: curTrStrainP(nintD, numel)
    real*8 :: tangentMatrix(nrowb2,nintD, numel)
    real*8 :: biotP(nrowb,nintD,numel)
    real*8 :: elementIsPlast(numel)
    real*8 :: pressure(1,numnp)
    real*8 :: pInit(1,numnp)
    integer :: isUndrained
    
    !variables
    logical quad
    integer :: j,k,l,nel,nitplas
    !
    !.... input/output vectors and matrizes of subroutine
    !
    real(8), dimension(nesd,nen)      :: xl
    real(8), dimension(ned2,nen)      :: uL
    !
    !.... local vectors and matrizes
    !
    real(8), dimension(nrowb,nesd)    :: bbarj
    real(8), dimension(nrowb,nrowb)   :: qixi, cbbar, cbbarm1
    real(8), dimension(nrowb) :: stressLEff, strainL, qixigrad
    real(8), dimension(nrowb) :: elasticStrain, eptrial, epinit, residuo
    !
    real(8) :: shld(nrowsh,nen,nintd), shgd(nrowsh,nen,nintd)
    real(8) :: shgbr(nrowsh,nen)
    real(8) :: detd(nintd), r(nintd), wd(nintd)

    real(8) :: young, poisson, fyield, fYield1(nrowb), fYield2(nrowb,nrowb)
    real*8 :: hgamma, hnorma
    real*8 :: dpK, dpAlpha
    real*8 :: pressureL(1,nen), pressureInitL(1,nen), pressureIntPoint, pressureInitIntPoint
    real*8 :: biotCoef
    logical :: converged
    real*8 :: identI(nrowb)
    real*8 :: Cep(nrowb,nrowb), NVect(nrowb,1), NxN(nrowb,nrowb)
    real*8 :: qixiGradf(nrowb), gradfQixiGradf, gradfIdentI, denN, tempResult(1,1)
    
    real*8 :: c1
    real*8 :: area
    
    real*8 :: intTrStrainP
    real*8 :: elTrStrain, elTrStrainP
    real*8 :: elP, elPInit

    real*8, external :: rowdot, coldot
    real*8, external :: traceTensor
    !
    !.... etotal      :  deformacao total
    !.... eplas       :  deformacao plastica
    !.... eetrial     :  deformacao elastica teste
    !.... epinit      :  deformacao plastica inicial
    !.... hgamma      :  multiplicador de lagrange na lei de
    !....                de escoamento da deformacao plastica

    !.... subrutinas usadas

    !.... fyield      : calcula a funcao de escoamento
    !.... grads       : gradiente da funcao de escoamento e tensor qixi
    !
    !.... generation of local shape functions and weight values
    !
    !--------------------------
    qixigrad = 0.0d0
    identI = (/ 1.d0, 1.d0, 0.d0, 1.d0 /)
    !      qixip    = 0.0d0
    !
    call shlq(shld,wd,nintd,nen)
    !
    stressLEff = 0.0d0
    !
    do nel=1,numel
        if (isUndrained == 1) then
            call calcMatrixUndrainedProperties(poisson,young, nel)
        else if (isUndrained == 0) then
            call calcMatrixProperties(poisson, young, nel)
        end if
            
        dpK = geoIndic('DPCOHES',geoform(nel))
        dpAlpha = geoIndic('DPALPHA',geoform(nel))
        
        !.... setup stochastic elasticity tensor for bbar method
        call setupc(cbbar,young,poisson,nrowb,iopt)
        
        !..... compute inverse of elastic material tensor cbbarm1
        call setcbbm1(cbbarm1,young ,poisson,nrowb)
        
        !..... ..localize coordinates and dirichlet b.c.
        call local(conecnodaiselem(1,nel),x,xl,nen,nsd,nesd)
        call local(conecnodaiselem(1,nel),u,uL,nen,ndofd,ned2)
        call local(conecnodaiselem(1,nel),pressure,pressureL,nen,1,1)
        call local(conecnodaiselem(1,nel),pInit, pressureInitL,nen,1,1)
        
        !.... ... calculates biot coefficients
        biotCoef = calcBiotCoefficient(nel)
        
        !
        quad = .true.
        if (conecnodaiselem(3,nel).eq.conecnodaiselem(4,nel)) quad = .false.
        !
        call shgq(xl,detd,shld,shgd,nintd,nel,quad,nen)
        !..
        !..... ..setup for axisymmetric option
        !..
        if (iopt.eq.2) then
            do l=1,nintd
                r(l)    = rowdot(shgd(nrowsh,1,l),xl,nrowsh,nesd,nen)
                detd(l) = detd(l)*r(l)
            end do
        endif
        
        !.... .. calculate mean values of shape function global derivatives
        !....    for mean-dilatational b-bar formulation
        call xmeansh(shgbr,wd,detd,r,shgd,nen,nintd,iopt,nesd,nrowsh)
        
        area = 0.d0
        elP = 0.d0
        elPInit = 0.d0
        elTrStrain = 0.d0
        elTrStrainP = 0.d0
        
        !.... ...loop over integration points
        do l=1,nintd
            ! get pressure value at integration point
            pressureIntPoint = dot_product(shgD(nrowsh,1:nen,l), pressureL(1,1:nen))
            pressureInitIntPoint = dot_product(shgD(nrowsh,1:nen,l), pressureInitL(1,1:nen))
            
            !..... ..clear initial strain
            strainL = 0.0d0
            do j=1,nen
                !
                !.... ..... upload b-bar matrix at node j
                call setbb(bbarj,shgd(1:nrowsh,j,l),shgbr(1:nrowsh,j), &
                    &                    r(l),nrowsh,nrowb,iopt,ibbar)
                
                !.... ..... compute strains within intrinsic b-bar formulation
                do k=1,nrowb
                    strainL(k) = strainL(k) + dot_product(bbarj(k,1:2),uL(1:2,j))
                end do
            end do
            
            !material iteration
            !initialize plastic strains
            do k=1,nrowb
                ePTrial(k) = prevStrainP(k,l,nel)
                ePInit(k)  = prevStrainP(k,l,nel)
            end do
            
            call tang2qx(tangentMatrix(1:16,l,nel),qixi)
            
            !.... part 1 from box 4.1 simo-hughes compute predictors
            
            !.... ... deform elast (trial): e_elast=e_total-e_plast
            do k=1,nrowb
                elasticStrain(k) = strainL(k) - eptrial(k)
            end do
            
            !.... ..... compute trial stress
            stressLEff = matmul(cbbar, elasticStrain)
            
            !.... ... compute yield function:
            call dpYieldEfStress(fYield,stressLEff, pressureIntPoint, biotCoef, dpK, dpAlpha)
            
            if (fyield.ge.tolyield) then
                elementIsPlast(nel) = 1.d0
                converged = .false.
                
                hgamma = 0.0d0
                !
                !.... ************* closest point-projection ********************
                !
                do nItPlas = 1, 1000
                    !
                    !.... ...deform elast (trial): e_elast=e_total-e_plast
                    !
                    do k=1,nrowb
                        elasticStrain(k) = strainL(k) - ePTrial(k)
                    end do
                    
                    !.... part 2a de box 4.1 simo-hughes compute residuals
                    
                    !.... ... compute elastic predictor stress: stress = c*e_elastico
                    stressLEff = matmul(cbbar, elasticStrain)
                    
                    !.... ... compute the yield function, its gradient and hessian (drucker-prager):
                    call dpGrads(fYield, fYield1, fYield2, stressLEff, pressureIntPoint, biotCoef, dpK, dpAlpha, .false.)
                    
                    !.... ... residual computation
                    do k=1,nrowb
                        residuo(k) = eptrial(k)-epInit(k)-hGamma*fYield1(k)
                    end do
                    
                    !.... ... compute norm of residuo
                    hNorma = dsqrt(dot_product(residuo,residuo))
                    
                    !.... part 2b de box 4.1 simo-hughes: convergence test:
                    !.... ... convergence test
                    if ((fYield.lt.tolYield).and.(hNormA.lt.tolePlas)) then
                        converged = .true.
                        !exit
                        
                        qixiGradf = matmul(qixi, fyield1)
                        call matMulmTn(fyield1, qixiGradf, tempResult, nrowb, 1, nrowb, 1, 1, 1)
                        gradfQixiGradf = tempResult(1,1)
                        denN = dsqrt(gradfQixiGradf)
                        NVect(1:nrowb,1) = qixiGradf / denN
            
                        NxN = matmul(NVect,transpose(NVect))
                        gradfIdentI = dot_product(fYield1, identI)
                        
                        Cep = qixi - NxN
                        !calc alphaP
                        biotP(1:nrowb,l,nel) = biotCoef*identI(1:nrowb) + (1-biotCoef) * gradfIdentI / gradfQixiGradf * qixiGradf(1:nrowb)

                        exit
                    end if

                    !... compute inverse definition of QIXI (simo pag.175. equation 4.3.15)
                    call calcQixiM1(qixi,cbbarm1,hGamma,fYield2)
                    
                    !.... ... update trial plastic deformation, gamma and qixigrad values
                    call trials(eptrial,hgamma,qixigrad,qixi, residuo, fYield1, fyield,cbbarm1,nrowb)
                end do
                
                if (converged .eqv. .false.) then
                    write(*,*) 'Plastic newton iteration not converged'
                    exit
                end if
            else
                elementIsPlast(nel) = 0.d0
                Cep = cbbar
                biotP(1:nrowb,l,nel) = biotCoef*identI(1:nrowb)
            end if
            
            do k=1,nrowb
                strainP(k,l,nel) = eptrial(k)
                curStress(k,l,nel) = stressLEff(k)
                curStressTotal(k,l,nel) = stressLEff(k) - biotCoef*pressureIntPoint*IdentI(k)
            end do

            intTrStrainP = traceTensor(strainP(:,l,nel),nrowb)
            
            curStressS(l,nel) = traceTensor(curStressTotal(:,l,nel),nrowb)/3.0d0
            curTrStrainP(l,nel) = intTrStrainP
            
            !updatePorosity
            c1 = detd(l)*wd(l)
            area = area + c1
            
            elTrStrain = elTrStrain + traceTensor(strainL,nrowb)*c1
            elTrStrainP = elTrStrainP + intTrStrainP*c1
            elP = elP + pressureIntPoint*c1
            elPInit = elPInit + pressureInitIntPoint*c1
            
            !
            !.... ... update tangent moduli
            !.... .... transfer 4x4-order matrix to global tangent array
            call qx2tang(Cep,tangentMatrix(1:16,l,nel))
        end do
        
        !update porosity
        elTrStrain = elTrStrain / area
        elTrStrainP = elTrStrainP / area
        elP = elP / area
        elPInit = elPInit / area
        
        call updatePorosity(nel, elTrStrain, elTrStrainP, elP, elPInit)
    end do
    
    return
    
    end subroutine pos4plast
    !**************************************************************************************
    !**************************************************************************************


    !
    !**** NEW **** FOR STOCHASTIC AND NON-LINEAR FORMULATION *****************
    !
    subroutine geosetup(numel,nrowb,nintd,iopt,isUndrained, biotP)
    ! function import
    use mPropGeofisica, only: calcMatrixUndrainedProperties, calcMatrixProperties, calcBiotCoefficient
    
    ! variable import
    
    implicit none
    !
    !.... PROGRAM TO SETUP INITIAL INELASTIC TANGENT MATRIX
    !
    !variables input
    integer :: numel,nrowb,nintd,iopt
    integer :: isUndrained
    real*8 :: biotP(nrowb,nintD,numel)
    
    !variables
    real(8) :: cbbar(nrowb, nrowb)
    real(8) :: young
    real(8) :: poisson
    real*8 :: biotCoef
    integer :: nel, l
    real*8 :: identI(nrowb)
    !
    identI = (/ 1.d0, 1.d0, 0.d0, 1.d0 /)

    DO NEL=1,NUMEL
        if (isUndrained == 1) then
            call calcMatrixUndrainedProperties(poisson, young, nel)
        else if (isUndrained == 0) then
            call calcMatrixProperties(poisson, young, nel)
        end if
        
        !
        !.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD
        !
        call setupc(cbbar,young,poisson,nrowb,iopt)
        biotCoef = calcBiotCoefficient(nel)
        !
        !.... SETUP INITIAL TANGENT MATRIX AT ELEMENT GAUSS POINT
        !
        do l=1,nintd
            ! set D
            call qx2tang(cbbar,hmTTG(1:16,l,nel))

            ! set alphaP
            biotP(1:nrowb,l,nel) = biotCoef * identI
        end do
    end do
    !
    RETURN
    !
    END SUBROUTINE
    !
    !**** NEW ******************
    !
    FUNCTION POWERLAW(X,NDERHO)
    !
    USE mPropGeoFisica, only: DTCREEP, CREEPZERO, POWERN, SIGMAREF
    !
    !.....PROGRAM TO COMPUTE POWER LAW LIKE FUCNTION
    !
    IMPLICIT NONE
    !
    REAL*8 :: X
    INTEGER NDERHO
    !
    REAL*8 :: POWERLAW, XIMAGE
    !
    XIMAGE=(DTCREEP*CREEPZERO)*(X/SIGMAREF)**POWERN
    !
    GOTO(100,200) NDERHO
    !
    !.... FUNCTION ONLY
    !
100 CONTINUE
    POWERLAW = XIMAGE
    RETURN
    !..
    !.... FIRST DERIVATIVE OF POWER LAW
    !..
200 CONTINUE
    POWERLAW = POWERN*XIMAGE/X
    RETURN
    !
    END FUNCTION
    !
    !**** NEW ******************
    !
    FUNCTION FNOTLIN(QTRIAL,THREEG,NDERHO,XINPUT)
    !
    !
    !.....PROGRAM TO COMPUTE NON LINEAR FUNCTION FROM POWER LAW
    !
    IMPLICIT NONE
    !
    INTEGER :: NDERHO
    REAL(8) :: QTRIAL, THREEG, X, XINPUT, FNOTLIN
    !
    X = QTRIAL-THREEG*XINPUT
    !
    GOTO(100,200) NDERHO
    !
    !... FUNCTION
    !
100 CONTINUE
    FNOTLIN = XINPUT-POWERLAW(X,NDERHO)
    RETURN
    !
    !... FIRST DERIVATIVE FUNCTION
    !
200 CONTINUE
    FNOTLIN = 1.0D0+POWERLAW(X,NDERHO)*THREEG
    RETURN
    !
    END FUNCTION
    !
    !**** NEW ******************
    !
    FUNCTION DEVNORM2(X,NROWB)
    !
    !.....PROGRAM TO COMPUTE SQUARE NORM OF DEVIATOR TENSOR (2D MODEL)
    !
    IMPLICIT NONE
    !
    INTEGER :: NROWB
    REAL*8  :: X(NROWB)
    !
    REAL*8  :: DEVNORM2
    !
    DEVNORM2 = X(1)*X(1)+X(2)*X(2)+X(4)*X(4)+2.0D0*X(3)*X(3)
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW ******************
    !
    SUBROUTINE COMPTRACE(P,STRS,SIGMAT,NROWB,NUMEL,numelReserv)
    !
    !.... PROGRAM TO COMPUTE TRACE OF TOTAL STRESS
    !.... TOTAL STRESS = SOLID_STRSS_TRACE+FLUID_PRESSURE
    !
    use mMalha,            only: NSD
    use mGlobaisEscalares, only: S3DIM
    use mPropGeoFisica,    only: POISVECT, GRAINBLK
    use mPropGeoFisica,    only: YOUNG, GEOFORM, BULK
    !
    IMPLICIT NONE
    !
    INTEGER :: NEL,NUMEL,NROWB,numelReserv,INDX
    !
    REAL(8), DIMENSION(numelReserv) :: P, SIGMAT
    REAL(8), DIMENSION(NROWB,NUMEL) :: STRS
    !
    REAL(8) :: BULKROCK, BIOTCOEF
    !
    INDX = 6-NSD
    !
    DO 500 NEL=1,numelReserv
        BULKROCK = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
        BIOTCOEF = 1.0D0 - BULKROCK/GRAINBLK(1)
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            SIGMAT(NEL)=STRS(1,NEL)+STRS(2,NEL)+STRS(INDX,NEL) &
                &                  -S3DIM*BIOTCOEF*P(NEL)

        ENDIF
500 CONTINUE


    !
    RETURN
    !
    END SUBROUTINE
    !
    !
    !**** NEW ***** FOR CREEP MODELING *************************************
    !
    SUBROUTINE DEVTENSOR(DEVIATOR,TRACED3,TENSORIN,NROWB)
    !..
    !.... PROGRAM TO COMPUTE DEVIATOR FOR PLANE DEFORMATIONS STATE
    !..
    IMPLICIT NONE
    !
    INTEGER :: NROWB
    REAL*8 :: DEVIATOR(NROWB), TENSORIN(NROWB), TRACED3
    !
    TRACED3 = (TENSORIN(1)+TENSORIN(2)+TENSORIN(4))/3.0D0
    !
    DEVIATOR(1)= TENSORIN(1)-TRACED3
    DEVIATOR(2)= TENSORIN(2)-TRACED3
    DEVIATOR(3)= TENSORIN(3)
    DEVIATOR(4)= TENSORIN(4)-TRACED3
    !
    RETURN
    !
    END SUBROUTINE
    !
    !*** NEW *** FOR STOCHASTIC YOUNG MODULUS *******************************
    !
    SUBROUTINE SETTGMX(QIXI,DEVSTRS,GTIMES2,GTIMES3,BULKROCK,&
        &                   QTRIAL,GAMMA, MECLAW)
    !
    !..... PROGRAM TO SETUP TANGENT MATRIX
    !
    !
    IMPLICIT NONE
    !
    REAL*8  :: GTIMES2, GTIMES3, BULKROCK, QTRIAL, GAMMA
    CHARACTER(5) :: MECLAW
    !
    REAL*8 :: UNITVECT(4),UNITTENS(4,4),DEVPROJ(4,4),QIXI(4,4)
    REAL*8 :: DEVSTRS(4), VECTNORM, X
    REAL*8 :: BLOCOA, BLOCO2G, BLOCO6G2
    INTEGER :: I, J
    !
    UNITVECT(1)=1.0D0
    UNITVECT(2)=1.0D0
    UNITVECT(3)=0.0D0
    UNITVECT(4)=1.0D0
    !
    UNITTENS(1,1)=1.0D0
    UNITTENS(1,2)=0.0D0
    UNITTENS(1,3)=0.0D0
    UNITTENS(1,4)=0.0D0
    UNITTENS(2,1)=0.0D0
    UNITTENS(2,2)=1.0D0
    UNITTENS(2,3)=0.0D0
    UNITTENS(2,4)=0.0D0
    UNITTENS(3,1)=0.0D0
    UNITTENS(3,2)=0.0D0
    UNITTENS(3,3)=0.5D0
    UNITTENS(3,4)=0.0D0
    UNITTENS(4,1)=0.0D0
    UNITTENS(4,2)=0.0D0
    UNITTENS(4,3)=0.0D0
    UNITTENS(4,4)=1.0D0
    !..
    !... MOUNT DEVIATOR PROJETOR
    !..
    DO 20 I=1,4
        DO 10 J=1,4
            DEVPROJ(I,J)=UNITTENS(I,J)-UNITVECT(I)*UNITVECT(J)/3.0D0
10      CONTINUE
20  CONTINUE
    !
    !      WRITE(102,2222) 1,1, ((DEVPROJ(I,J),I=1,4),J=1,4)
    !..
    !.... COMPUTE DEVIATOR NORM
    !
    VECTNORM = DEVNORM2(DEVSTRS,4)
    !..
    !.... COMPUTE FACTORS THAT MULTIPLY FOURTH-ORDER MATRICES
    !..
    X = QTRIAL-GTIMES3*GAMMA
    BLOCOA   = POWERLAW(X,2)/(1.0D0+GTIMES3*POWERLAW(X,2))
    !
    IF (MECLAW.EQ.'CREEP') THEN
        BLOCO2G  = GTIMES2*(1.0D0-GAMMA*GTIMES3/QTRIAL)
        BLOCO6G2 = GTIMES2*GTIMES3*(GAMMA/QTRIAL-BLOCOA)/VECTNORM
    ELSE
        BLOCO6G2 = 0.0D0
        BLOCO2G  = GTIMES2
    ENDIF
    !
    DO 400 I=1,4
        DO 300 J=1,4
            QIXI(I,J)=  BLOCO2G*DEVPROJ(I,J) &
                &               + BULKROCK*UNITVECT(I)*UNITVECT(J) &
                &               + BLOCO6G2*DEVSTRS(I)*DEVSTRS(J)
300     CONTINUE
400 CONTINUE
    !
    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !*** NEW *** FOR STOCHASTIC YOUNG MODULUS *******************************
    !
    SUBROUTINE QX2TANG(QMATR4X4,TANGENT)
    !
    !     PROGRAM TO TRANSFER FROM QIXI 4X4 MATRIX TO TANGENT ARRAY
    !
    IMPLICIT NONE
    !
    !     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION
    !
    REAL*8 :: TANGENT(16),QMATR4X4(4,4)
    !
    TANGENT(1)  = QMATR4X4(1,1)
    TANGENT(2)  = QMATR4X4(1,2)
    TANGENT(3)  = QMATR4X4(1,3)
    TANGENT(4)  = QMATR4X4(1,4)
    TANGENT(5)  = QMATR4X4(2,1)
    TANGENT(6)  = QMATR4X4(2,2)
    TANGENT(7)  = QMATR4X4(2,3)
    TANGENT(8)  = QMATR4X4(2,4)
    TANGENT(9)  = QMATR4X4(3,1)
    TANGENT(10) = QMATR4X4(3,2)
    TANGENT(11) = QMATR4X4(3,3)
    TANGENT(12) = QMATR4X4(3,4)
    TANGENT(13) = QMATR4X4(4,1)
    TANGENT(14) = QMATR4X4(4,2)
    TANGENT(15) = QMATR4X4(4,3)
    TANGENT(16) = QMATR4X4(4,4)
    !
    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW *** FOR PLASTICITY *********************************************
    !
    SUBROUTINE UPDATEINCR(IDDIS,BRHSD,NDOFD,NUMNP)
    !
    !.... UPDATE TRIAL DISPLACEMENT ARRAY
    !
    IMPLICIT NONE
    !
    !.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
    !
    INTEGER :: NDOFD,NUMNP
    INTEGER :: IDDIS(NDOFD,*)
    REAL*8  :: BRHSD(*)
    !
    INTEGER :: I, J, K
    !
    DO 200 I=1,NDOFD
        !
        DO 100 J=1,NUMNP
            K = IDDIS(I,J)
            IF (K.GT.0) then
                DTRL(I,J) = DTRL(I,J) + BRHSD(K)
            ENDIF
100     CONTINUE
        !
200 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !
    !*** NEW *** FOR STOCHASTIC YOUNG MODULUS *******************************
    !
    SUBROUTINE TANG2QX(TANGENT,QMATR4X4)
    !
    !     PROGRAM TO TRANSFER FROM TANGENT ARRAY TO QIXI 4X4 MATRIX
    !
    IMPLICIT NONE
    !
    !     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION
    !
    real*8 :: TANGENT(16),QMATR4X4(4,4)
    !
    QMATR4X4(1,1) = TANGENT(1)
    QMATR4X4(1,2) = TANGENT(2)
    QMATR4X4(1,3) = TANGENT(3)
    QMATR4X4(1,4) = TANGENT(4)
    QMATR4X4(2,1) = TANGENT(5)
    QMATR4X4(2,2) = TANGENT(6)
    QMATR4X4(2,3) = TANGENT(7)
    QMATR4X4(2,4) = TANGENT(8)
    QMATR4X4(3,1) = TANGENT(9)
    QMATR4X4(3,2) = TANGENT(10)
    QMATR4X4(3,3) = TANGENT(11)
    QMATR4X4(3,4) = TANGENT(12)
    QMATR4X4(4,1) = TANGENT(13)
    QMATR4X4(4,2) = TANGENT(14)
    QMATR4X4(4,3) = TANGENT(15)
    QMATR4X4(4,4) = TANGENT(16)
    !
    RETURN
    !
2000 FORMAT(5(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW **** FOR INCOMPRESSIBILITY ***************************************
    !
    SUBROUTINE PRINT_DX(DIS ,PORE  ,P     ,S     ,MASCN , &
        &                    VC  ,AVSTRS,NDOF2 ,NUMEL ,NROWB ,NUMNP)
    !
    USE mGlobaisEscalares, only: NUMDX, NNP
    use mLeituraEscritaSimHidroGeoMec,   only: PRINT_DXINFO, IFEDX, PATHDX
    use mPropGeoFisica,    only: GEOFORM, GEOINDIC
    !
    !.... PROGRAM TO PRINT DATA FROM GEOMECHANIC MODEL
    !
    IMPLICIT NONE
    !
    !.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
    !
    CHARACTER*30 NIDISP,NIPRSR,NIPORE,NICREP,NISATR,NIVELT
    CHARACTER*30 NISIGX,NISIGY,NISGTA,NISIGZ
    CHARACTER*30 NIS2S1,NIMASC
    !
    CHARACTER*3 ASTEP
    !
    !.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
    !
    INTEGER :: I,J,K,NEL,NDOF2,NUMEL,NROWB,NUMNP
    INTEGER :: IDISP, IPRSR, IPORE, ICREP, ISATR, IVELT
    INTEGER :: IMASC, ISIGX, ISIGY, ISGTA, ISIGZ, IS2S1
    !
    REAL(8), DIMENSION(NDOF2,NUMNP) :: DIS
    REAL(8), DIMENSION(NDOF2,*)     :: VC
    REAL(8), DIMENSION(NUMEL)       :: PORE,P,S,MASCN
    REAL(8), DIMENSION(NROWB,NUMEL) :: AVSTRS
    !
    REAL(8), DIMENSION(NROWB)       :: TENSAO, DEVSTRS
    !
    REAL(8) :: ROOT3D2,TRCSTRS,QTRIAL,STRIAL,SPORE
    REAL(8) :: PHIDRO, ROOT, S1, S2
    !
    DATA ROOT3D2/1.224744871391589D0/
    !
    !... OUTPUT DATA FILES
    !
    IF (NUMDX.EQ.0) RETURN
    !
    IDISP = 631
    IPRSR = 632
    IPORE = 633
    ICREP = 634
    ISATR = 635
    IVELT = 636
    IMASC = 637
    ISIGX = 638
    ISIGY = 639
    ISGTA = 640
    ISIGZ = 641
    IS2S1 = 642
    !
    !.... SETUP FILES COUNTER
    !
    !      WRITE(ASTEP,'(I3.3)') IDINT(TPRT_PHI/DTPRT_PHI)
    WRITE(ASTEP,'(I3.3)') NNP/NUMDX
    !
    !.... OUT-PUT FILES FOR OPEN-DX VISUALIZATION
    !.... DATA VALUES AT NODAL FEM POINTS
    !
    NIDISP = PATHDX//'/disp'//ASTEP//'.stoc'
    NIPRSR = PATHDX//'/prsr'//ASTEP//'.stoc'
    NIPORE = PATHDX//'/pore'//ASTEP//'.stoc'
    NICREP = PATHDX//'/crep'//ASTEP//'.stoc'
    NISATR = PATHDX//'/satr'//ASTEP//'.stoc'
    NIVELT = PATHDX//'/velt'//ASTEP//'.stoc'
    NIMASC = PATHDX//'/masc'//ASTEP//'.stoc'
    NISIGX = PATHDX//'/sigx'//ASTEP//'.stoc'
    NISIGY = PATHDX//'/sigy'//ASTEP//'.stoc'
    NISGTA = PATHDX//'/sigt'//ASTEP//'.stoc'
    NISIGZ = PATHDX//'/sigz'//ASTEP//'.stoc'
    NIS2S1 = PATHDX//'/s2s1'//ASTEP//'.stoc'
    !
    !.....OPEN NODAL DATA FILES
    !
    OPEN(UNIT=IDISP, FILE= NIDISP)
    OPEN(UNIT=IPRSR, FILE= NIPRSR)
    OPEN(UNIT=IPORE, FILE= NIPORE)
    OPEN(UNIT=ICREP, FILE= NICREP)
    OPEN(UNIT=ISATR, FILE= NISATR)
    OPEN(UNIT=IVELT, FILE= NIVELT)
    OPEN(UNIT=IMASC, FILE= NIMASC)
    OPEN(UNIT=ISIGX, FILE= NISIGX)
    OPEN(UNIT=ISIGY, FILE= NISIGY)
    OPEN(UNIT=ISGTA, FILE= NISGTA)
    OPEN(UNIT=ISIGZ, FILE= NISIGZ)
    OPEN(UNIT=IS2S1, FILE= NIS2S1)
    !
    !.... PRINT DISPLACEMENTS
    !
    DO 30 I=1,NUMNP
        WRITE(IDISP,4000) (DIS(J,I),J=1,NDOF2)
30  CONTINUE
    !
    DO 500 NEL=1,NUMEL
        !
        !.... PRINT PRESSURE
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') SPORE = P(NEL)
        !      IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -10.0D0
        IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -GEOPRSR(NEL)
        !      WRITE(IPRSR,4000) XC(1,NEL),XC(2,NEL),SPORE
        WRITE(IPRSR,4000) SPORE
        !
        !.... PRINT GEOMECHANIC POROSITY: "PORE"
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') SPORE = PORE(NEL)
        IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -10.0D0
        WRITE(IPORE,4000) SPORE
        !
        !.... PRINT STRESS CREEP VALUE:
        !
        DO 100 K=1,NROWB
            !         TENSAO(KK) = STRSS(NEL,1,KK)
            TENSAO(K) = AVSTRS(K,NEL)  ! +STRSS0(K,NEL)
100     CONTINUE
        !
        !      IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
        !            QTRIAL = -10.0D0
        !         ELSE
        !.... ... .. COMPUTE DEVIATOR STRESS TENSOR
        CALL DEVTENSOR(DEVSTRS,TRCSTRS,TENSAO,NROWB)
        !.... ... .. COMPUTE EFFECTIVE TRIAL STRESS
        QTRIAL=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
        !      ENDIF
        !
        WRITE(ICREP,4000) QTRIAL
        !
        !.... PRINT SATURATION
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') STRIAL =  S(NEL)
        IF (GEOFORM(NEL).NE.'RESERVATORIO') STRIAL = -10.0D0
        WRITE(ISATR,3000) STRIAL
        !
        !.... PRINT MASS CONTENT
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') STRIAL = MASCN(NEL)
        IF (GEOFORM(NEL).NE.'RESERVATORIO') STRIAL = -10.0D0
        WRITE(IMASC,4000) STRIAL
        !
        !.... PRINT TOTAL VELOCITY
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            STRIAL = VC(1,NEL)
            SPORE  = VC(2,NEL)
        ELSE
            STRIAL = 0.0D0
            SPORE  = 0.0D0
        ENDIF
        !       WRITE(IVELT,4000) (VC(J,NEL),J=1,NDOF2)
        WRITE(IVELT,4000)  STRIAL, SPORE
        !
        !.... PRINT STRESS_X COMPONENT
        !
        WRITE(ISIGX,4000) TENSAO(1)
        !
        !.... PRINT STRESS_Y COMPONENT
        !
        WRITE(ISIGY,4000) TENSAO(2)
        !
        !.... PRINT STRESS_XY COMPONENT
        !
        WRITE(ISGTA,4000) TENSAO(3)
        !      WRITE(ISGTA,4000) (TENSAO(1)+tensao(2))*GEOINDIC('POISSON',GEOFORM(NEL)), tensao(4)
        !.... COMPUTE PRINCIPAL STRESS
        !
        PHIDRO = 0.5D0*(TENSAO(1)+TENSAO(2))
        ROOT=DSQRT((0.5D0*(TENSAO(1)-TENSAO(2)))**2+TENSAO(3)**2)
        S1=PHIDRO+ROOT
        S2=PHIDRO-ROOT
        !
        !      WRITE(ISGTA,4000) S1
        !
        !      S3=GEOINDIC('POISSON',GEOFORM(NEL))*(s1+s2)
        !
        !      TAU = (S1-S2)**2+(S2-S3)**2+(S3-S1)**2
        !      TAU = DSQRT(TAU)/3.0D0
        !      TAU = ROOT
        !      SIGMAN = DABS(S1+S2+S3)/3.0D0
        !      SIGMAN = DABS(PHIDRO)

        !      WRITE(ISGTA,4000) s3- tensao(4)
        !.... PRINT STRESS_Z COMPONENT
        !
        WRITE(ISIGZ,4000) TENSAO(4)
        !
        !.... DILATANCY EXPERIMENT
        !
        !      IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
        !              STRIAL = 0.0D0
        !           ELSE
        !              STRIAL = TAU-0.8996*SIGMAN+0.01697*SIGMAN**2
        !      ENDIF
        !
        !      WRITE(ISIGZ,4000) STRIAL
        !
        !.... PRINT YOUNG MODULUS
        !
        !      IF (NNP.EQ.0) S1=1.0D0

        WRITE(IS2S1,4000) s2  !/s1 ! tensao(2)-spore
        !
500 CONTINUE
    
    CLOSE(IDISP)
    CLOSE(IPRSR)
    CLOSE(IPORE)
    CLOSE(ICREP)
    CLOSE(ISATR)
    CLOSE(IMASC)
    CLOSE(ISIGX)
    CLOSE(ISIGY)
    CLOSE(ISIGZ)
    CLOSE(IS2S1)
    CLOSE(IVELT)
    !
    !.... PRINT INFORMATION ON NODAL DX FILE "nodestoc.dx"
    !
    CALL PRINT_DXINFO('WRITE_FEDX_DATA',IFEDX,NUMNP,NUMEL)
    !
    RETURN
    !
2000 FORMAT(I5,2X,40(1PE15.8,2X))
3000 FORMAT(2X,5(F25.15,2x))
4000 FORMAT(2X,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW **** FORCES APPLIED ON IRREGULAR MESH ******************************
    !
    SUBROUTINE InSeaLoad(FDIS,NDOF2,NSD,NUMNP,NLVECT,XTERLOAD)
    !
    !**** PROGRAM TO MULTIPLY VERTICAL NODAL FORCES ON TOP POST-SALT BY SEALOAD
    !
    IMPLICIT NONE
    !
    INTEGER  :: NSD, NODE, NDOF2, NUMNP, NLVECT
    REAL(8), DIMENSION(NDOF2,NUMNP,NLVECT) :: FDIS
    REAL(8) :: XTERLOAD
    !
    IF (NSD.EQ.2) THEN
        DO 120 NODE=1,NUMNP
            !         IF (FDIS(1,NODE,1).NE.0.0D0) FDIS(1,NODE,1) = XTERLOAD*FDIS(1,NODE,1)
            IF (FDIS(2,NODE,1).NE.0.0D0) FDIS(2,NODE,1) = XTERLOAD*FDIS(2,NODE,1)
120     CONTINUE
    ENDIF
    !
    IF (NSD.EQ.3) THEN
        DO 130 NODE=1,NUMNP
            !         IF (FDIS(1,NODE,1).NE.0.0D0) FDIS(1,NODE,1) = XTERLOAD*FDIS(1,NODE,1)
            !         IF (FDIS(2,NODE,1).NE.0.0D0) FDIS(2,NODE,1) = XTERLOAD*FDIS(2,NODE,1)
            IF (FDIS(3,NODE,1).NE.0.0D0) FDIS(3,NODE,1) = XTERLOAD*FDIS(3,NODE,1)
130     CONTINUE
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE
    !
    !**** NEW ******************
    !
    SUBROUTINE PRSRINIT(GEOPRSR,P,NUMEL,numelReserv,LDRAINED)
    !
    !.... PROGRAM TO SETUP INITIAL HIDROSTATIC PRESSURE
    !
    use mLeituraEscritaSimHidroGeoMec,   only: SOLIDONLY
    use mGlobaisEscalares, only: TypeProcess
    use mMalha,            only: LEFTLINE, RGHTLINE
    use mPropGeoFisica,    only: GEOFORM, FNCMECLAW
    use mPropGeoFisica,    only: XTERLOAD
    !
    use mPropGeoFisica,    only: BULK
    use mPropGeoFisica,    only: GEOFORM, GEOINDIC
    !
    IMPLICIT NONE
    !
    INTEGER ::  NEL, NUMEL, NUMELRESERV
    !
    LOGICAl :: LDRAINED
    REAL(8), DIMENSION(NUMEL)         :: GEOPRSR
    REAL(8), DIMENSION(1,numelReserv) :: P
    REAL*8  :: THREEA, BSKEM, UNDRAINU, PLOAD
    !
    IF (LDRAINED) THEN
        DO 100 NEL=1,numel
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = GEOPRSR(NEL)
            !.... ...projection of reservoir pressure to boundary condition for test
100     CONTINUE
        !         nelFirst = nelXReserv-nelYReserv+1
        !         nelStep  = nelXReserv*nelYReserv
        !         write(2020,4500) numelReserv,nelXReserv
        !         DO 150 NEL=nelXReserv,numelReserv,nelXReserv
        !            k = nel/nelXReserv
        !            PWELL(NEL/nelXReserv) = P(1,NEL)
        !             write(2020,4500) nel, K, (xc(i,nel),i=1,nsd), pwell(K)
        !             write(2020,1000) (xc(i,nel),i=1,nsd)  !, pwell(K)
        !150   CONTINUE
        !
        !         DO 150 NEL=nelFirst,numelReserv,nelStep
        !            K = nel/nelstep + 1
        !            PWELL(K) = P(1,NEL)
        !150      CONTINUE
        RETURN
    ENDIF
    !
    IF (TRIM(TypeProcess).EQ.'TERZAGHI') THEN
        DO 200 NEL=1,numel
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = XTERLOAD
200     CONTINUE
        RETURN
    ENDIF
    !
    IF (TRIM(TypeProcess).EQ.'MANDEL') THEN
        THREEA   = 3.0D0*(RGHTLINE-LEFTLINE)
        BSKEM    = 1.0D0
        UNDRAINU = 0.5D0
        PLOAD    = XTERLOAD*BSKEM*(1.0D0+UNDRAINU)/THREEA
        DO 210 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = PLOAD
210     CONTINUE
        RETURN
    ENDIF
    !
    IF (SOLIDONLY) THEN
        DO 220 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = GEOPRSR(NEL)
220     CONTINUE
        RETURN
    ENDIF
    !
    RETURN
1000 FORMAT(2X,40(1PE15.8,2X))
4500 FORMAT(I8,X,i4,x,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW ******************
    !
    SUBROUTINE PRSRINIT3D(GEOPRSR,P, NUMEL,numelReserv,LDRAINED)
    !
    !.... PROGRAM TO SETUP INITIAL HIDROSTATIC PRESSURE
    !
    use mGlobaisEscalares, only: TypeProcess
    use mLeituraEscritaSimHidroGeoMec, only: SOLIDONLY
    use mMalha,            only: LEFTLINE, RGHTLINE
    use mPropGeoFisica,    only: GEOFORM, FNCMECLAW
    use mPropGeoFisica,    only: XTERLOAD
    !
    use mPropGeoFisica,    only: BULK
    use mPropGeoFisica,    only: GEOFORM, GEOINDIC
    !
    IMPLICIT NONE
    !
    INTEGER ::  NEL, NUMEL, NUMELRESERV
    !
    LOGICAl :: LDRAINED
    REAL(8), DIMENSION(NUMEL)         :: GEOPRSR
    REAL(8), DIMENSION(1,numelReserv) :: P
    REAL*8  :: THREEA, BSKEM, UNDRAINU, PLOAD
    !
    IF (LDRAINED) THEN
        DO 100 NEL=1,numel
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') then
                P(1,NEL) = GEOPRSR(NEL)
            endif
            !.... ...projection of reservoir pressure to boundary condition for test
100     CONTINUE

        RETURN
    ENDIF
    !
    IF (TRIM(TypeProcess).EQ.'TERZAGHI') THEN
        DO 200 NEL=1,numel
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') then
                P(1,NEL) = XTERLOAD
            endif
200     CONTINUE
        RETURN
    ENDIF
    !
    IF (TRIM(TypeProcess).EQ.'MANDEL') THEN
        THREEA   = 3.0D0*(RGHTLINE-LEFTLINE)
        BSKEM    = 1.0D0
        UNDRAINU = 0.5D0
        PLOAD    = XTERLOAD*BSKEM*(1.0D0+UNDRAINU)/THREEA
        DO 210 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = PLOAD
210     CONTINUE
        RETURN
    ENDIF
    !
    IF (SOLIDONLY) THEN
        DO 220 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = GEOPRSR(NEL)
220     CONTINUE
        RETURN
    ENDIF
    !
    RETURN
1000 FORMAT(2X,40(1PE15.8,2X))
4500 FORMAT(I8,X,i4,x,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION HIDROSTAT(Z1, Z2, RHOEQ, BULKEQ, G, TASK)
    !
    !.... PROGRAM TO COMPUTE HIDROSTATIC PRESSURE FIELD
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z1, Z2, HIDROSTAT, RHOEQ, BULKEQ, G
    CHARACTER(6) TASK
    !                  123456
    IF (TASK.EQ.'INCOMP') THEN
        HIDROSTAT = PRSRINCOMP(Z1,Z2,RHOEQ,G)
        RETURN
    ENDIF
    !
    IF (TASK.EQ.'EXACTL') THEN
        HIDROSTAT = PRSREXACT(Z1,Z2,RHOEQ,BULKEQ,G)
        RETURN
    ENDIF
    !
    IF (TASK.EQ.'LINEAR') THEN
        HIDROSTAT = PRSRLINEAR(Z1,Z2,RHOEQ,BULKEQ,G)
        RETURN
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION PRSRINCOMP(Z1, Z2, RHOEQ, G)
    !
    !.... PROGRAM TO COMPUTE PRESSURE OF INCOMPRESSIBLE FLUID
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z1, Z2, PRSRINCOMP, RHOEQ, G
    !
    PRSRINCOMP = RHOEQ*G*(Z1-Z2)
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION PRSREXACT(Z1, Z2, RHOEQ, BULKEQ, G)
    !
    !.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: EXACT FORMULATION
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z1, Z2, PRSREXACT, RHOEQ, BULKEQ, G
    !
    PRSREXACT = PRSRLOG(Z1,RHOEQ,BULKEQ,G)-PRSRLOG(Z2,RHOEQ,BULKEQ,G)
    !
    RETURN
    !
    END FUNCTION
    !
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION PRSRLOG(Z, RHOEQ, BULKEQ, G)
    !
    !.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: EXACT FORMULATION
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z, PRSRLOG, RHOEQ, BULKEQ, G
    !
    PRSRLOG = -(DLOG(1.0D0-RHOEQ*G*BULKEQ*Z))/BULKEQ
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION PRSRLINEAR(Z1, Z2, RHOEQ, BULKEQ, G)
    !
    !.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: LINEAR APROXIMATION
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z1, Z2, PRSRLINEAR, RHOEQ, BULKEQ, G
    !
    PRSRLINEAR = PRSREXP(Z1,RHOEQ,BULKEQ,G)-PRSREXP(Z2,RHOEQ,BULKEQ,G)
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW **** FOR INITIAL STRESS *****************************************
    !
    FUNCTION PRSREXP(Z, RHOEQ, BULKEQ, G)
    !
    !.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: LINEAR APROXIMATION
    !
    IMPLICIT NONE
    !
    REAL(8) :: Z, PRSREXP, RHOEQ, BULKEQ, G
    !
    PRSREXP = (DEXP((RHOEQ*BULKEQ)*(G*Z))-1.0D0)/BULKEQ
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW *****************
    !
    SUBROUTINE POS4_MASSCNT(DIVU,DIVU0, P, P0, &
        &         PORE,PORE0,YOUNG,MASCN,MASCN0,PHIEULER,NUMEL,NUMELRESERV)
    !
    !.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
    !
    use mLeituraEscritaSimHidroGeoMec,   only: CODERROR
    use mGlobaisEscalares, only: S3DIM
    use mPropGeoFisica,    only: POISVECT, GRAINBLK, BULKWATER
    use mPropGeoFisica,    only: GEOFORM, BULK
    !
    IMPLICIT NONE
    !
    INTEGER :: NEL, NUMEL, numelReserv
    !
    !.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
    !
    REAL*8,  INTENT(in)             :: P(1,numelReserv), P0(numelReserv)
    REAL(8), DIMENSION(numelReserv) :: PORE, PORE0, MASCN, MASCN0, PHIEULER
    REAL(8), DIMENSION(NUMEL)       :: DIVU, DIVU0, YOUNG
    !
    !.... ..LOCAL VECTORS AND MATRIZES
    !
    REAL(8) :: DIFFDIVU, DIFFPRES
    !
    REAL(8) :: BULKROCK, BIOTCOEF
    REAL(8) :: DEFNM1, DEFMM1, JACOBIAN
    !
    DO 500 NEL=1,NUMELRESERV
        !
        IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            !
            BULKROCK  = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
            BIOTCOEF  = 1.0D0 - BULKROCK/GRAINBLK(1)   !also called ALPHA
            !
            !... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35
            !
            DEFNM1    = (BIOTCOEF-PORE0(NEL))/GRAINBLK(1)
            !
            !... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
            !
            DEFMM1 = DEFNM1 + PORE0(NEL)/BULKWATER
            !
            !... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27
            !
            JACOBIAN  = 1.0D0 + DIVU(NEL)
            !
            DIFFDIVU  = DIVU(NEL)-DIVU0(NEL)
            DIFFPRES  = P(1,NEL)-P0(NEL)
            !
            ! ..COMPUTE LAGRANGIAN POROSITY
            !           PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-diffdivu)
            !
            !new line: linearized See: Coussy 2.Edt. Eqs. 4.19
            PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES

            !
            !...COMPUTE EULERAIN POROSITY
            PHIEULER(NEL) = PORE(NEL)/JACOBIAN
            !
            !...COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
            !
            MASCN(NEL)= MASCN0(NEL)+(BIOTCOEF*DIFFDIVU+DEFMM1*DIFFPRES)

            IF ((MASCN(NEL).LT.0.0D0).OR.(PORE(NEL).LT.0.0D0)) THEN
                write(2020,*) 'nel=',nel
                write(2020,*) 'alfa=',BIOTCOEF
                write(2020,*) 'divu=',divu(nel)
                write(2020,*) 'divu0=',divu0(nel)
                write(2020,*) 'biotmod=',defmm1
                write(2020,*) 'diffpress=',diffpres
                write(2020,*) '1/m(p-p0)=',defmm1*DIFFPRES
                write(2020,*) 'masscont=',mascn(nel), mascn0(nel)
                !                               123456789+12345678
                CALL CODERROR(4,'ALSO SEE fort.2020')
            ENDIF
            !
        ENDIF
        !
500 CONTINUE
    !
    RETURN
    !
4000 FORMAT(2X,40(1PE15.8,2X))
    !
    END SUBROUTINE
    !
    !**** NEW **** FOR VISCOELASTICITY ***************************************
    !
    SUBROUTINE POS4STRS(x, conecNodaisElem, STRESS, DIVU)
    !
    use mLeituraEscritaSimHidroGeoMec,   only: CODERROR
    use mMalha,            only: nen, LOCAL
    use mMalha,            only: nsd, numnp, numel
    use mGlobaisEscalares, only: nrowsh
    use mGlobaisEscalares, only: ibbar
    use mFuncoesDeForma,   only: shgq, shlq
    use mPropGeoFisica,    only: YOUNG
    use mPropGeoFisica,    only: GEOFORM, GEOINDIC
    !
    !.... PROGRAM TO COMPUTE SOLID ELEMENT: STRESS AND VOLUMETRIC DEFORMATION
    !....  ......  .....    COMPUTED FROM CREEP MODEL
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in)             :: conecNodaisElem(NEN,NUMEL)
    REAL(8), intent(in)             :: X(NSD,NUMNP)
    
    REAL(8), DIMENSION(NROWB,NUMEL) :: STRESS
    REAL(8), DIMENSION(NUMEL)       :: DIVU
    real*8, external :: rowdot, coldot
    !
    LOGICAL QUAD
    REAL(8), DIMENSION(NESD,NEN)         :: XL
    REAL(8), DIMENSION(NED2,NEN)         :: DISL
    REAL(8), DIMENSION(NINTD)            :: WD, DETD, R
    REAL(8), DIMENSION(NROWSH,NEN,NINTD) :: SHLD,SHGD
    REAL(8), DIMENSION(NROWSH,NEN)       :: SHGBR
    REAL(8), DIMENSION(NROWB,NROWB)      :: CBBAR
    !
    !.... LOCAL VECTORS AND MATRIZES
    !
    REAL(8), DIMENSION(NROWB,NESD)  :: BBARJ
    REAL(8), DIMENSION(NROWB)       :: STRAIN, DEVSTRS, TENSAO
    REAL(8), DIMENSION(4)           :: UNITVEC
    REAL(8), DIMENSION(NUMEL)       :: DIVULOC
    !
    REAL*8  :: POISSON, AREA, C1, ROOT3D2, TRCSTRS, QVM
    INTEGER :: NEL, L, J, K
    !
    UNITVEC(1) = 1.0D0
    UNITVEC(2) = 1.0D0
    UNITVEC(3) = 0.0D0
    UNITVEC(4) = 1.0D0
    ROOT3D2    = 1.224744871391589D0
    !
    !.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
    !
    CALL SHLQ(SHLD,WD,NINTD,NEN)
    !
    DO 500 NEL=1,NUMEL
        !
        POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
        !
        !.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD
        !
        CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
        !
        !...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
        !
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL), DIS, DISL, NEN, NDOFD, NED2)
        !
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
        !
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
        !
        DIVULOC(NEL) = 0.0D0
        !..
        !.... ...CLEAR INITIAL STRAIN
        !..
        CALL CLEAR(STRAIN,NROWB)
        CALL CLEAR(TENSAO,NROWB)
        !..
        !.... ...DEFINE ELEMENT AREA
        !..
        AREA = 0.0D0
        !..
        !.... ...SETUP FOR AXISYMMETRIC OPTION
        !..
        IF (IOPT.EQ.2) THEN
            DO 10 L=1,NINTD
                R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
                DETD(L) = DETD(L)*R(L)
10          CONTINUE
        ENDIF
        !
        !.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
        !.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
        !
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
        !
        !.... ..LOOP OVER INTEGRATIONN POINTS
        !
        DO 300 L=1,NINTD
            !
            C1=WD(L)*DETD(L)
            AREA = AREA + C1
            !
            DO 200 J=1,NEN
                !
                !.... ..UPLOAD B-BAR MATRIX AT NODE J
                !
                CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), R(L),NROWSH,NROWB,IOPT,IBBAR)
                !
                !.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION
                !
                DO 200 K=1,NROWB
                    STRAIN(K)=STRAIN(K) + COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
200         CONTINUE
            !..
            !.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION
            !..
            !.... ..TO COMPUTE MEAN VOLUMETRIC CREEP
            !.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
            !..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
            !
            !.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
            !
            QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
            !
            DO 220 K=1,NROWB
                ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+                    &
                    &                     ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
220         CONTINUE
            !
            DO 250 K=1,NROWB
                TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
250         CONTINUE
            !
300     CONTINUE
        !
        DO 350 K=1,NROWB
            STRAIN(K)=STRAIN(K)/AREA
            STRESS(K,NEL)=TENSAO(K)/AREA! - STRSS0(K,NEL)
350     CONTINUE
        !
        !.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT
        !
        DIVULOC(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
        !
500 CONTINUE
    !..
    DO 700 NEL=1,NUMEL
        DIVU(NEL) = DIVULOC(NEL)
700 CONTINUE
    !
    RETURN
    !
4000 FORMAT(2X,40(1PE15.8,2X))
    !
    END SUBROUTINE POS4STRS
    !**************************************************************************************
    !**************************************************************************************
    subroutine pos4inc(x, conecNodaisElem, u, stress, stressS, stressTotal, p, pInit, isUndrained)
    !
    !.... Program to compute the strain, stress and internal force given an incremental displacement
    !
    
    !function imports
    use mLeituraEscritaSimHidroGeoMec,   only: CODERROR
    use mMalha,            only: nen, LOCAL
    use mMalha,            only: nsd, numnp, numel
    use mGlobaisEscalares, only: nrowsh
    use mGlobaisEscalares, only: ibbar
    use mFuncoesDeForma,   only: shgq, shlq
    use mPropGeoFisica, only: calcBiotCoefficient
    use mPropGeoFisica, only: calcMatrixUndrainedProperties, calcMatrixProperties
    use mPropGeoFisica, only: updatePorosity
    !variables input
    implicit none
    integer, intent(in)             :: conecNodaisElem(nen,numel)
    real(8), intent(in)             :: x(nsd,numnp)
    real*8, intent(in) :: u(ndofD, numnp)
    real*8 :: stress(nrowb, nintd, numel), stressS(nintD, numel), stressTotal(nrowB, nintD,numel)
    real*8 :: p(1,numnp), pInit(1,numnp)
    integer :: isUndrained
    
    !variables
    real*8, external :: rowdot, coldot
    !
    logical quad
    real(8), dimension(nesd,nen)         :: xl
    real(8), dimension(ned2,nen)         :: disl
    real(8), dimension(nintd)            :: wd, detd, r
    real(8), dimension(nrowsh,nen,nintd) :: shld,shgd
    real(8), dimension(nrowsh,nen)       :: shgbr
    real(8), dimension(nrowb,nrowb)      :: cbbar
    real(8), dimension(nrowb,nesd)  :: bbarj
    real*8 :: strainL(nrowb)
    real*8 :: pL(1,nen), pInitL(1,nen)
    real*8 :: pressureIntPoint, pressureInitIntPoint, biotCoeff, identI(nrowb)
    !
    real*8  :: poisson, young
    integer :: nel, l, j, k
    
    real*8 :: c1
    real*8 :: area
    
    real*8 :: elTrStrain
    real*8 :: elP, elPInit
    
    real*8, external :: traceTensor
    !
    !---------------------------------------------------------------------
    identI = (/ 1., 1., 0., 1. /)
    
    !clean stress, strain and initial force
    stress = 0.d0
    
    !.... generation of local shape functions and weight values
    !
    call shlq(shld,wd,nintd,nen)
    !
    do nel=1,numel
        !
        if (isUndrained == 1) then
            call calcMatrixUndrainedProperties(poisson,young, nel)
        else if (isUndrained == 0) then
            call calcMatrixProperties(poisson, young, nel)
        end if
        !
        !.... ..setup stochastic elasticity tensor for bbar method
        !
        call setupc(cbbar, young, poisson,nrowb,iopt)
        !
        !...  ..localize coordinates and dirichlet b.c.
        !
        call local(conecNodaisElem(1,nel), x, xl, nen, nsd, nesd)
        call local(conecNodaisElem(1,nel), u, disl, nen, ndofd, ned2)
        call local(conecNodaisElem(1,nel), p, pL, nen, 1, 1)
        call local(conecNodaisElem(1,nel), pInit, pInitL, nen, 1, 1)
        
        !.... ... calculates biot coefficient
        biotCoeff = calcBiotCoefficient(nel)
        !
        quad = .true.
        if (conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
        !
        call shgq(xl,detd,shld,shgd,nintd,nel,quad,nen)
        !..
        !.... ...setup for axisymmetric option
        !..
        if (iopt.eq.2) then
            do l=1,nintd
                r(l)    = rowdot(shgd(nrowsh,1,l),xl,nrowsh,nesd,nen)
                detd(l) = detd(l)*r(l)
            end do
        endif
        !
        !.... ..calculate mean values of shape function global derivatives
        !.... ..for mean-dilatational b-bar formulation
        !
        call xmeansh(shgbr,wd,detd,r,shgd,nen,nintd,iopt,nesd,nrowsh)
        
        area= 0.d0
        elP = 0.d0
        elPInit = 0.d0
        elTrStrain = 0.d0
        !
        !.... ..loop over integrationn points
        !
        do l=1,nintd
            ! get pressure value at integration point
            pressureIntPoint = dot_product(shgD(nrowsh,1:nen,l), pL(1,1:nen))
            pressureInitIntPoint = dot_product(shgD(nrowsh,1:nen,l), pInitL(1,1:nen))
            
            !..... ..clear initial strain
            strainL = 0.d0
            !
            do j=1,nen
                !
                !.... ..upload b-bar matrix at node j
                !
                call setbb(bbarj,shgd(1:nrowsh,j,l),shgbr(1:nrowsh,j), r(l),nrowsh,nrowb,iopt,ibbar)
                !
                !.... ..compute strains within intrinsic b-bar formulation
                !
                ! B \cdot a
                do k=1,nrowb
                    strainL(k) = strainL(k) + dot_product(bbarj(k,1:2),disL(1:2,j))
                end do
                
            end do
            
            stress(1:nrowb,l,nel) = matmul(cbbar,strainL)
            
            do k=1,nrowb
                stressTotal(k,l,nel) = stress(k,l,nel) - biotCoeff*pressureIntPoint*IdentI(k)
            end do
            stressS(l,nel) = traceTensor(stressTotal(:,l,nel),nrowb)/3.0d0
            
            !updatePorosity
            c1 = detd(l)*wd(l)
            area = area + c1
            
            elTrStrain = elTrStrain + traceTensor(strainL,nrowb)*c1
            elP = elP + pressureIntPoint*c1
            elPInit = elPInit + pressureInitIntPoint*c1
        end do
        
        !update porosity
        elTrStrain = elTrStrain / area
        elP = elP / area
        elPInit = elPInit / area
        
        call updatePorosity(nel, elTrStrain, 0.d0, elP, elPInit)
    end do
    
    end subroutine pos4inc
    !**************************************************************************************
    !**************************************************************************************
    
    SUBROUTINE UPDTSTRS(STRSS,STRSINIT,NROWB,NUMEL)
    !
    !.... PROGRAM TO UPDATE EFFECTIVE STRESS WITHIN INITIAL EFFECTIVE STRESS
    !
    IMPLICIT NONE
    !
    INTEGER :: NEL,NUMEL,NROWB, K
    !
    REAL(8), DIMENSION(NROWB,NUMEL) :: STRSS, STRSINIT
    !
    DO 200 NEL=1,NUMEL
        DO 100 K=1,NROWB
            STRSS(K,NEL) = STRSS(K,NEL)+STRSINIT(K,NEL)
100     CONTINUE
200 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !
    !***** %%% ***** %%% ***** %%% ***** %%% ***** %%% ***** %%% ****
    subroutine yield(fyield,sigmab,taub,tensao,poisson, misesyield)
    !
    !
    !..... program to compute gradient (fnabla) and hessian (fhessi)
    !..... of mohr-columb in function for normal and shear stress
    !..... as indendent variables.
    !
    real(8) :: fyield, sigmab, taub, poisson, vmpoisson, misesyield
    real(8), dimension(nrowb) :: tensao
    !
    sigmab = 0.5d0*(tensao(1)+tensao(2))
    taub   = dsqrt((0.5d0*(tensao(1)-tensao(2)))**2+tensao(3)**2)
    !
    !.... von mises criteria
    !
    vmpoisson = 2.0d0*(1.0d0-2.0d0*poisson)**2
    fyield    = vmpoisson*sigmab**2+6.0d0*taub**2-6.0d0*misesyield**2
    !
    return

    !.... mohr-coulomb criteria
    !
    !.... iopt=4---> linear classical mohr-coulomb yield
    !....           for (sigma_n,tau) stress formulation
    !
    !      fyield = taub + sigmab*tanfri-cohesion
    !
    !return
    !
    end subroutine yield
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE GRADS(GRAD,QIXI,TENSAO,HGAMMA,CBBARM1,POISSON,nrowb)
    !
    !.... PROGRAM TO CALCULATE YIELD'S GRADIENT: GRAD,
    !.... AND THE HINVERSE OF QIXI
    !
    IMPLICIT NONE
    !
    !.... MATRIX OF IN/OUT FUNCTION
    !
    INTEGER :: I, J, NROWB
    REAL(8) :: HGAMMA, POISSON
    REAL(8), DIMENSION(NROWB)       :: GRAD, TENSAO
    REAL(8), DIMENSION(NROWB,NROWB) :: QIXI
    !
    !.... MATRIX OF LOCAL COMPUTATIONS
    !
    REAL(8), DIMENSION(NROWB,NROWB) :: CBBARM1, GRAD2, HINVQIXI
    !
    GRAD  = 0.0D0
    GRAD2 = 0.0D0
    !
    !.... COMPUTE GRADIENT AND HESSIAN OF BASIC FUNCTIONS (SIGMA_B,TAU_B)
    !
    CALL SGTAGRAD(GRAD,GRAD2,TENSAO,POISSON,NROWB)
    !
    !.... QIXI INVERSE DEFINITION (SIMO PAG.175. EQUATION 4.3.15)
    !
    DO 200 I=1,NROWB
        DO 100 J=1,NROWB
            HINVQIXI(I,J) = CBBARM1(I,J) + HGAMMA*GRAD2(I,J)
100     CONTINUE
200 CONTINUE
    !
    !.... COMPUTE INVERSE OF SIMETRIC MATRIX HINVQIXI
    !
    CALL COMPQIXI(QIXI,HINVQIXI,NROWB)
    !
    !
    RETURN
    !
2222 FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X  &
        &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                  &
        &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                  &
        &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                  &
        &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
    !
    END SUBROUTINE GRADS
    !**************************************************************************************
    !**************************************************************************************
    subroutine dpYieldEfStress(fYield,efStress,p, biotCoef, k, alpha)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: fYield, efStress(nrowb)
    real*8 :: p, biotCoef
    real*8 :: k, alpha
    
    !variables
    real*8 :: dummyFYield1(nrowb), dummyFYield2(nrowb,nrowb)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    call dpGrads(fYield, dummyFYield1, dummyFYield2, efStress, p, biotCoef, k, alpha, .true.)
    
    end subroutine dpYieldEfStress
    !**************************************************************************************
    !**************************************************************************************
    function dpYield(I1, J2, k, alpha)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: dpYield, I1, J2, k, alpha
    
    !variables
    
    !--------------------------------------------------------------------------------------
    dpYield = dsqrt(J2) + alpha * I1/3.d0 - k
    
    end function dpYield
    !**************************************************************************************
    !**************************************************************************************
    function dpYield1(devStress, J2, alpha)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: devStress(nrowb), J2, alpha
    
    !variables
    real*8 :: dpYield1(nrowb)
    
    !-----------------------------------------
    dpYield1(1) = devStress(1)/2.d0/dsqrt(J2) + alpha/3.d0
    dpYield1(2) = devStress(2)/2.d0/dsqrt(J2) + alpha/3.d0
    dpYield1(3) = devStress(3)/dsqrt(J2)
    dpYield1(4) = devStress(4)/2.d0/dsqrt(J2) + alpha/3.d0
    
    end function dpYield1
    !*****************************************
    !*****************************************
    function dpYield2(plasticStress,I1,J2)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: plasticStress(nrowb), I1, J2
    
    !variables
    real*8 dpYield2(nrowb,nrowb)
    
    !--------------------------------------------------------------------------------------
    dpYield2(1,1) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(1)*plasticStress(4)-I1*plasticStress(4)+3*plasticStress(1)*plasticStress(2)-I1*plasticStress(2)-6*plasticStress(1)**2+2*I1*plasticStress(1)+12*J2))/3.6d+1
    dpYield2(1,2) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(1)*plasticStress(4)-I1*plasticStress(4)-6*plasticStress(1)*plasticStress(2)+2*I1*plasticStress(2)+3*plasticStress(1)**2-I1*plasticStress(1)-6*J2))/3.6d+1
    dpYield2(1,3) = -(J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(1)-I1)*plasticStress(3))/6.0d+0
    dpYield2(1,4) = -(J2**((-3.0d+0)/2.0d+0)*(6*plasticStress(1)*plasticStress(4)-2*I1*plasticStress(4)-3*plasticStress(1)*plasticStress(2)+I1*plasticStress(2)-3*plasticStress(1)**2+I1*plasticStress(1)+6*J2))/3.6d+1

    dpYield2(2,1) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(2)*plasticStress(4)-I1*plasticStress(4)+3*plasticStress(2)**2-6*plasticStress(1)*plasticStress(2)-I1*plasticStress(2)+2*I1*plasticStress(1)-6*J2))/3.6d+1
    dpYield2(2,2) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(2)*plasticStress(4)-I1*plasticStress(4)-6*plasticStress(2)**2+3*plasticStress(1)*plasticStress(2)+2*I1*plasticStress(2)-I1*plasticStress(1)+12*J2))/3.6d+1
    dpYield2(2,3) = -(J2**((-3.0d+0)/2.0d+0)*plasticStress(3)*(3*plasticStress(2)-I1))/6.0d+0
    dpYield2(2,4) = -(J2**((-3.0d+0)/2.0d+0)*(6*plasticStress(2)*plasticStress(4)-2*I1*plasticStress(4)-3*plasticStress(2)**2-3*plasticStress(1)*plasticStress(2)+I1*plasticStress(2)+I1*plasticStress(1)+6*J2))/3.6d+1

    dpYield2(3,1) = (J2**((-3.0d+0)/2.0d+0)*plasticStress(3)*(plasticStress(4)+plasticStress(2)-2*plasticStress(1)))/6.0d+0
    dpYield2(3,2) = (J2**((-3.0d+0)/2.0d+0)*plasticStress(3)*(plasticStress(4)-2*plasticStress(2)+plasticStress(1)))/6.0d+0
    dpYield2(3,3) = -J2**((-3.0d+0)/2.0d+0)*(plasticStress(3)**2-J2)
    dpYield2(3,4) = -(J2**((-3.0d+0)/2.0d+0)*plasticStress(3)*(2*plasticStress(4)-plasticStress(2)-plasticStress(1)))/6.0d+0

    dpYield2(4,1) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(2)*plasticStress(4)-6*plasticStress(1)*plasticStress(4)+I1*plasticStress(4)-3*plasticStress(2)**2+I1*plasticStress(2)-6*plasticStress(3)**2-3*plasticStress(1)**2+4*I1*plasticStress(1)-I1**2))/3.6d+1
    dpYield2(4,2) = -(J2**((-3.0d+0)/2.0d+0)*(6*plasticStress(2)*plasticStress(4)-3*plasticStress(1)*plasticStress(4)-I1*plasticStress(4)+3*plasticStress(2)**2-4*I1*plasticStress(2)+6*plasticStress(3)**2+3*plasticStress(1)**2-I1*plasticStress(1)+I1**2))/3.6d+1
    dpYield2(4,3) = -(J2**((-3.0d+0)/2.0d+0)*plasticStress(3)*(3*plasticStress(4)-I1))/6.0d+0
    dpYield2(4,4) = (J2**((-3.0d+0)/2.0d+0)*(3*plasticStress(2)*plasticStress(4)+3*plasticStress(1)*plasticStress(4)-2*I1*plasticStress(4)+6*plasticStress(2)**2-5*I1*plasticStress(2)+12*plasticStress(3)**2+6*plasticStress(1)**2-5*I1*plasticStress(1)+2*I1**2))/3.6d+1

    end function dpYield2
    !**************************************************************************************
    !**************************************************************************************
    subroutine dpGrads(fYield, fYield1, fYield2, efStress, p, biotCoef, k, alpha, onlyF)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: fYield, fYield1(nrowb), fYield2(nrowb,nrowb)
    real*8 :: efStress(nrowb), p, biotCoef, k, alpha
    logical :: onlyF
    
    !variables
    real*8 :: plasticStress(nrowb), devStress(nrowb), traceD3
    real*8 :: I1, J2
    
    !--------------------------------------------------------------------------------------
    ! calculates the effective plastic stress for the drucker prager model (dormieux(2006) pag.225)
    plasticStress = calcDPPlasticStress(efStress, p, biotCoef)
    
    !calculates the deviator of the stress tensor
    call devTensor(devStress,traceD3,plasticStress,nrowb)
    
    !calculates auxiliary values for the f function
    I1 = calcInvariantI1(plasticStress)
    J2 = calcDevInvariantJ2(devStress)

    if (dabs(J2) < 1d-38) J2 = 1d-38
    
    fYield = dpYield(I1, J2, k, alpha)
    if (onlyF) return
    fYield1 = dpYield1(devStress, J2, alpha)
    fYield2 = dpYield2(plasticStress, I1, J2)
    
    end subroutine dpGrads
    !**************************************************************************************
    !**************************************************************************************
    subroutine calcQixiM1(qixi,cbbarm1,hGamma,fYield2)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: qixi(nrowb,nrowb), cbbarm1(nrowb,nrowb), hGamma, fYield2(nrowb,nrowb)
    
    !variables
    real*8 :: hInvQIxi(nrowb,nrowb)
    integer :: i,j
    
    !--------------------------------------------------------------------------------------
    do j=1,nrowb
        do i=1,nrowb
            hInvQIxi(I,J) = cbbarm1(i,j) + hgamma*fYield2(i,j)
        end do
    end do
    call compQixi(qixi,hInvQIxi,nrowb)
    
    end subroutine calcQixiM1
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE SGTAGRAD(GRAD,GRAD2,TENSAO,POISSON,NROWB)
    !
    !..... PROGRAM TO COMPUTE GRADIENT (grad) AND HESSIAN (grad2)
    !..... FOR MISES AND MOHR-COULOMB CAP MODEL.
    !
    IMPLICIT NONE
    !
    INTEGER :: I, J,  NROWB
    REAL(8) :: POISSON
    REAL(8), DIMENSION(NROWB)  :: GRAD, TENSAO
    REAL(8), DIMENSION(NROWB,NROWB) :: GRAD2
    !
    REAL(8), DIMENSION(2,3) :: BNABLA
    REAL(8), DIMENSION(3,3) :: BHESSI
    REAL(8), DIMENSION(2)   :: FNABLA
    REAL(8), DIMENSION(2,2) :: FHESSI
    !
    !.... COMPUTE BASIC GRADIENTS (SEE INFORMATION FILE: contas2.pdf)
    !
    CALL BASGRAD(BNABLA,BHESSI,TENSAO,NROWB)
    !
    !.... COMPUTE GRADIENTS OF YIELD FUNCTION IN (SIGMA_B,TAU_B) COORDINATES
    !
    CALL FMCGRAD(FNABLA,FHESSI,TENSAO,POISSON,NROWB)
    !
    !..... COMPUTE GRADIENTS YIELD CRITERION AS CHAIN RULE PRODUCT
    !.....     SEE contas_new_code.pdf FILE
    !
    DO 12 I = 1,3
        DO 11 J = 1,2
            GRAD(I) = GRAD(I) + FNABLA(J)*BNABLA(J,I)
11      CONTINUE
12  CONTINUE
    !
    GRAD(4) = 0.0D0
    !
    GRAD2   = 0.0D0
    !
    !..... COMPUTE HESSIAN MATRIX OF YIELD FUNCTION
    !
    GRAD2(1,1)=BNABLA(1,1)*FHESSI(1,1)*BNABLA(1,1)+ &
        &           BNABLA(2,1)*FHESSI(2,2)*BNABLA(2,1)+ &
        &           FNABLA(2)*BHESSI(1,1)
    !
    GRAD2(1,2)=BNABLA(1,1)*FHESSI(1,1)*BNABLA(1,2)+ &
        &           BNABLA(2,1)*FHESSI(2,2)*BNABLA(2,2)+ &
        &           FNABLA(2)*BHESSI(1,2)
    !
    GRAD2(1,3)=BNABLA(2,1)*FHESSI(2,2)*BNABLA(2,3)+ &
        &           FNABLA(2)*BHESSI(1,3)
    !
    GRAD2(2,1)=GRAD2(1,2)
    !
    GRAD2(2,2)=BNABLA(1,2)*FHESSI(1,1)*BNABLA(1,2)+ &
        &           BNABLA(2,2)*FHESSI(2,2)*BNABLA(2,2)+ &
        &           FNABLA(2)*BHESSI(2,2)
    !
    GRAD2(2,3)=BNABLA(2,2)*FHESSI(2,2)*BNABLA(2,3)+ &
        &           FNABLA(2)*BHESSI(2,3)
    !
    GRAD2(3,1)=GRAD2(1,3)
    !
    GRAD2(3,2)=GRAD2(2,3)
    !
    GRAD2(3,3)=BNABLA(2,3)*FHESSI(2,2)*BNABLA(2,3)+ &
        &           FNABLA(2)*BHESSI(3,3)
    !
    RETURN
    !
    END SUBROUTINE SGTAGRAD
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE BASGRAD(BNABLA,BHESSI,TENSAO,NROWB)
    !
    !..... PROGRAM TO COMPUTE GRADIENT (BNABLA) AND HESSIAN (BHESSI)
    !..... OF BASIC FUNCTIONS CALLED \SIGMA_B ANR \TAU_B
    !.....
    !..... SUCH DERIVATIVES WE CALLED BASIC GRADIENTS
    !..... BECAUSE THEY ARE INDEPENDT OF THE YIELD CRITERIA.
    !
    IMPLICIT NONE
    !
    INTEGER NROWB
    REAL(8) :: ROOT, ROOTM1, P25RM13, SXMSY
    REAL(8), DIMENSION(2,3)   :: BNABLA
    REAL(8), DIMENSION(3,3)   :: BHESSI
    REAL(8), DIMENSION(NROWB) :: TENSAO
    !
    SXMSY   = TENSAO(1)-TENSAO(2)
    ROOT    = DSQRT(0.25D0*(SXMSY)**2 + TENSAO(3)**2)
    ROOTM1  = 1.0D0/ROOT
    P25RM13 = 0.25D0*(ROOTM1)**3
    !
    !.... BNABLA(1,*) GRADIENT OF TRACE STRESS
    !           (1,1) Drho /Drho X
    BNABLA(1,1) = 0.5D0
    !           (1,2) Drho /Drho Y
    BNABLA(1,2) = 0.5D0
    !           (1,3) Drho /Drho XY
    BNABLA(1,3) = 0.0D0
    !
    !.... BNABLA(2,*) GRADIENT OF SQUARE ROOT
    !           (2,1) Drho /Drho X
    BNABLA(2,1) =  0.25D0*SXMSY*ROOTM1
    !           (2,2) Drho /Drho Y
    BNABLA(2,2) = -BNABLA(2,1)
    !           (2,3) Drho /Drho XY
    BNABLA(2,3) =  TENSAO(3)*ROOTM1
    !
    !.... HESSIAN OF TRACE == ZERO
    !.... HESSIAN OF SQUARE ROOT :
    !           (1,1) Drho2/(Drho X Drho X)
    BHESSI(1,1) =  P25RM13*(TENSAO(3))**2
    !           (1,2) Drho2/(Drho X Drho Y)
    BHESSI(1,2) = -BHESSI(1,1)
    !           (1,3) Drho2/(Drho X Drho XY)
    BHESSI(1,3) = -P25RM13*(SXMSY)*TENSAO(3)
    !           (2,1) Drho2/(Drho Y Drho X)
    BHESSI(2,1) =  BHESSI(1,2)
    !           (2,2) Drho2/(Drho Y Drho Y)
    BHESSI(2,2) =  BHESSI(1,1)
    !           (2,3) Drho2/(Drho Y Drho XY)
    BHESSI(2,3) = -BHESSI(1,3)
    !           (3,1) Drho2/(Drho XY Drho X)
    BHESSI(3,1) =  BHESSI(1,3)
    !           (3,2) Drho2/(Drho XY Drho Y)
    BHESSI(3,2) =  BHESSI(2,3)
    !           (2,3) Drho2/(Drho XY Drho XY)
    BHESSI(3,3) =  P25RM13*SXMSY**2
    !
    RETURN
    !
    END SUBROUTINE BASGRAD
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE FMCGRAD(FNABLA,FHESSI,TENSAO,POISSON,NROWB)
    !
    !..... PROGRAM TO COMPUTE GRADIENT (FNABLA) AND HESSIAN (FHESSI)
    !..... OF MOHR-COLUMB IN FUNCTION OF NORMAL AND SHEAR STRESS
    !..... AS INDENDENTS VARIABLES.
    !
    IMPLICIT NONE
    !
    INTEGER :: NROWB
    REAL(8) :: SIGMAB, TAUB, POISSON, VMPOISSON
    REAL(8), DIMENSION(2)     :: FNABLA
    REAL(8), DIMENSION(2,2)   :: FHESSI
    REAL(8), DIMENSION(NROWB) :: TENSAO
    !
    !     GO TO (200,300,400,500), IOPT-1
    !
    SIGMAB = 0.5D0*(TENSAO(1)+TENSAO(2))
    TAUB   = DSQRT((0.5D0*(TENSAO(1)-TENSAO(2)))**2+TENSAO(3)**2)
    !
    !
    ! 300   CONTINUE
    !
    !.... IOPT=3---> VON MISES CRITERION (SIGMA_B,TAU_B) SPACE
    !
    VMPOISSON = 4.0D0*(1.0D0-2.0D0*POISSON)**2
    !
    !..... GRADIENT
    !
    FNABLA(1) = VMPOISSON*SIGMAB
    FNABLA(2) = 12.0D0*TAUB
    !
    !..... HESSIAN
    !
    FHESSI(1,1) = VMPOISSON
    FHESSI(1,2) = 0.0D0
    FHESSI(2,1) = 0.0D0
    FHESSI(2,2) = 12.0D0
    !
    RETURN
    !
    END SUBROUTINE FMCGRAD
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE COMPQIXI(QIXI,HINVQIXI,NROWB)
    !
    !.... PROGRAM TO COMPUTE INVERSE OF HINVQIXI MATRIX
    !.... OUTPUT IS MATRIX QIXI
    !
    IMPLICIT NONE
    !
    INTEGER :: I , J, NORDER, NROWB
    REAL(8) :: XDET
    !
    !.... INPUT/OUTPUT MATRIZES
    !
    REAL(8), DIMENSION(NROWB,NROWB) ::  QIXI, HINVQIXI
    !
    !.... LOCAL ARRAYS
    !
    REAL(8), DIMENSION(NROWB,NROWB)     :: ATEMP, COFATOR
    REAL(8), DIMENSION(NROWB-1,NROWB-1) :: AMINOR
    !
    ! external
    real*8, external :: detm3x3
    !
    !.... MOVE HINVQIXI TO ATEMP
    !
    NORDER = NROWB
    !
    DO 20 I=1,NORDER
        DO 10 J=1,NORDER
            ATEMP(I,J)=HINVQIXI(I,J)
10      CONTINUE
20  CONTINUE
    !
    DO 60 I=1,NORDER
        DO 50 J=1,NORDER
            !
            !.... .... COMPUTE MATRIX OF MINORS
            !
            CALL COMPMINOR(AMINOR,ATEMP,I,J,NORDER)
            !
            !.... .... COMPUTE MATRIX CO-FATORS
            !
            COFATOR(I,J)=((-1)**(I+J))*DETM3X3(AMINOR)
50      CONTINUE
60  CONTINUE
    !
    !..... COMPUTE DETERMINANT
    !
    XDET = 0.0D0
    !
    DO 70 I=1,NORDER
        XDET=XDET+ATEMP(1,I)*COFATOR(1,I)
70  CONTINUE
    !
    DO 90 I=1,NORDER
        DO 80 J=1,NORDER
            QIXI(I,J)=COFATOR(J,I)/XDET
80      CONTINUE
90  CONTINUE
    !
    RETURN
    !
2000 FORMAT(A15,4(1PE10.2,2X))
    !
    END SUBROUTINE COMPQIXI
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE TRIALS(EPTRIAL,HGAMMA,QIXIGRAD,QIXI,RESIDUO, GRAD,FYIELD,CBBARM1,NROWB)
    !
    use mMalha,            only: multab
    !
    !.... PROGRAM TO UPDATE PLASTIC AND GAMMA TRIAL'S
    !.... ALSO COMPUTE QIXIGRAD VETOR (NAME N SIMO VETOR)
    !.... TO USE IN CONSISTENT TANGENT MODULI.
    !
    !.... INPUT DATA: QIXI, RESIDUO, GRAD, FMOHR
    !.... OUTPUT DATA: EPTRIAL, HGAMMA, QIXIGRAD
    !
    IMPLICIT NONE
    !
    INTEGER :: I,NROWB
    REAL(8) :: HGAMMA, FYIELD, DHGAMMA, DENOMIN

    real*8, external :: coldot
    !
    !.... GLOBAL VETORS AND MATRIZES
    !
    REAL(8), DIMENSION(NROWB) :: EPTRIAL, GRAD, QIXIGRAD, RESIDUO
    REAL(8), DIMENSION(NROWB,NROWB) :: QIXI, CBBARM1
    !
    !.... LOCAL VETORS AND MATRIZES
    !
    REAL(8), DIMENSION(NROWB) :: QIXIRES, RESLOC, DPLAS, GRADMRES
    !
    !.... PART: 2D FROM BOX 4.1 "COMPUTE INCREMENTS"
    !
    QIXIRES  = 0.0D0
    DPLAS    = 0.0D0
    RESLOC   = 0.0D0
    QIXIGRAD = 0.0D0
    GRADMRES = 0.0D0
    !
    !.... ... COMPUTATION OF QIXIRES=QIXI*RESIDUO
    !
    qixires = matmul(qixi,residuo)
    !
    !.... ... COMPUTATION OF QIXIGRD=QIXI*GRAD
    !
    qixigrad = matmul(qixi,grad)

    DENOMIN = 1.0D0/COLDOT(GRAD,QIXIGRAD,NROWB)
    !
    !.... .. HGAMMA INCREMENT:
    !
    DHGAMMA = (FYIELD+COLDOT(GRAD,QIXIRES,NROWB))*DENOMIN
    !
    !.... .. PLASTIC INCREMENT
    !
    DO 20 I=1,NROWB
        GRADMRES(I) = DHGAMMA*GRAD(I)-RESIDUO(I)
20  CONTINUE

    
    resloc = matmul(qixi, gradmres)
    
    !.... ...SOLVE C*DPLAS = RESLOC, FOR DPLAS --> DPLAS = Cm1*RESLOC
    dPlas = matmul(cbbarm1, resloc)
    !
    !.... PARTE 2E BOX 4.1 UPDATE PLASTIC AND HGAMMA TRIALS
    !
    HGAMMA = HGAMMA + DHGAMMA
    !      WRITE(*,*)'HGAMMA= ',HGAMMA,' DHGAMMA = ',DHGAMMA
    !
    DO 100 I=1,NROWB
        EPTRIAL(I) = EPTRIAL(I) + DPLAS(I)
100 CONTINUE

    RETURN
    !
    END SUBROUTINE TRIALS
    !**************************************************************************************
    !**************************************************************************************
    SUBROUTINE STRSS4PLAST(x, conecNodaisElem)
    !
    !.... PROGRAM TO UPDATE STRESS FOR NON-LINEAR CREEP MODEL
    !
    use mMalha,            only: nsd, numnp, numel, nen, LOCAL
    use mMalha,            only: multab
    use mGlobaisEscalares, only: nrowsh
    use mFuncoesDeForma,   only: shgq, shlq
    use mGlobaisEscalares, only: ibbar
    use mPropGeoFisica,    only: YOUNG
    use mPropGeoFisica,    only: GEOFORM, GEOINDIC, FNCMECLAW
    !
    IMPLICIT NONE
    !
    REAL*8,  intent(in)    :: x(nsd,numnp)
    INTEGER, intent(in)    :: conecNodaisElem(nen,numel)
    CHARACTER(5)           :: MECLAW
    
    real*8, external :: coldot, rowdot
    !
    LOGICAL QUAD
    !
    INTEGER :: J,K,L,NEL
    !
    !.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
    !
    REAL(8), DIMENSION(NESD,NEN)      :: XL
    REAL(8), DIMENSION(NED2,NEN)      :: DLTRL
    !
    !.... LOCAL VECTORS AND MATRIZES
    !
    REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ
    REAL(8), DIMENSION(NROWB,NROWB)   :: CBBAR
    REAL(8), DIMENSION(NROWB)         :: TENSAO, STRAIN, EELAS
    !
    REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
    REAL(8) :: SHGBR(NROWSH,NEN)
    REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)

    REAL(8) :: POISSON
    !
    !.... NEXT FOR PLASTICITY
    !
    !
    !.... STRAIN      :  TOTAL DEFORMATION
    !.... EPLAS       :  deformacao plastica
    !.... Einelas     :  PLASTICA
    !
    !.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES
    !
    CALL SHLQ(SHLD,WD,NINTD,NEN)
    !
    TENSAO = 0.0D0
    !
    DO 500 NEL=1,NUMEL
        !
        MECLAW  = FNCMECLAW(GEOFORM(NEL))
        !
        POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
        !
        !.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD
        !
        CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
        !
        !..... ..LOCALIZE COORDINATES AND DIRICHLET B.C.
        !
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),DTRL,DLTRL,NEN,NDOFD,NED2)
        !
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
        !
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
        !..
        !..... ..SETUP FOR AXISYMMETRIC OPTION
        !..
        IF (IOPT.EQ.2) THEN
            DO 100 L=1,NINTD
                R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
                DETD(L) = DETD(L)*R(L)
100         CONTINUE
        ENDIF
        !
        !.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
        !....     FOR MEAN-DILATATIONAL B-BAR FORMULATION
        !
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
        !..
        !.... ...LOOP OVER INTEGRATION POINTS
        !..
        DO 400 L=1,NINTD
            !..
            !..... ..CLEAR INITIAL STRAIN
            !..
            STRAIN = 0.0D0
            !
            DO 200 J=1,NEN
                !...
                !.... ..... UPLOAD B-BAR MATRIX AT NODE J
                !..
                CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
                    &                    R(L),NROWSH,NROWB,IOPT,IBBAR)
                !..
                !.... ..... COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION
                !..
                DO 200 K=1,NROWB
                    STRAIN(K)=STRAIN(K)+ &
                        &                    COLDOT(BBARJ(K,1:2),DLTRL(1:2,J),2)
200         CONTINUE
            !..
            !.... ..... COMPUTE ELASTIC DEFORMATION
            !..
            DO 210 K=1,NROWB
                EELAS(K) = STRAIN(K) - EINELAS(NEL,L,K)
210         CONTINUE
            !
            !.... ..... COMPUTE STRESS
            !..
            CALL MULTAB(CBBAR,EELAS,TENSAO,4,4,4,4,4,1,1)
            !
            !.... ..... LOCATE ON STRESS FOR GLOBAL FORCE BALANCE
            !
            DO 250 K=1,NROWB
                STRSSP(NEL,L,K) = TENSAO(K)
250         CONTINUE
            !
400     CONTINUE
        !
500 CONTINUE
    !
    !      stop
    RETURN
    !
2222 FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X  &
        &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                  &
        &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                  &
        &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                  &
        &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
4000 FORMAT(2X,40(1PE15.8,2X))
4001 FORMAT('nel=',I5,x,'gauss=',I1,x,'I',i1,2X,40(1PE15.8,2X))
5000 FORMAT(I4,2X,I1,2X,40(1PE15.8,2X))
    !
    END SUBROUTINE STRSS4PLAST
    !**************************************************************************************
    !**************************************************************************************
    subroutine incrementMechanicElasticSolution(conecNodaisElem, nen, nel, nnp, nsd, x, curT, u, stress, stressS, stressTotal, p, pInit, isUndrained)
    !function imports
    use mSolverGaussSkyline, only: solverGaussSkyline
    
    !variables import
    
    implicit none

    !variables input
    integer :: nen, nel, nnp, nsd, conecNodaisElem(nen, nel)
    real*8 :: x(nsd, nnp), curT, u(ndofD,nnp)
    real*8 :: stress(nrowB,nintD,nel), stressS(nintD, nel), stressTotal(nrowB, nintD,nel)
    real*8 :: p(1,nnp), pInit(1,nnp)
    integer :: isUndrained

    ! Variables
    real*8, allocatable :: fExtT(:,:), fIntJ(:,:) !fExtT - external force at time T, fIntJ - internal Force at newton iteration J
    real*8, allocatable :: fExtNeumann(:,:), fExtDirichlet(:,:) !fExtT with zeros on dirichlet/neumann nodes
    real*8, allocatable :: dDis(:,:)
    real*8 :: error
    
    real*8 :: g1(nlvectD)
    real*8 :: g2(nlvectD)
    
    real*8, external :: matrixNorm
    integer :: j
    
    real*8 :: tolNewton
    logical :: converged

    !--------------------------------------------------------------------------------------

    ! allocate the required matrices
    if(.not.allocated(fExtT)) allocate(fExtT(ndofD, nnp))
    if(.not.allocated(fExtNeumann)) allocate(fExtNeumann(ndofD, nnp))
    if(.not.allocated(fExtDirichlet)) allocate(fExtDirichlet(ndofD, nnp))
    if(.not.allocated(fIntJ)) allocate(fIntJ(ndofD, nnp))
    fExtT = 0.d0
    fIntJ = 0.d0
    g2 = -1.d0

    if(.not.allocated(dDis)) allocate(dDis(ndofD, nnp))
    dDis = 0.d0
    
    tolNewton = 1.0D-6
    write (*,*) "Incremento elástico"

    !compute the new external force vector, and split into two, a dirichlet and a neumann one
    call splitBoundaryCondition(idDesloc,gmF,fExtDirichlet,fExtNeumann,ndofD,nnp,nlvectD)
    
    !computes the strain and stress
    call pos4inc(x, conecNodaisElem, u, stress, stressS, stressTotal, p, pInit, 0)

    !computes the internal force
    call calcInternalForce(x, conecNodaisElem, nen, nel, nnp, nsd, stress, fIntJ)

    converged = .false.
    do j = 1, 15
        alhsD = 0.d0
        brhsD = 0.d0
        dDis = 0.d0
        
        ! load force vector in the right side
        if (nlvectD.gt.0) then
            if (nltftnD >= 1) then
                call lfac(gmG, curT, g1, nltftnD, nptslfD)
            else
                g1 = 1.d0
            endif
            call load(idDesloc, fExtNeumann, brhsd, g1, ndofD, nnp, nlvectD)
            call load(idDesloc, fIntJ, brhsd, g2, ndofD, nnp, nlvectD) !g2 is always -1
            if (j == 1) then !only do this on the first iteration (Non-linear Finite Element Analysis of Solids and Structures pg.51)
                call ftodDif(idDesloc, dDis, fExtDirichlet, u, g1, ndofD, nnp, nlvectD)
            end if
        end if

        ! compute the tangential stiffness matrix
        call bbarmtrx_elast(x, conecNodaisElem, dDis, alhsD, brhsD, idiagD, lmD, p, isUndrained)
        
        error = dsqrt(dot_product(brhsD,brhsD))/neqD
        write(*,*) "erro do passo mecanico: ", error
        if (error < tolNewton .and. j > 1) then
            converged = .true.
            exit
        end if

        ! solve using LDU decomposition, and store the result in the right array
        call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'full')
        call btod(idDesloc, dDis, brhsD, ndofD, nnp)

        ! add correction to the incremental displacement vector
        u = u + dDis

        !computes the strain and stress
        call pos4inc(x, conecNodaisElem, u, stress, stressS, stressTotal, p, pInit, 0)

        !computes the internal force
        call calcInternalForce(x, conecNodaisElem, nen, nel, nnp, nsd, stress, fIntJ)
    end do
    if (converged.eqv..false.) write(*,*) 'Newton method not converged.'

    
    deallocate(fExtT)
    deallocate(fExtNeumann)
    deallocate(fIntJ)
    deallocate(dDis)

    end subroutine incrementMechanicElasticSolution
    !**************************************************************************************
    !**************************************************************************************
    subroutine incrementMechanicPlasticSolution(conecNodaisElem, nen, nel, nnp, nsd, x, curT, u, strainP, prevStrainP, stress, stressS, trStrainP, stressTotal, p, pInit, biotP, elementIsPlast, isUndrained)
    !function imports
    use mSolverGaussSkyline, only: solverGaussSkyline
    use mLeituraEscritaSimHidroGeoMec, only: escreverArqParaviewIntermed_CampoVetorial, escreverArqParaviewIntermed_CampoEscalar
    
    implicit none
    !variables input
    integer :: nen, nel, nnp, nsd, conecNodaisElem(nen, nel)
    real*8 :: x(nsd, nnp), curT, u(ndofD,nnp), strainP(nrowb, nintD, nel), prevStrainP(nrowb,nintD,nel)
    real*8 :: stress(nrowB,nintD,nel), stressS(nintD, nel), trStrainP(nintD,nel), stressTotal(nrowB, nintD,nel)
    real*8 :: p(1,nnp), pInit(1,nnp)
    real*8 :: biotP(:,:,:)
    real*8 :: elementIsPlast(nel)
    integer :: isUndrained

    ! Variables
    real*8, allocatable :: fExtT(:,:), fIntJ(:,:) !fExtT - external force at time T, fIntJ - internal Force at newton iteration J
    real*8, allocatable :: fExtNeumann(:,:), fExtDirichlet(:,:) !fExtT with zeros on dirichlet/neumann nodes
    real*8, allocatable :: dDis(:,:)
    real*8 :: error
    
    real*8 :: g1(nlvectD)
    real*8 :: g2(nlvectD)
    
    real*8, external :: matrixNorm

    integer :: j
    
    real*8 :: tolNewton
    logical :: converged

    !--------------------------------------------------------------------------------------
    

    ! allocate the required matrices
    if(.not.allocated(fExtT)) allocate(fExtT(ndofD, nnp))
    if(.not.allocated(fExtNeumann)) allocate(fExtNeumann(ndofD, nnp))
    if(.not.allocated(fExtDirichlet)) allocate(fExtDirichlet(ndofD, nnp))
    if(.not.allocated(fIntJ)) allocate(fIntJ(ndofD, nnp))
    fExtT = 0.d0
    fIntJ = 0.d0
    g2=-1.0d0

    if(.not.allocated(dDis)) allocate(dDis(ndofD, nnp))
    
    tolNewton = 1.0d-6

    write(*,*) "Incremento plastico"
    call geoSetup(nel,nrowb,nintd,iopt, isUndrained, biotP)

    !compute the new external force vector, and split into two, a dirichlet and a neumann one
    call splitBoundaryCondition(idDesloc,gmF,fExtDirichlet,fExtNeumann,ndofD,nnp,nlvectD)
    
    !Init stress tensor and other variables
    call pos4plast(x, conecNodaisElem, u, strainP, prevStrainP, stress, stressS, trStrainP, stressTotal, hmTTG, biotP, elementIsPlast, p, pInit, isUndrained)

    !computes the internal force
    call calcInternalForce(x, conecNodaisElem, nen, nel, nnp, nsd, stress, fIntJ)

    converged = .false.
    do j = 1, 50
        alhsD = 0.d0
        brhsD = 0.d0
        dDis = 0.d0

        ! load force vector in the right side
        if (nlvectD.gt.0) then
            if (nltftnD >= 1) then
                call lfac(gmG, curT, g1, nltftnD, nptslfD)
            else
                g1 = 1.d0
            endif
            call load(idDesloc, fExtNeumann, brhsd, g1, ndofD, nnp, nlvectD)
            call load(idDesloc, fIntJ, brhsd, g2, ndofD, nnp, nlvectD)
            if (j==1) then !only do this on the first iteration (Non-linear Finite Element Analysis of Solids and Structures pg.51)
                call ftodDif(idDesloc, dDis, fExtDirichlet, u, g1, ndofD, nnp, nlvectD)
            end if
        end if
        
        ! compute the tangential stiffness matrix, and fill brhsd with the dirichlet node info
        call bbarmtrx_plast(x, conecNodaisElem, fExtDirichlet, alhsD, brhsD, idiagD, lmD, hmTTG, p)
        
        error = dsqrt(dot_product(brhsD,brhsD))/neqD
        write(*,*) "erro do passo mecanico: ", error
        if (error < tolNewton .and. j > 1) then
            converged = .true.
            exit
        end if

        ! solve using LDU decomposition, and store the result in the right array
        call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'full')
        call btod(idDesloc, dDis, brhsD, ndofD, nnp)

        ! add correction da to the incremental displacement vector
        u = u + dDis

        !updates the plastic strain, and the stress at each gauss point
        call pos4plast(x, conecNodaisElem, u, strainP, prevStrainP, stress, stressS, trStrainP, stressTotal, hmTTG, biotP, elementIsPlast, p, pInit, isUndrained)

        !computes the internal force
        call calcInternalForce(x, conecNodaisElem, nen, nel, nnp, nsd, stress, fIntJ)
    end do
    if (converged.eqv..false.) write(*,*) 'Newton method not converged.'


1100 format ("u(", I1,",",I1")")
1200 format ("u(", I1,",",I2")")
1300 format ("u(", I2,",",I1")")
1400 format ("u(", I2,",",I2,")")

    end subroutine incrementMechanicPlasticSolution
    !**************************************************************************************
    !**************************************************************************************
    subroutine calcInternalForce(x, conecNodaisElem, nen, nel, nnp, nsd, stress, fIntJ)
    !function import
    use mFuncoesDeForma, only: shgq, shlq
    use mMalha, only: local
    
    !variables import
    use mGlobaisEscalares, only: nrowsh
    use mGlobaisEscalares, only: ibbar
    
    implicit none
    !variables input
    integer :: conecNodaisElem(nen,nel), nen, nel, nnp, nsd
    real*8 :: x(nsd, nnp)
    real*8 :: stress(nrowb,nintD, nel), fIntJ(ndofD, nnp)
    
    ! variables
    real*8 :: xl(nesd, nen), disL(ned2, nen)
    real*8 :: shG(nrowsh,nen,nintD), shL(nrowsh,nen,nintD) !global shape function with derivative, local shape function with derivative
    real*8 :: det(nintD), w(nintD) ! jacobian determinant, gauss integration weight
    real*8 :: shgbr(nrowsh, nen), bbarJ(nrowb, nesD)
    real*8, allocatable :: fIntLoc(:,:)
    real*8, allocatable :: tempVect(:,:)
    
    real*8 :: r(nintD)
    real*8 :: cl
    
    logical :: quad
    integer :: curElement, l, i, j
    integer :: gnn
    
    !--------------------------------------------------------------------------------------
    if (.not.allocated(fIntLoc)) allocate(fIntLoc(ndofD, nen))
    if (.not.allocated(tempVect)) allocate(tempVect(ndofD, nen))
    fIntJ = 0.d0
    tempVect = 0.d0
    
    call shlq(shL,w,nintd,nen) ! generation of local shape functions
    
    do curElement = 1, nel ! foreach element
        !localize coordinates and dirichlet b.c.
        call local(conecnodaiselem(1,curElement),x,xl,nen,nsd,nesd)
        call local(conecnodaiselem(1,curElement),dis,disl,nen,ndofd,ned2)
        
        ! check if element has any coalesced nodes
        quad = .true.
        if (nen.eq.4.and.conecNodaisElem(3,curElement).eq.conecNodaisElem(4,curElement)) quad = .false.
        
        ! calculates global derivatives of shape functions and jacobian determinants
        if(nen==3) call shgq  (xl,det,shL,shG,nintD,curElement,quad,nen)
        if(nen==4) call shgq  (xl,det,shL,shG,nintD,curElement,quad,nen)
        
        !.... ..calculate mean values of shape function global derivatives
        !.... ..for mean-dilatational b-bar formulation
        call xmeansh(shgbr,w,det,r,shg,nen,nintd,iopt,nesd,nrowsh)
        
        fIntLoc = 0.d0 ! resets the local internal force
        do l=1, nintD
            cl = w(l)*det(l)
            
            bbarj=0.d0
            do j=1, nen
                call setbb(bbarj,shg(1:nrowsh,j,l),shgbr(1:nrowsh,j), r(l),nrowsh,nrowb,iopt,ibbar)

                !multiply B'\sigma
                call matMulmTn(bbarj, stress(1:nrowb,l,curElement), tempVect(1:ndofD, j), nrowb, nesD, nrowb, 1, ndofD, 1)
                do i = 1, ndofD
                    fIntLoc(i, j) = fIntLoc(i, j) + tempVect(i, j) * cl
                end do
            end do
        end do
        
        !assemble
        do j=1, nen
            gnn = conecNodaisElem(j, curElement) !get global node coordinate
            do i=1, ndofD
                fIntJ(i, gnn) = fIntJ(i, gnn) + fIntLoc(i, j)
            end do
        end do
    end do
    
    end subroutine calcInternalForce
    !**************************************************************************************
    !**************************************************************************************
    subroutine montarEstruturasDadosDeslocamentoSkyline(conecNodaisElem, nen, nel)
    !function imports
    use mMalha, only:formlm
    
    ! variable input
    integer :: nen, nel, conecNodaisElem(nen, nel)
    
    !--------------------------------------------------------------------------------------
    !generation of lm array
    call formlm(idDesloc,conecNodaisElem,lmD,ndofD,ned2,nen,nel)
    
    !compute column heights in global left-hand-side matrix
    idiagD = 0.d0
    call colht(idiagD,lmD,ned2,nen,nel,neqD)
    
    !compute diagonal addresses of left-hand-side matrix
    call diag(idiagD, neqD, nalhsD)
    
    if(.not.allocated(alhsD)) allocate(alhsd(nalhsD))
    if(.not.allocated(brhsD)) allocate(brhsd(neqD))
    
    alhsD = 0.d0
    brhsD = 0.d0
    
    end subroutine montarEstruturasDadosDeslocamentoSkyline
    !**************************************************************************************
    !**************************************************************************************
    function calcTotalStress(efStress, p, biotCoef)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: efStress(nrowB)
    real*8 :: p, biotCoef
    
    !variables
    real*8 :: calcTotalStress(nrowb)
    real*8 :: identI(nrowb)
    integer :: k
    
    !--------------------------------------------------------------------------------------
    identI = (/ 1., 1., 0., 1. /)
    
    do k = 1, nrowb
        calcTotalStress(k) = efStress(k) - biotCoef * p * identI(k)
    end do
    
    end function calcTotalStress
    !**************************************************************************************
    !**************************************************************************************
    function calcDPPlasticStress(efStress, p, biotCoef)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: efStress(nrowB)
    real*8 :: p, biotCoef
    
    !variables
    real*8 :: calcDPPlasticStress(nrowb)
    real*8 :: identI(nrowb)
    integer :: k
    
    !--------------------------------------------------------------------------------------
    identI = (/ 1., 1., 0., 1. /)
    
    ! terzagui, i.e., plastic incompressibility.
    do k = 1, nrowb
        !calcDPPlasticStress(k) = (efStress(k) + (1 - biotCoef) * p * identI(k))/ (1 + p/h)
        calcDPPlasticStress(k) = (efStress(k) + (1 - biotCoef) * p * identI(k))
    end do
    
    end function calcDPPlasticStress
    !**************************************************************************************
    !**************************************************************************************
    function calcInvariantI1(stress)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: stress(nrowb)
    
    !variables
    real*8 :: calcInvariantI1
    
    !--------------------------------------------------------------------------------------
    calcInvariantI1 = stress(1) + stress(2) + stress(4)
    
    end function calcInvariantI1
    !**************************************************************************************
    !**************************************************************************************
    function calcDevInvariantJ2(devStress)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: devStress(nrowb)
    
    !variables
    real*8 :: calcDevInvariantJ2
    
    !--------------------------------------------------------------------------------------
    calcDevInvariantJ2 = (devStress(1)*devStress(1)+devStress(2)*devStress(2) + devStress(4)*devStress(4) + 2.0D0*devStress(3)*devStress(3)) / 2.d0
    
    end function calcDevInvariantJ2
    !**************************************************************************************
    !**************************************************************************************
    end module mGeomecanica