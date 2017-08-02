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
    module mGlobaisArranjos

    integer, allocatable :: npar(:)
    real*8               :: etime(6)
    character(len=20)   :: title

    real*8,  allocatable :: uTempoN(:)

    integer, allocatable :: mat(:)
    real*8, allocatable  :: grav(:), bf(:,:), c(:,:)

    real*8, allocatable  :: beta(:)
    logical :: listaSolverDisponivel(3)

    end module mGlobaisArranjos
    !
    !
    module mGlobaisEscalares

    integer :: ntype
    integer :: iprtin
    integer :: numParElem=15
    real*8, parameter  :: zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0
    real*8, parameter  :: four=4.0d0, five=5.0d0, six=6.0d0
    real*8, parameter  :: pt1667=0.1666666666666667d0, pt25=0.25d0, pt5=0.5d0
    real*8  :: coef

    integer :: numat
    integer :: nrowsh,npint

    integer :: nvel,nnp,ns

    real(8) :: tzero,tTransporte,tc,tt,dtBlocoTransp
    real(8) :: tempoNucTrans

    real(8) :: tempoTotalVelocidade,tempoTotalPressao
    real(8) :: tempoTotalTransporte,tempoTotalGeomecanica
    !
    real*8  :: tempoMontagemVel, tempoSolverVel
    real*8  :: tempoMontagemGeo, tempoSolverGeo

    logical :: ligarBlocosHetBeta=.false.
    integer :: iflag_beta
    logical :: novaMalha

    CHARACTER(LEN=6)  :: SPLITT
    REAL(8)           :: YEARINJ, S3DIM
    CHARACTER*25      :: TypeProcess
    INTEGER           :: NITGEO, NUMDX, NCREEP
    INTEGER           :: NITHIDRO, MAXITERC, IBBAR
    LOGICAL           :: INITS3, SALTCREEP, LDRAINED
    REAL(8)           :: TOLSIGMA, TOLVELOC
    real*8 :: tolYield, toleplas

    !
    end module mGlobaisEscalares

