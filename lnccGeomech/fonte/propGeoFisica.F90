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
    !=============================================================================
    !
    module mPropGeoFisica
    !
    integer :: iflag_linear
    real(8) :: xmio,xmiw,sro,srw
    real(8) :: rhow, rhoo
    real(8) :: rmi
    integer :: nsw, nso
    real(8) :: xcbloco,ycbloco,zcbloco,xlbloco,ylbloco,zlbloco
    real(8) :: sbloco,sinicial
    real(8) :: sinj
    real(8) :: gf1,gf2,gf3
    real(8) :: eps_df=1.0d+0
    real(8) :: zer_df=1.0d-2
    !
    real*8, allocatable :: perm(:,:), charLenPorNetworkSquared(:,:)
    real*8, allocatable  :: poro0(:), poroL(:), poroE(:)
    real*8, allocatable  :: YOUNG(:)
    real*8, allocatable  :: MASCN(:),MASCN0(:)
    !
    CHARACTER*12, ALLOCATABLE :: GEOFORM(:)
    REAL(8), ALLOCATABLE :: PWELL(:)
    
    real*8 :: permInicial(3) !2d exclusive, todo: deixar gen�rico para dimens�o
    !
    real(8) :: xcbloco_perm,ycbloco_perm,zcbloco_perm
    real(8) :: xlbloco_perm,ylbloco_perm,zlbloco_perm
    !
    real(8) :: xcbloco_YNG, ycbloco_YNG, zcbloco_YNG
    real(8) :: xlbloco_YNG, ylbloco_YNG, zlbloco_YNG
    REAL(8) :: YNGINICIAL, YNGBLOCO
    !
    real(8) :: xcbloco_phi,ycbloco_phi,zcbloco_phi
    real(8) :: xlbloco_phi,ylbloco_phi,zlbloco_phi
    real(8) :: phibloco,phiinicial
    !
    !     Descricao: dados dos campos randomicos
    !
    integer :: nr
    !
    integer :: ncont_mix,ncont_prod
    integer :: iflag_mix,iflag_prod
    integer :: iflag_read_phi,iflag_read_perm
    integer :: IFLAG_READ_YNG,IFLAG_READ_YNG2,IFLAG_READ_YNG3

    !
    real(8) :: tprt_mix,tprt_prod,dtprt_mix,dtprt_prod
    real(8) :: tprt_prodF, dtprt_prodF, tprt_conc, dtprt_conc
    integer :: np_rand_prod,np_rand_mix,npmix,npprod,np_rand_conc
    integer :: np_rand_prodF, ninit_prodF, npprodF, npconc
    !
    character(len=128) :: perm_inx,perm_iny,perm_inz,phi_in, beta_in
    character(len=128) :: mixing_out,prod_out
    CHARACTER(len=128) :: YNG_IN,YNG2_IN,YNG3_IN
    !
    real(8) :: kg,kgphi, kgbeta    ! media geometrica
    real(8) :: rho,rhophi, rhobeta ! coeficiente da variancia (strenght)
    real(8) :: KGYNG, RHOYNG
    real(8) :: KGYNG2, RHOYNG2
    real(8) :: KGYNG3, RHOYNG3
    real(8) :: YNGX1,YNGX2
    real(8) :: YNG3X1,YNG3X2
    !
    integer :: nkozenycarman,normalphi
    real(8) :: sup_KC,const_KC

    !
    !.... Geometria
    !
    integer :: nelx, nely, nelz
    integer :: nelxReserv, nelyReserv, nelzReserv
    real(8) :: dimx, dimy, dimz
    real(8) :: hx, hy, hz
    !
    !.... Geomecanica
    !
    REAL*8 :: CREEPZERO, STRSZERO, XTERLOAD
    REAL*8 :: DTCREEP, POWEREPS, XTEST, SEADEPTH
    !
    REAL(8) :: SIGMAREF, POWERN, TOLCREEP
    REAL(8) :: BULKWATER, BULKOIL
    !
    REAL(8), DIMENSION(9) :: POISVECT
    REAL(8), DIMENSION(9) :: YUNGVECT
    REAL(8), DIMENSION(9) :: RHODVECT
    REAL(8), DIMENSION(9) :: LAMBFACT
    REAL(8), DIMENSION(9) :: GRAINBLK
    REAL(8), DIMENSION(9) :: PORELAGR
    REAL(8), DIMENSION(9) :: MEANRHOW
    REAL(8), DIMENSION(9) :: MEANRHOO
    REAL(8), DIMENSION(9) :: MEANSATR
    REAL(8), DIMENSION(9) :: MEANBLKW
    REAL(8), DIMENSION(9) :: MEANBLKO
    REAL(8), DIMENSION(9) :: MEANBULK
    REAL(8), DIMENSION(9) :: MEANDENS
    real*8 :: waterViscosity(9)
    real*8 :: misesYield(9)
    real*8 :: mcFriction(9), mcC(9)
    real*8 :: dpK(9), dpAlpha(9)
    !
    CHARACTER(5), DIMENSION(9) :: GEOMECLAW

    REAL(8) :: XLEFT,XRIGHT,YTOP,YBOTTOM
    !
    !funcoes e subrotinas
    public :: calcdim, readperm, mapeando
    public :: atribuirPhiBloco, readphi
    public :: estatistica,  xkkc, xlt, xlw,  xlo, xkrw, xkro
    public :: GEOINDIC, GEOYLOC, GEOYLOC3D, FNCMECLAW
    public :: rhoCell
    !
    contains
    !
    !=======================================================================
    !
    SUBROUTINE MAPEIARESERVATORIO(permk,perm,conecNodaisElem,x,xlx,xly,xlz, &
        nelemx,nelemy,nelemz,nelem,nen,nsd)
    !
    !     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
    !
    use mMalha,            only: numel
    !
    implicit none
    !
    integer :: nelemx,nelemy,nelemz,nelem,nen,nsd
    integer,dimension(nen,*)  :: conecNodaisElem
    real(8),dimension(nsd,*)  :: x
    real(8),dimension(*)      :: perm,permk
    !
    integer :: nel,i,NBLOCKX,NBLOCKY
    real(8) :: xlx,xly,xlz,TOL
    real(8) :: cdx,cdy,cdz,aux,aux1,aux2,X1,X2,Y1,Y2
    REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP,CENTRO
    REAL(8) :: XPOS,YPOS
    !
    TOL=1E-6
    XPOS=0.0D0
    YPOS=0.0D0
    !
    if(nsd==2) then
        nelzReserv=1
        nelemz=1
    endif
    !
    !     definicacao do tamanlhos dos elementos
    !     Campo de permeabilidades
    !
    cdx=xlx/nelemx
    cdy=xly/nelemy
    if(nsd==3)cdz=xlz/nelemz
    !
    !     Impressao dos parametros
    !
    IF(nsd==2)THEN
        write(*,800)dimx,dimy,nelxReserv,nelyReserv
        write(*,890)xlx,xly,nelemx,nelemy
        write(*,891)hx,hy,cdx,cdy
    ELSE
        write(*,900)dimx,dimy,dimz,nelxReserv,nelyReserv,nelzReserv
        write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
        write(*,991)hx,hy,hz,cdx,cdy,cdz
    END IF
    !
    !     verificacao da consistencia
    !
    aux = (hx-cdx)
    aux1= (hy-cdy)
    if(nsd==3) aux2= (hz-cdz)
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    end if
    !
    aux = (dimx-xlx)
    aux1= (dimy-xly)
    if(nsd==3) aux2= (dimz-xlz)
    tol = 1.0e-08
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,111)
            stop 305
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol.or.aux2.gt.tol)then
            write(*,111)
            stop 305
        end if
    endif
    CALL ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
    !
    !     loop nos elementos
    !
    nelem=nelxReserv*nelyReserv
    !
    DO NEL=1,NUMEL
        X1=X(1,conecNodaisElem(1,NEL))+TOL
        X2=X(1,conecNodaisElem(3,NEL))-TOL
        Y1=X(2,conecNodaisElem(1,NEL))+TOL
        Y2=X(2,conecNodaisElem(3,NEL))-TOL
        IF(X1.GE.POINTINF(1))THEN
            IF(X2.LE.POINTSUP(1))THEN
                IF(Y1.GE.POINTINF(2))THEN
                    IF(Y2.LE.POINTSUP(2))THEN
                        CENTRO(1)=(X2+X1)*0.5D0
                        CENTRO(2)=(Y1+Y2)*0.5D0
                        I=0
                        DO NBLOCKY=1,nelemy
                            YPOS = POINTINF(2)+REAL(NBLOCKY-1)*CDY
                            DO NBLOCKX=1,nelemx
                                I=I+1
                                XPOS = POINTINF(1)+REAL(NBLOCKX-1)*CDX
                                IF((CENTRO(1).GE.XPOS).AND.(CENTRO(1).LE.(XPOS+CDX)))THEN
                                    IF((CENTRO(2).GE.YPOS).AND.(CENTRO(2).LE.(YPOS+CDY)))THEN
                                        PERMK(NEL) = PERM(I)
                                        exit
                                    END IF
                                END IF
                            END DO
                        END DO
                    END IF
                END IF
            END IF
        END IF
    END DO
    !      STOP
    !
    !
800 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f10.5,2x,'Ly=',f10.5,/,             &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
        '##############################')
    !
890 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f10.5,2x,'Ly=',f10.5,/,               &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
        '##############################')
    !
891 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,/,&
        '##############################',/)

900 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
        '##############################')
    !
990 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
        '##############################')
    !
991 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
        '##############################',/)

110 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
        'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
        'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
        '####################################',/)

111 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS       ',/, &
        'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
        'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
        '####################################',/)
    !
    return
    !
    END SUBROUTINE MAPEIARESERVATORIO

    !
    !=======================================================================
    !
    SUBROUTINE ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
    !
    USE mMalha, only: x,nsd, conecNodaisElem
    USE mMalha, only: numelReserv
    !     USE mPropGeoFisica, only: GEOFORM
    IMPLICIT NONE
    !
    INTEGER :: NO,NELEMI,NELEMF
    REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP
    REAL(8) :: AUXMAX,AUXMIN
    REAL(8) :: DIMRESERVX,DIMRESERVY

    !
    AUXMAX=-1E20
    AUXMIN=1E20
    !
    NELEMI=1
    NELEMF=numelReserv
    !
    NO=conecNodaisElem(1,NELEMI)
    POINTINF(1)=X(1,NO)
    POINTINF(2)=X(2,NO)
    NO=conecNodaisElem(3,NELEMF)
    POINTSUP(1)=X(1,NO)
    POINTSUP(2)=X(2,NO)
    !
    DIMRESERVX=POINTSUP(1)-POINTINF(1)
    DIMRESERVY=POINTSUP(2)-POINTINF(2)
    WRITE(*,*)'COORDENADA INICIAL DO RESERVATORIO:',POINTINF
    WRITE(*,*)'COORDENADA FINAL DO RESERVATORIO..:',POINTSUP
    WRITE(*,*)'DIMENSOES DO RESERVATORIO.........:',DIMRESERVX,DIMRESERVY
    !
    END SUBROUTINE ENCONTRARLIMITESRESERVATORIO
    !
    !=======================================================================
    !
    SUBROUTINE MAPEIADOMINIO(permk,perm,conecNodaisElem,x,xlx,xly,xlz, &
        nelemx,nelemy,nelemz,nelem,nen,nsd,YX1,YX2)
    !
    !     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
    !
    use mMalha,            only: numel
    !       use MMCMC,             only: XLEFT,XRIGHT,YTOP,YBOTTOM
    !
    implicit none
    !
    integer :: nelemx,nelemy,nelemz,nelem,nen,nsd
    integer,dimension(nen,*)  :: conecNodaisElem
    real(8),dimension(nsd,*)  :: x
    real(8),dimension(*)      :: perm,permk
    !
    integer :: nel,j,k,NBLOCK
    real(8) :: xlx,xly,xlz,TOL
    real(8) :: cdx,cdy,cdz,aux,aux1,aux2,X1,X2,Y1,Y2
    REAL(8),DIMENSION(NSD) :: POINTINF,POINTSUP,CENTRO
    REAL(8) :: XPOS,YPOS,AREA,YX1,YX2
    !
    ! PONTO INICIAL SERA A POSICAO (XLEFT,YTOP)
    !      WRITE(*,*)'()()()MAPEANDO()()()'
    !      WRITE(*,*) 'XLEFT = ', XLEFT,'XRIGHT  = ',XRIGHT
    !      WRITE(*,*) 'YTOP  = ',  YTOP,'YBOTTOM = ',YBOTTOM
    !
    TOL=1E-6
    XPOS=0.0D0
    YPOS=0.0D0
    !
    if(nsd==2) then
        nelzReserv=1
        nelemz=1
    endif
    !
    !     definicacao do tamanlhos dos elementos
    !     Campo de permeabilidades
    !
    cdx=xlx/nelemx
    cdy=xly/nelemy
    if(nsd==3)cdz=xlz/nelemz
    !
    !     Impressao dos parametros
    !
    IF(nsd==2)THEN
        !         write(*,800)dimx,dimy,nelxReserv,nelyReserv
        !         write(*,890)xlx,xly,nelemx,nelemy
        write(*,891)hx,hy,cdx,cdy
    ELSE
        write(*,900)dimx,dimy,dimz,nelxReserv,nelyReserv,nelzReserv
        write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
        write(*,991)hx,hy,hz,cdx,cdy,cdz
    END IF
    !
    !     verificacao da consistencia
    !
    AUX = -10.0
    AUX1= -10.0
    DO NEL=1,NUMEL
        X1=X(1,conecNodaisElem(1,NEL))
        X2=X(1,conecNodaisElem(3,NEL))
        Y1=X(2,conecNodaisElem(1,NEL))
        Y2=X(2,conecNodaisElem(3,NEL))
        IF(AUX.LT.ABS(X2-X1))AUX=ABS(X2-X1)
        IF(AUX1.LT.ABS(Y2-Y1))AUX1=ABS(Y2-Y1)
    END DO
    WRITE(*,*)'###################################'
    WRITE(*,*)'###################################'
    WRITE(*,*)'## MAPEANDO MODULO DE YOUNG #######'
    WRITE(*,*)'###################################'
    WRITE(*,*)'### PONTO INICIAL DE MAPEAMENTO ###'
    WRITE(*,*)'###################################'
    WRITE(*,*)'POSICAO X:',YX1
    WRITE(*,*)'POSICAO Y:',YX2
    WRITE(*,*)'###################################'
    WRITE(*,*)'DIMENSAO MAXIMA DE UM ELEMENTO EM X'
    WRITE(*,*)AUX
    WRITE(*,*)'DIMENSAO MAXIMA DE UM ELEMENTO EM Y'
    WRITE(*,*)AUX1
    WRITE(*,*)'###################################'
    WRITE(*,*)'###################################'
    !
    aux = (AUX-cdx)
    aux1= (AUX1-cdy)
    if(nsd==3) aux2= (hz-cdz)
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    end if
    !
    aux = (dimx-xlx)
    aux1= (dimy-xly)
    if(nsd==3) aux2= (dimz-xlz)
    tol = 1.0e-08
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,111)
            stop 305
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol.or.aux2.gt.tol)then
            write(*,111)
            stop 305
        end if
    endif
    !      CALL ENCONTRARLIMITESRESERVATORIO(POINTINF,POINTSUP)
    POINTINF(1)=XLEFT
    POINTINF(2)=YBOTTOM
    POINTSUP(1)=XRIGHT
    POINTSUP(2)=YTOP
    !
    !     loop nos elementos
    !
    nelem=nelemx*nelemy
    !
    DO NEL=1,NUMEL
        AREA=(X(1,conecNodaisElem(4,NEL))*X(2,conecNodaisElem(1,NEL))- &
            X(1,conecNodaisElem(1,NEL))*X(2,conecNodaisElem(4,NEL)))
        DO J=1,NEN-1
            AREA = AREA+(X(1,conecNodaisElem(J,NEL))*X(2,conecNodaisElem(J+1,NEL))- &
                X(1,conecNodaisElem(J+1,NEL))*X(2,conecNodaisElem(J,NEL)))
        END DO
        AREA=AREA*0.5
        AREA=1.0D0/(6.0*AREA)
        AUX=(X(1,conecNodaisElem(nen,NEL))*X(2,conecNodaisElem(1,NEL))-  &
            X(1,conecNodaisElem(1,NEL))*X(2,conecNodaisElem(nen,NEL)))
        CENTRO(1)=(X(1,conecNodaisElem(nen,NEL))+X(1,conecNodaisElem(1,NEL)))*AUX
        CENTRO(2)=(X(2,conecNodaisElem(nen,NEL))+X(2,conecNodaisElem(1,NEL)))*AUX
        DO J=1,NEN-1
            AUX = (X(1,conecNodaisElem(J,NEL))*X(2,conecNodaisElem(J+1,NEL))-  &
                X(1,conecNodaisElem(J+1,NEL))*X(2,conecNodaisElem(J,NEL)))
            CENTRO(1) = CENTRO(1)+ &
                (X(1,conecNodaisElem(J,NEL))+X(1,conecNodaisElem(J+1,NEL)))*AUX
            aux1 = (X(1,conecNodaisElem(J,NEL))+X(1,conecNodaisElem(J+1,NEL)))
            CENTRO(2) = CENTRO(2)+ &
                (X(2,conecNodaisElem(J,NEL))+X(2,conecNodaisElem(J+1,NEL)))*AUX
        END DO
        CENTRO=CENTRO*AREA
        NBLOCK=0
        DO K=1,NELEMY
            YPOS=YX2+REAL(K-1)*CDY
            DO J=1,NELEMX
                NBLOCK=NBLOCK+1
                XPOS=YX1+REAL(J-1)*CDX
                IF(CENTRO(1).GE.XPOS)THEN
                    IF(CENTRO(1).LE.XPOS+CDX)THEN
                        IF(CENTRO(2).GE.YPOS)THEN
                            IF(CENTRO(2).LE.YPOS+CDY)THEN
                                PERMK(NEL) = PERM(NBLOCK)
                                exit
                            END IF
                        END IF
                    END IF
                END IF
            END DO
        END DO
    END DO
    !
800 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f10.5,2x,'Ly=',f10.5,/,             &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
        '##############################')
    !
890 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f10.5,2x,'Ly=',f10.5,/,               &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
        '##############################')
    !
891 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,/,&
        '##############################',/)

900 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
        '##############################')
    !
990 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
        '##############################')
    !
991 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
        '##############################',/)
    !
100 FORMAT('POSICAO X:',F10.5,/,'POSICAO Y:',F10.5)
    !
110 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
        'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
        'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
        '####################################',/)
    !
111 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS       ',/, &
        'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
        'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
        '####################################',/)
    !
    return
    !
    END SUBROUTINE MAPEIADOMINIO
    !
    !=======================================================================
    !
    FUNCTION xkkcGeo(phi,perm)
    !
    use mGlobaisEscalares, only: TypeProcess
    !      use mPropGeoFisica, only: nkozenycarman,sup_KC,const_KC
    !
    IMPLICIT NONE
    !
    !.... Permeabilidade de Kozeny-Carman
    !
    LOGICAL :: LFLAG1, LFLAG2
    REAL(8) :: XKKCGEO, PHI, perm
    REAL(8) :: CS
    !
    LFLAG1 = .FALSE.
    LFLAG2 = .FALSE.
    !
    CS = sup_KC
    !
    IF (TypeProcess.EQ.'TERZAGHI')  LFLAG1 = .TRUE.
    IF (TypeProcess.EQ.'MANDEL')    LFLAG2 = .TRUE.
    !
    IF (LFLAG1.OR.LFLAG2) THEN
        XKKCGEO = perm
    ELSE
        IF(nkozenycarman.eq.1)THEN
            XKKCGEO = (PHI**3)/(CS*CS*(1.0D0-PHI)**2)
        ELSE
            XKKCGEO = PERM
        END IF
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !=======================================================================
    !
    function xkkc(phi)
    implicit none
    !
    !.... Permeabilidade de Kozeny-Carman
    !
    real(8) :: xkkc,phi
    real(8) :: c0,cs
    !
    c0=1.d+2
    cs=1.d0
    !
    xkkc=1.d0*c0*phi**3/((1.d0-phi)*cs)**2
    !
    end function
    !
    !=======================================================================
    !
    function xlt(uu)
    !
    !     calcula a mobilidade total
    !
    implicit none
    real(8) :: xlt,uu
    !
    !     caso linear
    !
    !      xlt=1.d0
    !
    !       xlt=xlw(uu)+xlo(uu)
    !
    xlt=xkrw(uu)/xmiw+xkro(uu)/xmio
    !
    end function
    !
    !=======================================================================
    !
    function xltGeo(uu)
    !
    !     calcula a mobilidade total
    !
    implicit none
    real(8) :: xltGeo,uu
    !
    !     caso linear
    !
    !      xlt=1.d0
    !
    xltGeo=xlw(uu)+xlo(uu)
    !
    !       xlt=xkrw(uu)/xmiw+xkro(uu)/xmio
    !
    end function
    !
    !=======================================================================
    !
    function xlw(uu)
    !
    !     calcula a mobilidade da agua
    !
    implicit none
    real(8) :: xlw,uu
    !
    !.... Agua-Oleo
    !
    !     mobilidades
    !
    xlw=xkrw(uu)/xmiw
    !
    end function
    !
    !=======================================================================
    !
    function xlo(uu)
    !
    !     calcula a mobilidade total
    !
    implicit none
    real(8) :: xlo,uu
    !
    !.... Agua-Oleo
    !
    !     mobilidades
    !
    xlo=xkro(uu)/xmio
    !
    end function
    !
    !=======================================================================
    !
    function xkrw(uu)
    !
    !     calcula a permeabilidade relativa da agua
    !
    implicit none
    !
    real(8) :: xkrw,uu
    real(8) :: uus,zero
    !
    zero=0.d0
    xkrw=0.d0
    !
    !.... Permeabilidades relativas
    !
    !     Agua-Oleo
    !
    uus=uu-srw
    if(uus.lt.zero) uus=zero
    !      write(*,*) ' iflag_linear = ',  iflag_linear
    select case(iflag_linear)
    case(1)
        xkrw=uus/(1.d0-srw)
    case(2)
        xkrw=uus**2/(1.d0-srw)**nsw
    end select
    !
    end function
    !
    !=======================================================================
    !
    function xkro(uu)
    !
    !     calcula a permeavilidade relativa do oleo
    !
    implicit none
    !
    real(8) :: xkro,uu
    real(8) :: uus,zero
    !
    zero=0.d0
    xkro=0.d0
    !
    !.... Permeabilidades relativas
    !
    !     Agua-Oleo
    !
    uus=1.d0-sro-uu
    !      write(*,*) ' iflag_linear = ',  iflag_linear
    if(uus.lt.zero) uus=zero
    select case(iflag_linear)
    case(1)
        xkro=uus/(1.d0-sro)
    case(2)
        xkro=uus**2/(1.d0-sro)**nso
    end select
    !
    end function
    !
    !======================================================================
    !
    subroutine calcdim(nsd,numnp,x)
    !
    implicit none
    !
    integer :: nsd,numnp
    real(8), dimension(nsd,*) :: x
    !
    dimx=dabs(x(1,numnp)-x(1,1))
    dimy=dabs(x(2,numnp)-x(2,1))
    if(nsd==3)dimz=dabs(x(3,numnp)-x(3,1))

    hx  =dimx/nelxReserv
    hy  =dimy/nelyReserv
    if(nsd==3)hz  =dimz/nelzReserv

    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine lerPropriedadesFisicas()
    use mMalha, only: numel,nsd
    use mMalha, only: numelReserv
    !
    implicit none
    !
    integer :: tensDim
    integer :: i, j

    if (nsd == 2) then
        tensDim = 3
    else
        tensDim = 6
    end if
    
    !     cria o vetor de permeabilidades
    if(.not.allocated(perm)) allocate(perm(tensDim,numelReserv));
    if(.not.allocated(charLenPorNetworkSquared)) allocate(charLenPorNetworkSquared(tensDim,numelReserv));
    if(.not.allocated(poro0)) allocate(poro0(numelReserv));
    if(.not.allocated(poroL)) allocate(poroL(numelReserv));
    if(.not.allocated(poroE)) allocate(poroE(numelReserv));
    !
    if(.not.allocated(YOUNG)) allocate(YOUNG(NUMEL));          YOUNG    = 0.0D0
    IF (NSD.EQ.2) then
        if(.not.allocated(PWELL)) ALLOCATE(PWELL(nelYReserv)); PWELL    = 0.0D0
    ELSE
        if(.not.allocated(PWELL)) ALLOCATE(PWELL(nelZReserv)); PWELL    = 0.0D0
    ENDIF
    
    do j=1,numelReserv
        poro0(j) = PORELAGR(1)
        poroL(j) = PORELAGR(1)
        poroE(j) = PORELAGR(1)
        
        do i=1,tensDim
            perm(i,j) = permInicial(i)
            charLenPorNetworkSquared(i,j) = permInicial(i) / calcKozenyCarmanCoefficient(poro0(j))
        end do
        
        
    end do
    !
    end subroutine lerPropriedadesFisicas
    !
    !=======================================================================
    !
    subroutine readperm(fname,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
    !
    !     Objetivo: le as permeabilidades de um arquivo
    !
    !----------------------------------------------------------------------
    !
    implicit none
    !
    real(8), dimension(*) :: perm
    integer :: nsd
    !
    integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j,k
    real(8) :: xlx,xly,xlz,beta
    character(len=128) :: NAME, FILE_IN
    character(len=128) :: fname
    character(len=4)   :: EXT,tipo
    character(len=5)   :: C
    integer :: cont, contK, fim, numel_
    !
    !
    inperm=201
    ! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    !
    file_in=fname
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
    open(inperm, file= NAME)
    !
    !     dimensoes do dominio
    !
    read(inperm,*) xlx
    read(inperm,*) xly
    if(nsd==3) read(inperm,*) xlz
    !
    !     numero de elementos em cada direcao
    !
    read(inperm,*) nelemx
    read(inperm,*) nelemy
    if(nsd==3)read(inperm,*) nelemz
    !
    !     verificacao do tamanho do vetor para leitura
    !
    if(nsd==2) then
        if(nelemx*nelemy.gt.nelxReserv*nelyReserv*100)then
            write(*,115)(nr-1)
            stop 400
        end if
    else
        if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
            write(*,115)(nr-1)
            stop 400
        end if
    endif
    !
    !     ntype = 1; campo exponencial
    !     ntype = 2; campo fractal
    !
    read(inperm,*) ntype
    !
    !     beta: coeficiente de Hurst
    !
    read(inperm,*) beta
    !
    !     leituras vazias
    !
    read(inperm,*)
    read(inperm,*)
    !
    !     inicio da leitura do campo
    !
    if(nsd==2) fim=1
    if(nsd==3) fim=nelemz
    contK=0
    !
    do k=1,fim
        !
        contK=1+(k-1)*nelemx*nelemy
        !
        if(nsd==3) then
            read(inperm,*) nline
            if(nline+1.ne.k) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop 225
            end if
        endif
        !
        do j=1,nelemy
            !
            read(inperm,*) nline
            if(nline+1.ne.j) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop 225
            end if
            !
            cont=contK+(j-1)*nelemx
            read(inperm,*) (perm(i),i=cont,nelemx+cont-1)
            read(inperm,*) nflag
            if(nflag.ne.192837465) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nflag'
                stop 225
            end if
            !
        end do ! nelemy
    end do ! nelemz
    !
    close(inperm)
    !
    !     calculando a permeabilidade lognormal
    !
    if(nsd==2) numel_=nelemy*nelemx
    if(nsd==3) numel_=nelemz*nelemy*nelemx

    do j=1,numel_
        perm(j)=kg*dexp(rho*perm(j))
    end do
    !
    tipo='K'
    call estatistica(perm,nelemx,nelemy,nelemz,tipo)
    !
111 FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113 FORMAT(I5)
115 FORMAT(                                      &
        '######################################',/, &
        'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
        'LEITURA DO CAMPO DE PERMEABILIDADES   ',/, &
        'NUMERO DO CAMPO:',I5,/,                    &
        '######################################',/)
    !
    end subroutine readperm

    !======================================================================
    subroutine estatistica(prm,nmx,nmy,nmz,sc)
    !
    use mMalha, only: nsd
    !     calculando a media e variancia
    implicit none
    real(8), dimension(*) :: prm
    integer :: nmx,nmy,nmz,i,numel_
    real(8) :: xm,xmm,xv,XMAX,XMIN
    character*4 sc
    !
    xm  = 0.d0
    xmm = 0.d0
    xv  = 0.d0
    XMAX = -1.d18
    XMIN = 1.d18
    !
    if(nsd==2) numel_=nmy*nmx
    if(nsd==3) numel_=nmz*nmy*nmx
    do i=1,numel_
        xm = xm+prm(i)
        xmm= xmm+prm(i)*prm(i)
        if(prm(i).gt.XMAX) XMAX = prm(i)
        if(prm(i).lt.XMIN) XMIN = prm(i)
    end do
    xm  = xm/numel_
    xmm = xmm/numel_
    xv  = xmm - xm*xm
    !
    !
    write(*,123)sc,xm,xv,XMAX,XMIN
123 FORMAT(A,                          &
        & ' ######################### ',/,  &
        & 'MEDIA     =',1PE15.8,/,          &
        & 'VARIANCIA =',1PE15.8,/,          &
        & 'MAXIMO    =',1PE15.8,/,          &
        & 'MINIMO    =',1PE15.8,/,          &
        & '##############################')
    !
    end subroutine estatistica
    !
    !===========================================================================
    !
    subroutine mapeando(permk,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
    !
    !     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
    !
    !
    implicit none
    !
    integer :: nelemx,nelemy,nelemz,nelem,nen,nsd
    integer,dimension(nen,*):: conecNodaisElem
    real(8),dimension(nsd,*)  :: x
    real(8),dimension(*)    :: perm,permk
    !
    integer :: nel,i,no,j,k,ncont
    real(8) :: xi,xf,yi,yf,zi,zf,xx,yy,zz,xlx,xly,xlz
    real(8) :: cdx,cdy,cdz,aux,aux1,aux2,tol
    !
    if(nsd==2) then
        nelz=1
        nelemz=1
    endif
    !
    !     definicacao do tamanlhos dos elementos
    !     Campo de permeabilidades
    !
    cdx=xlx/nelemx
    cdy=xly/nelemy
    if(nsd==3)cdz=xlz/nelemz
    !
    !     Impressao dos parametros
    !
    if(nsd==2) then
        write(*,800)dimx,dimy,nelxReserv,nelyReserv
        write(*,890)xlx,xly,nelemx,nelemy
        write(*,891)hx,hy,cdx,cdy
    else
        write(*,900)dimx,dimy,dimz,nelx,nely,nelz
        write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
        write(*,991)hx,hy,hz,cdx,cdy,cdz
    endif
    !
    !     verificacao da consistencia
    !
    aux = (hx-cdx)
    aux1= (hy-cdy)
    if(nsd==3) aux2= (hz-cdz)
    tol = 1.0e-08
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,110)
            stop 300
        end if
    end if
    !
    aux = (dimx-xlx)
    aux1= (dimy-xly)
    if(nsd==3) aux2= (dimz-xlz)
    tol = 1.0e-08
    !
    if(nsd==2) then
        if(aux.gt.tol.or.aux1.gt.tol)then
            write(*,111)
            stop 305
        end if
    else
        if(aux.gt.tol.or.aux1.gt.tol.or.aux2.gt.tol)then
            write(*,111)
            stop 305
        end if
    endif
    !
    !     loop nos elementos
    !
    nelem=nelxReserv*nelyReserv
    !
    if((nelx==nelemx).and.(nely==nelemy).and.(nelz==nelemz)) then
        do nel=1,nelem
            xx=0.d0
            yy=0.d0
            if(nsd==3)zz=0.d0
            do i=1,nen
                no = conecNodaisElem(i,nel)
                xx=xx+x(1,no)
                yy=yy+x(2,no)
                if(nsd==3)zz=zz+x(3,no)
            end do
            xx=xx/nen
            yy=yy/nen
            if(nsd==3)zz=zz/nen
            permk(nel)=perm(nel)
        end do
        return
    end if
    !
    do nel=1,nelem
        !
        !     calculo do centroide do elemento
        !
        xx=0.d0
        yy=0.d0
        if(nsd==3)zz=0.d0
        do i=1,nen
            no = conecNodaisElem(i,nel)
            xx=xx+x(1,no)
            yy=yy+x(2,no)
            if(nsd==3)zz=zz+x(3,no)
        end do
        xx=xx/nen
        yy=yy/nen
        if(nsd==3)zz=zz/nen
        !
        !     mapear o centroide no campo de permeabilidade
        !
        ncont=0d0
        !
        if(nsd==2) then
            !
            ncont=0d0
            do j=1,nelemy
                yi = (j-1)*cdy
                yf = yi+cdy
                !
                do i=1,nelemx
                    ncont=ncont+1
                    xi = (i-1)*cdx
                    xf = xi+cdx
                    if(xx.ge.xi.and.xx.le.xf)then
                        if(yy.ge.yi.and.yy.le.yf)then
                            permk(nel)=perm(ncont)
                            go to 100
                        end if
                    end if
                end do !i
            end do !j
            !
        else
            !
            do k=1,nelemz
                zi = (k-1)*cdz
                zf = zi+cdz
                do j=1,nelemy

                    yi = (j-1)*cdy
                    yf = yi+cdy
                    !
                    do i=1,nelemx
                        ncont=ncont+1
                        xi = (i-1)*cdx
                        xf = xi+cdx
                        if(xx.ge.xi.and.xx.le.xf)then
                            if(yy.ge.yi.and.yy.le.yf)then
                                if(zz.ge.zi.and.zz.le.zf)then
                                    permk(nel)=perm(ncont)
                                    go to 100
                                end if
                            end if
                        end if
                    end do !i
                end do !j
            end do !k

        endif
        !
100     continue
        !
    end do !nel:nelem
    !
800 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f15.7,2x,'Ly=',f15.7,/,             &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
        '##############################')
    !
890 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f15.7,'Ly=',f15.7,/,               &
        'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
        '##############################')
    !
891 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,/,&
        '##############################',/)

900 format(/,                            &
        '##############################',/, &
        'TAMANHO DO RESERVATORIO:',/,'Lx=', &
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
        '##############################')
    !
990 format(/,                              &
        '##############################',/,   &
        'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
        f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
        'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
        '##############################')
    !
991 format(/,                            &
        '##############################',/, &
        'Tamanho do elemento:',/,'hx=',     &
        f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
        '##############################',/, &
        'Tamanho do bloco geologico:',/,'gx=', &
        f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
        '##############################',/)

110 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
        'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
        'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
        '####################################',/)

111 format('####################################',/, &
        'INCONSISTENCIA NO TAMANHO DOS       ',/, &
        'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
        'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
        '####################################',/)
    !
    return
    !
    end subroutine
    
    !
    !======================================================================
    !
    subroutine readphi(perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
    !
    !     Objetivo: le as porosidades de um arquivo
    !
    !----------------------------------------------------------------------
    !
    implicit none
    !
    real(8), dimension(*) :: perm
    integer :: nsd
    !
    integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j,k
    real(8) :: xlx,xly,xlz,beta
    character(len=128) :: NAME, FILE_IN
    character(len=4)   :: EXT, tipo
    character(len=5)   :: C
    REAL(8) :: TOL
    !
    TOL = 1E-2
    !
    if(nsd==2) then
        nelemz=1
        nelz=1
        xlz=1.0
    endif
    !
    !     abre o arquivo de porosidades
    !
    inperm=203
    ! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    FILE_IN=phi_in
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
    open(inperm, file= NAME)
    !
    !     dimensoes do dominio
    !
    read(inperm,*) xlx
    read(inperm,*) xly
    if(nsd==3)read(inperm,*) xlz
    !
    !     numero de elementos em cada direcao
    !
    read(inperm,*) nelemx
    read(inperm,*) nelemy
    if(nsd==3)read(inperm,*) nelemz
    !
    !     verificacao do tamanho do vetor para leitura
    !
    if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
        write(*,115)(nr-1)
        stop 405
    end if
    !
    !     ntype = 1; campo exponencial
    !     ntype = 2; campo fractal
    !
    read(inperm,*) ntype
    !
    !     beta: coeficiente de Hurst
    !
    read(inperm,*) beta
    !
    !     leituras vazias
    !
    read(inperm,*)
    read(inperm,*)
    !
    !     inicio da leitura do campo
    !
    do k=1,nelemz
        !
        if(nsd==3) then
            read(inperm,*) nline
            if(nline+1.ne.k) then
                write(*,*) 'Erro na leitura do campo de porosidade, nline'
                stop 230
            end if
        endif
        !
        do j=1,nelemy
            !
            read(inperm,*) nline
            if(nline+1.ne.j) then
                write(*,*) 'Erro na leitura do campo de porosidade, nline'
                stop 230
            end if
            !
            read(inperm,*) (perm(i+(j-1)*nelemx+(k-1)*nelemx*nelemy),i=1,nelemx)
            !
            read(inperm,*) nflag
            if(nflag.ne.192837465) then
                write(*,*) 'Erro na leitura do campo de porosidade, nflag'
                stop 230
            end if
            !
        end do ! nelemy
    end do ! nelemz
    !
    close(inperm)
    !
    !     calculando a porosidade lognormal
    !
    if(normalphi.eq.1)then
        do j=1,nelemz*nelemy*nelemx
            perm(j)=kgphi+rhophi*perm(j)
            IF(perm(j)<TOL)perm(j)=TOL
        end do
    else
        do j=1,nelemz*nelemy*nelemx
            perm(j)=kgphi*dexp(rhophi*perm(j))
        end do
    end if
    !
    tipo='PHI'
    call estatistica(perm,nelemx,nelemy,nelemz,tipo)
    !
111 FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113 FORMAT(I5)
115 FORMAT(                                     &
        '######################################',/, &
        'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
        'LEITURA DO CAMPO DE POROSIDADES       ',/, &
        'NUMERO DO CAMPO:',I5,/,                    &
        '######################################',/)
    !
    end subroutine readphi
    !
    !=======================================================================
    !
    subroutine atribuirPhiBloco(nsd,numel,u,x)
    !
    implicit none
    !
    integer :: nsd,numel
    real*8 :: x(nsd,*)
    real(8), dimension(*)   :: u
    !
    integer :: i
    real(8) :: xx,yy,zz,xi,xf,yi,yf,zi,zf
    !
    !.... Bloco
    !
    xi=xcbloco_phi-xlbloco_phi/2.d0
    xf=xcbloco_phi+xlbloco_phi/2.d0
    yi=ycbloco_phi-ylbloco_phi/2.d0
    yf=ycbloco_phi+ylbloco_phi/2.d0
    if(nsd==3) zi=zcbloco_phi-zlbloco_phi/2.d0
    if(nsd==3) zf=zcbloco_phi+zlbloco_phi/2.d0


    do i=1,numel
        !
        xx=x(1,i)
        yy=x(2,i)
        if(nsd==3) zz=x(3,i)
        !
        u(i)=phiinicial

        if(nsd==2) then
            if(xx.gt.xi.and.xx.lt.xf) then
                if(yy.gt.yi.and.yy.lt.yf) then
                    u(i)=phibloco
                end if
            end if
        else
            if(xx.gt.xi.and.xx.lt.xf) then
                if(yy.gt.yi.and.yy.lt.yf) then
                    if(zz.gt.zi.and.zz.lt.zf) then
                        u(i)=phibloco
                    end if
                end if
            end if
        endif
        !
    end do
    !
    end subroutine atribuirPhiBloco
    !
    !*** *** NEW: ** READPERM SUBROUTINE MODIFIED FOR YOUNG MODULUS ***********
    SUBROUTINE READYNG(perm,xlx,xly,nelemx,nelemy,YNGIN,RHO,KG)
    !
    !....  READ STOCHASTIC FIELD
    !
    implicit none
    !
    real(8), dimension(*) :: perm
    !
    integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j
    real(8) :: xlx,xly,beta,RHO,KG
    character(len=128) :: NAME, FILE_IN,YNGIN
    character(len=4)   :: EXT , tipo
    character(len=5)   :: C
    !

    nelemz=0
    !     abre o arquivo de permeabilidades
    !
    inperm=204
    ! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    !
    file_in=YNGIN
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
    open(inperm, file= NAME)
    !
    !     dimensoes do dominio
    !
    read(inperm,*) xlx
    read(inperm,*) xly
    !
    !     numero de elementos em cada direcao
    !
    read(inperm,*) nelemx
    read(inperm,*) nelemy
    !
    !     verificacao do tamanho do vetor para leitura
    !
    if(nelemx*nelemy.gt.nelx*nely*100)then
        write(*,115)(nr-1)
        stop 410
    end if
    !
    !     ntype = 1; campo exponencial
    !     ntype = 2; campo fractal
    !
    read(inperm,*) ntype
    !
    !     beta: coeficiente de Hurst
    !
    read(inperm,*) beta
    !
    !     leituras vazias
    !
    read(inperm,*)
    read(inperm,*)
    !
    !     inicio da leitura do campo
    !
    do j=1,nelemy
        !
        read(inperm,*) nline
        if(nline+1.ne.j) then
            write(*,*) 'Erro na leitura do campo de YOUNG, nline'
            stop 235
        end if
        !
        read(inperm,*) (perm(i+(j-1)*nelemx),i=1,nelemx)
        !
        read(inperm,*) nflag
        if(nflag.ne.192837465) then
            write(*,*) 'Erro na leitura do campo de YOUNG, nflag'
            stop 235
        end if
        !
    end do ! nelemy
    !
    close(inperm)
    !     calculando a permeabilidade lognormal
    !
    DO J=1,NELEMY*NELEMX
        PERM(j)=KG*DEXP(RHO*PERM(J))
    END DO
    !
    tipo='EE'
    call estatistica(perm,nelemx,nelemy,nelemz,tipo)
    !
111 FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113 FORMAT(I5)
115 FORMAT( &
        & '######################################',/, &
        & 'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
        & 'LEITURA DO CAMPO DE MODULO DE YOUNG   ',/, &
        & 'NUMERO DO CAMPO:',I5,/, &
        & '######################################',/)
    !
    end SUBROUTINE READYNG
    !
    !
    !======================================================================
    !
    SUBROUTINE READYNGBLOCO1(U,GEOFORM,NUMEL)
    !
    IMPLICIT NONE
    !
    INTEGER :: NEL,NUMEL
    REAL(8),      DIMENSION(NUMEL) :: U
    CHARACTER*12, DIMENSION(NUMEL) :: GEOFORM
    !
    DO 500 NEL=1,NUMEL
        U(NEL) = GEOINDIC('YOUNGMD',GEOFORM(NEL))
500 CONTINUE
    !
    RETURN
    !
    END SUBROUTINE
    !
    !*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
    !
    function bulk(young,poisson,s3dim)
    !
    !... compute bulk modulus of rock bulk modulus
    !
    implicit none
    !
    real(8) :: amu2, alam, young, poisson, bulk, s3dim
    !
    !... stress dimension:  s3dim = 1.0d0 or 3.0d0
    !
    amu2  = young/(1.0d0 + poisson)
    !
    alam  = amu2*poisson/(1.0d0-2.0d0*poisson)
    !
    bulk  = alam + amu2/s3dim
    !
    end function
    !
    !**** NEW ***** FOR SISMIC REPRESENTATION ******************************
    !
    function geoindic(task,geof)
    !
    use mmalha, only: nsd
    !
    !.... function to call mechanical parameters and densities of
    !....         geological formation
    !
    implicit none
    integer      :: i, lastregion
    character*7  :: task
    character*12 :: geof
    character*12, dimension(9) :: refgeo
    real*8                     :: geoloc, geoindic
    !
    data refgeo(1)     ,     refgeo(2),     refgeo(3)/ &
        &     'RESERVATORIO','RIFT_CAPROCK','RIGHT_RESERV'/ &
        &     refgeo(4)     ,     refgeo(5),     refgeo(6)/ &
        &     'SALT_CAPROCK','LEFT__RESERV','POS_SAL_OVER'/ &
        &     refgeo(7)     ,     refgeo(8),     refgeo(9)/ &
        &     'FRONT_RESERV','BASE_CAPROCK','BACK__RESERV'/

    !--------------------------------------------------------------------------------------
    geoindic = 0.0d0
    lastregion = 3*nsd

    select case (task)
    case ('POISSON')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = poisvect(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('YOUNGMD')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = yungvect(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('ROKDENS')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = rhodvect(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('BLKGRIN')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = grainblk(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('POROSTY')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = porelagr(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('MEANSAT')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = meansatr(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('FLUDENS')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = meandens(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('BLKFLUD')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = meanbulk(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('BLKWATR')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = meanblkw(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('MISESYI')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = misesYield(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('MCFRICT')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = mcFriction(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('MCCOHES')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = mcC(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('DPCOHES')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = dpK(i)
                exit
            end if
        end do
        geoindic = geoloc
    case ('DPALPHA')
        do i=1,lastregion
            if (refgeo(i).eq.geof) then
                geoloc = dpAlpha(i)
                exit
            end if
        end do
        geoindic = geoloc
        case default
        
        write(*,*)  "Hmmmm, I don't know"
    end select
    !
    end function
    !
    !**** NEW ***** FOR GEOLOGICAL FORMATION ******************************
    !                        123456789+12
    FUNCTION FNCMECLAW(GEOFORMATION)
    !
    use mMalha, only: NSD
    !
    !.... FUNCTION TO SET STRESS--STRAIN CONSTITUTIVE LAW OF
    !....         GEOLOGICAL FORMATION
    !
    IMPLICIT NONE
    !
    INTEGER      :: I, LASTREGION
    CHARACTER*5  :: GEOLAW, FNCMECLAW
    CHARACTER*12 :: GEOFORMATION
    CHARACTER*12, DIMENSION(9) :: REFGEO
    !
    DATA REFGEO(1)     ,     REFGEO(2),     REFGEO(3)/ &
        &     'RESERVATORIO','RIFT_CAPROCK','RIGHT_RESERV'/ &
        &     REFGEO(4)     ,     REFGEO(5),     REFGEO(6)/ &
        &     'SALT_CAPROCK','LEFT__RESERV','POS_SAL_OVER'/ &
        &     REFGEO(7)     ,     REFGEO(8),     REFGEO(9)/ &
        &     'FRONT_RESERV','BASE_CAPROCK','BACK__RESERV'/
    !
    LASTREGION = 3*NSD
    !
    DO 100 I=1,LASTREGION
        IF (REFGEO(I).EQ.GEOFORMATION) GEOLAW = GEOMECLAW(I)
100 CONTINUE
    !
    FNCMECLAW = GEOLAW
    !
    END FUNCTION
    !
    !**** NEW ***** FOR SISMIC REPRESENTATION ******************************
    !
    FUNCTION GEOYLOC(XLOC,YLOC)
    !
    use mMalha, only: POSTLINE, SALTLINE, RTOPLINE, RBTTLINE, RIFTLINE
    use mMalha, only: LEFTLINE, RGHTLINE, IDome
    !
    !.... FUNCTION TO CALL MECHANICAL PARAMETERS AND DENSITIES OF
    !....         GEOLOGICAL FORMATION
    !
    IMPLICIT NONE
    !
    CHARACTER*12          :: GEOYLOC
    REAL*8                :: YMIN, YMAX
    REAL(8)               :: XLOC, YLOC
    REAL(8), DIMENSION(4) :: YMARK
    !
    YMARK(1) = SALTLINE
    YMARK(2) = RTOPLINE
    YMARK(3) = RBTTLINE
    YMARK(4) = RIFTLINE
    !
    !.... POST-SALT GEOFORMATION:
    !
    YMIN     = POSTLINE
    YMAX     = DOMESURF(XLOC,YMARK(1),IDome) ! YMARK(1)  !
    !
    IF ((YLOC.LT.YMIN).AND.(YLOC.GT.YMAX)) THEN
        GEOYLOC = 'POS_SAL_OVER'
        RETURN
    ENDIF
    !
    !.... SALT GEOFORMATION
    !
    YMIN = YMAX
    YMAX = YMARK(2)
    !
    IF ((YLOC.LT.YMIN).AND.(YLOC.GT.YMAX)) THEN
        GEOYLOC = 'SALT_CAPROCK'
        RETURN
    ENDIF
    !
    !.... SIDE BURDENS OF RESERVOIR GEOFORMATIONS
    !
    YMIN = YMAX
    YMAX = YMARK(3)
    !
    IF ((YLOC.LT.YMIN).AND.(YLOC.GT.YMAX)) THEN
        !.... .. RESERVOIR LEFT SIDE GEOFORMATION
        IF (XLOC.LT.LEFTLINE) THEN
            GEOYLOC = 'LEFT__RESERV'
            RETURN
        ENDIF
        !.... .. RESERVOIR GEOFORMATION
        IF ((XLOC.GE.LEFTLINE).AND.(XLOC.LE.RGHTLINE)) THEN
            GEOYLOC = 'RESERVATORIO'
            RETURN
        ENDIF
        !.... .. RESERVOIR RIGHT SIDE GEOFORMATION
        IF (XLOC.GT.RGHTLINE) THEN
            GEOYLOC = 'RIGHT_RESERV'
            RETURN
        ENDIF
    ENDIF
    !
    !.... RIFT ROCK GEOFORMATION
    !
    YMIN = YMAX
    YMAX = YMARK(4)
    !
    IF ((YLOC.LT.YMIN).AND.(YLOC.GT.YMAX)) THEN
        GEOYLOC = 'RIFT_CAPROCK'
        RETURN
    ENDIF
    !
    END FUNCTION
    !
    !**** NEW ***** FOR SISMIC REPRESENTATION ******************************
    !
    FUNCTION GEOYLOC3D(XLOC,YLOC,ZLOC)
    !
    use mMalha, only: POSTLINE, SALTLINE, RTOPLINE, RBTTLINE, RIFTLINE
    use mMalha, only: LEFTLINE, RGHTLINE, FRNTLINE, BACKLINE, BASELINE
    use mMalha, only: IDome
    !
    !.... FUNCTION TO CALL MECHANICAL PARAMETERS AND DENSITIES OF
    !....         GEOLOGICAL FORMATION
    !
    IMPLICIT NONE
    !
    CHARACTER*12          :: GEOYLOC3D
    REAL*8                :: ZMIN, ZMAX
    REAL(8)               :: XLOC, YLOC, ZLOC
    REAL(8), DIMENSION(5) :: ZMARK
    !
    ZMARK(1) = SALTLINE
    ZMARK(2) = RTOPLINE
    ZMARK(3) = RBTTLINE
    ZMARK(4) = RIFTLINE
    ZMARK(5) = BASELINE
    !
    !.... POST-SALT GEOFORMATION:
    !
    ZMIN     = POSTLINE
    ZMAX     = DOMESURF(XLOC,ZMARK(1),IDome) ! YMARK(1)  !
    !
    IF ((ZLOC.LT.ZMIN).AND.(ZLOC.GT.ZMAX)) THEN
        GEOYLOC3D = 'POS_SAL_OVER'
        RETURN
    ENDIF
    !
    !.... SALT GEOFORMATION
    !
    ZMIN = ZMAX
    ZMAX = ZMARK(2)
    !
    IF ((ZLOC.LT.ZMIN).AND.(ZLOC.GT.ZMAX)) THEN
        GEOYLOC3D = 'SALT_CAPROCK'
        RETURN
    ENDIF
    !
    !.... SIDE BURDENS OF RESERVOIR GEOFORMATIONS
    !
    ZMIN = ZMAX
    ZMAX = ZMARK(3)
    !
    IF ((ZLOC.LT.ZMIN).AND.(ZLOC.GT.ZMAX)) THEN
        IF ((YLOC.GE.BACKLINE).AND.(YLOC.LE.FRNTLINE)) THEN
            !.... ..... RESERVOIR LEFT SIDE GEOFORMATION
            IF (XLOC.LT.LEFTLINE) THEN
                GEOYLOC3D = 'LEFT__RESERV'
                RETURN
            ENDIF
            !.... ..... RESERVOIR GEOFORMATION
            IF ((XLOC.GE.LEFTLINE).AND.(XLOC.LE.RGHTLINE)) THEN
                GEOYLOC3D = 'RESERVATORIO'
                RETURN
            ENDIF
            !.... ..... RESERVOIR RIGHT SIDE GEOFORMATION
            IF (XLOC.GT.RGHTLINE) THEN
                GEOYLOC3D = 'RIGHT_RESERV'
                RETURN
            ENDIF
        ELSE
            !.... ..... RESERVOIR FRONT FACE SIDE GEOFORMATION
            IF (YLOC.GT.FRNTLINE) THEN
                GEOYLOC3D = 'FRONT_RESERV'
                RETURN
            ENDIF
            !.... .. RESERVOIR RIGHT SIDE GEOFORMATION
            IF (YLOC.LT.BACKLINE) THEN
                GEOYLOC3D = 'BACK__RESERV'
                RETURN
            ENDIF
        ENDIF
    ENDIF
    !
    !.... RIFT ROCK GEOFORMATION
    !
    ZMIN = ZMAX
    ZMAX = ZMARK(4)
    !
    IF ((ZLOC.LT.ZMIN).AND.(ZLOC.GT.ZMAX)) THEN
        GEOYLOC3D = 'RIFT_CAPROCK'
        RETURN
    ENDIF
    !
    !.... BASEMENT ROCK GEOFORMATION
    !
    ZMIN = ZMAX
    ZMAX = ZMARK(5)
    !
    IF ((ZLOC.LT.ZMIN).AND.(ZLOC.GT.ZMAX)) THEN
        GEOYLOC3D = 'BASE_CAPROCK'
        RETURN
    ENDIF
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR PLANE SALT TOP BOUNDARY  ***********************
    !
    FUNCTION DOMESURF(X,TRNL,IDome)
    !
    !.... FUNCTION TO DETERMINE TOP BOUNDARY OF SALT FORMATION
    !
    IMPLICIT NONE
    !
    INTEGER     :: IDome
    REAL(8)     :: DOMESURF  !,CURVSURF
    REAL(8)     :: X,TRNL
    !
    IF (IDome.EQ.0) THEN
        DOMESURF = TRNL
    ELSE
        DOMESURF = CURVSURF(X,TRNL)
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION CURVSURF(X,TRNL)
    !
    !.... FUNCTION TO DETERMINE TOP BOUNDARY OF SALT FORMATION
    !....   THE POINTS WERE OBTAINED FROM TRAMOS.DAT
    !
    IMPLICIT NONE
    !
    REAL(8) :: CURVSURF
    REAL(8) :: X,TRNL
    !
    IF (X.LT.0.0D0) THEN
        CURVSURF = TRAMO0(TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GE.0.0D0).AND.(X.LE.900.0D0)) THEN
        CURVSURF = TRAMO1(X,TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GE.900.0D0).AND.(X.LE.1260.0D0)) THEN
        CURVSURF = TRAMO2(X,TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GT.1260.0D0).AND.(X.LE.4320.0D0)) THEN
        CURVSURF = TRAMO3(X,TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GT.4320.0D0).AND.(X.LT.7380.0D0)) THEN
        CURVSURF = TRAMO4(X,TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GE.7380.0D0).AND.(X.LE.9540.0D0)) THEN
        CURVSURF = TRAMO5(TRNL)
        RETURN
    ENDIF
    !
    IF ((X.GT.9540.0D0).AND.(X.LT.12600.1D0)) THEN
        CURVSURF = TRAMO6(X,TRNL)
        RETURN
    ENDIF
    !
    IF (X.GE.12600.1D0) THEN
        CURVSURF = TRAMO7(TRNL)
        RETURN
    ENDIF
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO0(TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: TRAMO0,TRNL
    !
    TRAMO0 =  2000.0D0 + TRNL
    !
    END FUNCTION
    !
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO1(X,TRNL)
    !
    !.... FIRST SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: X,TRAMO1,TRNL
    !
    TRAMO1 = 2000.0D0 +  1.111111D0*X + TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO2(X,TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: X,TRAMO2,TRNL
    !
    TRAMO2 = -3500.0D0 + 12.7777780D0*X - 0.00617D0*X**2 + TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO3(X,TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: X,TRAMO3,TRNL
    !
    TRAMO3 = 3478.20D0 - 0.619D0*X + 0.0000641D0*X**2 + TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO4(X,TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: X,TRAMO4,TRNL
    !
    TRAMO4 =  4681.3805D0 - 1.238D0*X + 0.000145D0*X**2 + TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO5(TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: TRAMO5,TRNL
    !
    TRAMO5 =  4500.0D0 + TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO6(X,TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: X,TRAMO6,TRNL
    !
    TRAMO6 = 35707.39745D0-5.3615D0*X+0.000207D0*X**2+TRNL
    !
    END FUNCTION
    !
    !*** NEW FUNCTIONS FOR SALT TOP BOUNDARY  ***********************
    !
    FUNCTION TRAMO7(TRNL)
    !
    !.... SECOND SECTION
    !
    IMPLICIT NONE
    !
    REAL(8) :: TRAMO7,TRNL
    !
    TRAMO7 =  1000.0D0 + TRNL
    !
    END FUNCTION
    !
    !=======================================================================
    !
    subroutine readbeta(fname,beta,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
    !
    !     Objetivo: le as permeabilidades de um arquivo
    !
    !----------------------------------------------------------------------
    !
    implicit none
    !
    real(8), dimension(*) :: beta
    integer :: nsd
    !
    integer :: nelemx,nelemy,nelemz,ntype,inbeta,nline,nflag,i,j,k
    real(8) :: xlx,xly,xlz,betaH
    character(len=128) :: NAME, FILE_IN
    character(len=128) :: fname
    character(len=4)   :: EXT,tipo
    character(len=5)   :: C
    integer :: cont, contK, fim, numel_
    !
    !
    inbeta=209
    ! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(C,113)(nr-1)
    C=ADJUSTL(C)
    !
    file_in=fname
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
    WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
    open(inbeta, file= NAME)
    !
    !     dimensoes do dominio
    !
    read(inbeta,*) xlx
    read(inbeta,*) xly
    if(nsd==3) read(inbeta,*) xlz
    !
    !     numero de elementos em cada direcao
    !
    read(inbeta,*) nelemx
    read(inbeta,*) nelemy
    if(nsd==3)read(inbeta,*) nelemz
    !
    !     verificacao do tamanho do vetor para leitura
    !
    if(nsd==2) then
        if(nelemx*nelemy.gt.nelxReserv*nelyReserv*100)then
            write(*,115)(nr-1)
            stop
        end if
    else
        if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
            write(*,115)(nr-1)
            stop
        end if
    endif
    !
    !     ntype = 1; campo exponencial
    !     ntype = 2; campo fractal
    !
    read(inbeta,*) ntype
    !
    !     beta: coeficiente de Hurst
    !
    read(inbeta,*) betaH
    !
    !     leituras vazias
    !
    read(inbeta,*)
    read(inbeta,*)
    !
    !     inicio da leitura do campo
    !
    if(nsd==2) fim=1
    if(nsd==3) fim=nelemz
    contK=0
    !
    do k=1,fim
        !
        contK=1+(k-1)*nelemx*nelemy
        !
        if(nsd==3) then
            read(inbeta,*) nline
            if(nline+1.ne.k) then
                write(*,*) 'Erro na leitura do campo de beta, nline'
                stop
            end if
        endif
        !
        do j=1,nelemy
            !
            read(inbeta,*) nline
            if(nline+1.ne.j) then
                write(*,*) 'Erro na leitura do campo de beta, nline'
                stop
            end if
            !
            cont=contK+(j-1)*nelemx
            read(inbeta,*) (beta(i),i=cont,nelemx+cont-1)
            read(inbeta,*) nflag
            if(nflag.ne.192837465) then
                write(*,*) 'Erro na leitura do campo de beta, nflag'
                stop
            end if
            !
        end do ! nelemy
    end do ! nelemz
    !
    close(inbeta)
    !
    !     calculando a permeabilidade lognormal
    !
    if(nsd==2) numel_=nelemy*nelemx
    if(nsd==3) numel_=nelemz*nelemy*nelemx

    do j=1,numel_
        beta(j)=kgbeta*dexp(rhobeta*beta(j))
    end do
    !
    tipo='K'
    call estatistica(beta,nelemx,nelemy,nelemz,tipo)
    !
111 FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113 FORMAT(I5)
115 FORMAT(                                      &
        '######################################',/, &
        'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
        'LEITURA DO CAMPO DE BETA   ',/, &
        'NUMERO DO CAMPO:',I5,/,                    &
        '######################################',/)
    !
    end subroutine readbeta
    !
    !=======================================================================
    !
    subroutine atribuirBetaBloco(nsd,numel,u,x)
    !
    implicit none
    !
    integer :: nsd,numel
    real*8 :: x(nsd,*)
    real(8), dimension(numel)   :: u
    !
    integer :: i
    real(8) :: xx,yy,zz,xi,xf,yi,yf

    do i=1,numel
        !
        xx=x(1,i)
        yy=x(2,i)
        if(nsd==3)zz=x(3,i)
        !
        !.... Bloco
        !
        u(i)=1.0e-14

        ! ! bloquinho1
        xi=100.00d0-0.001; xf=250d0+0.001
        yi=5.00d0-0.001;   yf=15.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-12
                !         print*, i, u(i)
            end if
        end if

        ! ! bloquinho2
        xi=100.0d0-0.001; xf=250d0+0.001
        yi=30.0d0-0.001;   yf=45.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-10
                !        print*, i, u(i)
            end if
        end if

        ! ! bloquinho3
        xi=475.0d0-0.001; xf=625d0+0.001
        yi=5.0d0-0.001;   yf=15.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-10
                !        print*, i, u(i)
            end if
        end if

        ! ! bloquinho4
        xi=475.0d0-0.001; xf=625d0+0.001
        yi=30.0d0-0.001;  yf=45.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-11
                !        print*, i, u(i)
            end if
        end if!
        !

        ! ! bloquinho5
        xi=850.0d0-0.001; xf=1000.0d0+0.001
        yi=5.0d0-0.001;   yf=15.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-18
                !        print*, i, u(i)
            end if
        end if

        ! ! bloquinho6
        xi=850.0d0-0.001; xf=1000.0d0+0.001
        yi=30.0d0-0.001;  yf=45.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-09
                !        print*, i, u(i)
            end if
        end if
        !

        ! ! bloquinho7
        xi=1225.0d0-0.001; xf=1375.0d0+0.001
        yi=5.0d0-0.001;    yf=15.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-09
                !        print*, i, u(i)
            end if
        end if

        ! ! bloquinho8
        xi=1225.0d0-0.001; xf=1375.0d0+0.001
        yi=30.0d0-0.001;   yf=45.0d0+0.001

        if(xx.gt.xi.and.xx.lt.xf) then
            if(yy.gt.yi.and.yy.lt.yf) then
                u(i)=1.0e-10
                !        print*, i, u(i)
            end if
        end if

        !
    end do
    !
    end subroutine atribuirBetaBloco
    !
    !*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
    !
    FUNCTION BSKEMPTON(YOUNG,NU,ELEMNT)
    !
    !.... COMPUTE SKEMPTON COEFICIENT ON ELEMENT
    !
    IMPLICIT NONE
    !
    INTEGER :: ELEMNT
    REAL(8) :: YOUNG, NU
    REAL*8  :: bulkgrain, bulkfluid, porosity, bulkrock
    REAL*8  :: ALPHARS, ALPHAFS, ALPHARF
    !
    REAL(8) :: BSKEMPTON
    !
    !.... STRESS DIMENSION:  S3DIM = 1.0D0 OR 3.0D0
    !
    bulkgrain = GEOINDIC('BLKGRIN',GEOFORM(ELEMNT))
    bulkfluid = GEOINDIC('BLKFLUD',GEOFORM(ELEMNT))
    porosity  = GEOINDIC('POROSTY',GEOFORM(ELEMNT))
    !!!
    IF (POROSITY.LE.1.0D-5) THEN
        BSKEMPTON = 0.0D0
    ELSE
        bulkrock  = BULK(YOUNG,NU,3.0D0)
        ALPHARS   = 1.0D0 - bulkrock/bulkgrain
        ALPHAFS   = 1.0D0 - 1.0D0/(bulkgrain*bulkfluid)
        ALPHARF   = bulkrock*porosity*bulkfluid
        BSKEMPTON = ALPHARS/(ALPHARS+ALPHARF*ALPHAFS)
        !           write(*,*)ELEMNT,BSKEMPTON
    ENDIF
    !
    END FUNCTION
    !
    !*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
    !
    FUNCTION ALFABIOT(YOUNG,NU,ELEMNT)
    !
    !.... COMPUTE ALFA--BIOT COEFICIENT ON ELEMENT
    !
    IMPLICIT NONE
    !
    INTEGER :: ELEMNT
    REAL(8) :: YOUNG, NU
    REAL(8) :: bulkgrain, BULKROCK, POROSITY
    !
    REAL(8) :: ALFABIOT
    !
    !.... STRESS DIMENSION:  S3DIM = 1.0D0 OR 3.0D0
    !
    POROSITY  = GEOINDIC('POROSTY',GEOFORM(ELEMNT))
    !
    IF (POROSITY.LE.1.0D-5) THEN
        ALFABIOT = 0.0D0
    ELSE
        bulkgrain = GEOINDIC('BLKGRIN',GEOFORM(ELEMNT))
        bulkrock = BULK(YOUNG,NU,3.0D0)
        ALFABIOT = 1.0D0 - bulkrock/bulkgrain
    ENDIF
    !
    END FUNCTION
    !
    !
    !=======================================================================
    !
    SUBROUTINE KOZENYCARMAN(permkx,permky,phi,nelemx,nelemy,nelemz)
    !
    !     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO
    !
    use mMalha,            only: numelReserv
    !
    implicit none
    !
    integer              :: nelemx,nelemy,nelemz
    real(8),dimension(*) :: permkx,permky,phi
    !     real(8)              :: xkkcGeo
    character(len=4)     :: tipo
    !
    integer :: nel
    !
    DO nel=1,numelReserv
        permkx(nel)=xkkcGeo(phi(nel),permkx(nel))
        permky(nel)=permkx(nel)
    END DO
    tipo='K'
    call estatistica(permkx,nelemx,nelemy,nelemz,tipo)
    !
    return
    !
    END SUBROUTINE KOZENYCARMAN
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function getRegionIndex(curElement)
    
    !variables import
    use mMalha, only: nsd
    
    implicit none
    
    !variables input
    integer :: curElement
    
    !variables
    integer :: getRegionIndex
    integer :: lastRegion
    character*12, dimension(9) :: refgeo
    integer :: i
    
    
    data refgeo(1)     ,     refgeo(2),     refgeo(3)/ &
        &     'RESERVATORIO','RIFT_CAPROCK','RIGHT_RESERV'/ &
        &     refgeo(4)     ,     refgeo(5),     refgeo(6)/ &
        &     'SALT_CAPROCK','LEFT__RESERV','POS_SAL_OVER'/ &
        &     refgeo(7)     ,     refgeo(8),     refgeo(9)/ &
        &     'FRONT_RESERV','BASE_CAPROCK','BACK__RESERV'/
    !------------------------------------------------------------------------------------------------------------------------------------
    lastRegion = 3*nsd
    do i=1,lastRegion
        if (refgeo(i).eq.geoform(curElement)) then
            getRegionIndex = i
            exit
        end if
    end do
    
    end function
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function totalCompressibility(curElement, matrixBulk)
    
    implicit none
    
    !variables input
    integer :: curElement
    real*8 :: matrixBulk
    
    ! Variables
    integer :: geoFormIndx
    real*8 :: initialPorosityL
    real*8 :: nInv, mInv, biotCoefficient, grainBulk
    real*8 :: totalCompressibility
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    initialPorosityL = poro0(curElement)
    
    grainBulk = GRAINBLK(geoFormIndx)
    biotCoefficient = 1 - matrixBulk/grainBulk

    nInv = calcNInv(curElement, biotCoefficient)
    mInv = nInv + initialPorosityL/MEANBLKW(geoFormIndx)
    totalCompressibility = mInv + biotCoefficient*biotCoefficient/matrixBulk
    
    end function totalCompressibility
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function rhoCell(curElement)
    ! variable input
    integer :: curElement
    
    ! variables
    real*8 rhoCell
    integer :: geoFormIndx
    real*8 :: poroL, fluidDensity, grainDensity
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    poroL = PORELAGR(geoFormIndx)
    fluidDensity = MEANRHOW(geoFormIndx)
    grainDensity = RHODVECT(geoFormIndx)
    
    rhoCell = poroL * fluidDensity + (1.0d0-poroL) * grainDensity
    end function rhoCell
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcBiotCoefficient(curElement)
    implicit none
    !variables input
    integer :: curElement
    
    !variables
    real*8 :: calcBiotCoefficient
    integer :: geoFormIndx
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    calcBiotCoefficient = 1 - bulk(YUNGVECT(geoFormIndx),POISVECT(geoFormIndx),3.d0)/GRAINBLK(geoFormIndx)
    end function calcBiotCoefficient
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcMatrixBulk(curElement)
    implicit none
    !variables input
    integer :: curElement
    
    !variables
    real*8 :: calcMatrixBulk
    integer :: geoFormIndx
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    calcMatrixBulk = bulk(YUNGVECT(geoFormIndx),POISVECT(geoFormIndx),3.d0)
    
    end function calcMatrixBulk
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcBulkFromMatrix(curTangentMatrix)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: curTangentMatrix(16)
    
    !variables
    real*8 :: calcBulkFromMatrix
    real*8 :: identI(1,4)
    real*8 :: qixi(4,4)
    real*8 :: temp(1,4)
    real*8 :: scalTemp(1,1)
    
    !------------------------------------------------------------------------------------------------------------------------------------
    identI = reshape((/ 1., 1., 0., 1. /), shape(identI))
    
    qixi = transpose(reshape(curTangentMatrix, shape(qixi)))
    
    temp = matmul(identI,qixi)
    scalTemp = matmul(temp,transpose(identI))
    calcBulkFromMatrix = scalTemp(1,1) / 9
    
    end function calcBulkFromMatrix
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcKozenyCarmanPerm(curElement, tensDim)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    integer :: curElement, tensDim
    
    !variables
    real*8 :: calcKozenyCarmanPerm(tensDim)
    real*8 :: deltaN
    integer :: geoFormIndx
    integer :: i
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    deltaN = calcKozenyCarmanCoefficient(poroE(curElement))
    
    do i = 1,tensDim
        calcKozenyCarmanPerm(i) = charLenPorNetworkSquared(i,curElement) * deltaN / waterViscosity(geoFormIndx)
    end do
    
    end function calcKozenyCarmanPerm
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcKozenyCarmanPermInitialPoro(curElement, tensDim)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    integer :: curElement, tensDim
    
    !variables
    real*8 :: calcKozenyCarmanPermInitialPoro(tensDim)
    real*8 :: deltaN
    integer :: geoFormIndx
    integer :: i
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    deltaN = calcKozenyCarmanCoefficient(poro0(curElement))
    
    do i = 1,tensDim
        calcKozenyCarmanPermInitialPoro(i) = charLenPorNetworkSquared(i,curElement) * deltaN / waterViscosity(geoFormIndx)
    end do
    
    end function calcKozenyCarmanPermInitialPoro
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcKozenyCarmanCoefficient(porosity)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: porosity
    
    !variables
    real*8 :: calcKozenyCarmanCoefficient
    
    !------------------------------------------------------------------------------------------------------------------------------------
    calcKozenyCarmanCoefficient = porosity**3/(1-(porosity**2))
    
    end function calcKozenyCarmanCoefficient
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcPoissonUndrained(poisson, biotCoeff, skemptonCoeff)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: poisson, biotCoeff, skemptonCoeff
    
    !variables
    real*8 :: calcPoissonUndrained
    real*8 :: tempVar
    
    !------------------------------------------------------------------------------------------------------------------------------------
    tempVar = biotCoeff * skemptonCoeff * (2*poisson - 1)
    calcPoissonUndrained = (3.0d0*poisson - tempVar)/(tempVar + 3)
    
    end function calcPoissonUndrained
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcBulkUndrained(bulkModulus, biotCoeff, biotMModulus)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: bulkModulus, biotCoeff, biotMModulus
    
    !variables
    real*8 :: calcBulkUndrained
    
    !------------------------------------------------------------------------------------------------------------------------------------
    calcBulkUndrained = bulkModulus + biotCoeff * biotCoeff * biotMModulus
    
    end function calcBulkUndrained
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcSkemptonCoefficient(biotCoeff, biotMModulus, bulkModulusUndrained)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: biotCoeff, biotMModulus, bulkModulusUndrained
    
    !variables
    real*8 :: calcSkemptonCoefficient
    
    !------------------------------------------------------------------------------------------------------------------------------------
    calcSkemptonCoefficient = biotCoeff * biotMModulus / bulkModulusUndrained
    
    end function calcSkemptonCoefficient
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine calcMatrixUndrainedProperties(poissonUndrained, youngUndrained, curElement)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: poissonUndrained, youngUndrained
    integer :: curElement
    
    !variables
    integer :: geoFormIndx
    real*8 :: initialPorosityL
    real*8 :: nInv, mInv, biotMModulus
    real*8 :: bulkMatrix, bulkUndrained
    real*8 :: biotCoeff, skemptonCoeff
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    !calculates the required moduli biotCoeff, biotMModulus and bulkModulus
    biotCoeff = calcBiotCoefficient(curElement)
    
    initialPorosityL = poro0(curElement)
    nInv = (biotCoeff - initialPorosityL)/GRAINBLK(geoFormIndx)
    mInv = nInv + initialPorosityL/MEANBLKW(geoFormIndx)
    biotMModulus = 1/mInv
    
    bulkMatrix = calcMatrixBulk(curElement)
    
    !calculate the undrained bulk modulus
    bulkUndrained = calcBulkUndrained(bulkMatrix, biotCoeff, biotMModulus)
    
    !calculate the skempton coefficient
    skemptonCoeff = calcSkemptonCoefficient(biotCoeff, biotMModulus, bulkUndrained)
    
    !calculate the undrained poisson modulus
    poissonUndrained = calcPoissonUndrained(POISVECT(geoFormIndx), biotCoeff, skemptonCoeff)
    
    !calculate the undrained young modulus
    youngUndrained = 3*bulkUndrained * (1 - 2.0d0*poissonUndrained) 
    
    end subroutine calcMatrixUndrainedProperties
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine calcMatrixProperties(poisson, young, curElement)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: poisson, young
    integer :: curElement
    
    !variables
    integer :: geoFormIndx
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    poisson = POISVECT(geoFormIndx)
    young = YUNGVECT(geoFormIndx)
    
    end subroutine calcMatrixProperties
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    function calcNInv(curElement, biotCoefficient)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    real*8 :: biotCoefficient
    integer :: curElement
    
    !variables
    integer :: geoFormIndx
    real*8 :: calcNInv
    
    !------------------------------------------------------------------------------------------------------------------------------------
    geoFormIndx = getRegionIndex(curElement)
    
    calcNInv = (biotCoefficient - poro0(curElement))/GRAINBLK(geoFormIndx)
    
    end function calcNInv
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    subroutine updatePorosity(curElement, trStrain, trStrainP, pressure, pressureInit)
    !function imports
    
    !variables import
    
    implicit none
    !variables input
    integer :: curElement
    real*8 :: trStrain, trStrainP, pressure, pressureInit
    
    !variables
    real*8 :: nInv
    real*8 :: biotCoefficient
    
    !------------------------------------------------------------------------------------------------------------------------------------
    biotCoefficient = calcBiotCoefficient(curElement)
    nInv = calcNInv(curElement, biotCoefficient)
    
    poroL(curElement) = poro0(curElement) &
                    & + trStrainP & ! plastic porosity
                    & + biotCoefficient*(trStrain - trStrainP) + (pressure-pressureInit) * nInv !elastic component
    poroE(curElement) = poroL(curElement) / (1 + trStrain)
    
    end subroutine updatePorosity
    !************************************************************************************************************************************
    !************************************************************************************************************************************
    end module
