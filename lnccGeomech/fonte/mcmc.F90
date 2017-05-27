!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE leituraEscritaMCMC
  use mHidrodinamicart, only   :  ELEM_CONDP_Hidro=>ELEM_CONDP

  integer :: IFLAG_PRODF,IFLAG_DISF,IFLAG_CONC,IFLAG_PRESF
  REAL(8) :: TPRT_PRODF,DTPRT_PRODF,TPRT_DISF,DTPRT_DISF
  REAL(8) :: TPRT_PRESF,DTPRT_PRESF,TPRT_CONC,DTPRT_CONC
  REAL(8), DIMENSION(2,400) :: PCONDC,PCONDD,PCONDPR
  INTEGER, DIMENSION(400) :: ELEM_CONDC,ELEM_CONDD,ELEM_CONDPR
  CHARACTER(LEN=128) :: DISF_IN,DISF_OUT,PRODF_IN,PRODF_OUT
  CHARACTER(LEN=128) :: PRESF_IN,PRESF_OUT,CONC_OUT,CONC_IN
  INTEGER :: MCMCCONC,MCMCPRODF,MCMCDISF,NCONDC,NCONDD
  INTEGER :: MCMCPRESF,NCONDPR
  
   contains
  
    SUBROUTINE INITMCMC(NCONDP,PCONDP)

    use mGlobaisEscalares, only: tt,tzero,NVEL
!
    IMPLICIT NONE

    INTEGER :: NCONDP
    REAL(8), DIMENSION(2,400) :: PCONDP
!
    CHARACTER(LEN=128) :: FNOME,FLAG,FLAG1
    INTEGER :: ISTAT
    REAL(8) :: TEMPO, TOL=1E-6
!    
    character(len=50) :: keyword_name
    integer :: ierr
!

    keyword_name='mcmc_concentracao'
    call readRandMCMC(keyword_name,IFLAG_CONC, CONC_IN, CONC_OUT, MCMCCONC, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC CONC......................:',ADJUSTL(TRIM(CONC_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCCONC
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(CONC_IN))   
    IF(IFLAG_CONC==1) CALL COND_POINTS(CONC_IN,NCONDC,PCONDC)
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_producao'
    call readRandMCMC(keyword_name,IFLAG_PRODF, PRODF_IN, PRODF_OUT, MCMCPRODF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC PROD......................:',ADJUSTL(TRIM(PRODF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCPRODF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(PRODF_IN))
!     CALL COND_POINTS(PRODF_IN,NCONDP,PCONDP) !está sendo lido no driver
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_deslocamentos'
    call readRandMCMC(keyword_name,IFLAG_DISF, DISF_IN, DISF_OUT, MCMCDISF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC DISLOC....................:',ADJUSTL(TRIM(DISF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCDISF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(DISF_IN))
    IF(IFLAG_DISF==1) CALL COND_POINTS(DISF_IN,NCONDD,PCONDD)
    WRITE(*,*)'#####################'
!
    keyword_name='mcmc_pressao'
    call readRandMCMC(keyword_name,IFLAG_PRESF, PRESF_IN, PRESF_OUT, MCMCPRESF, ierr)
    WRITE(*,*)'#####################'
    WRITE(*,*)'ARQUIVO MCMC PROD.....................:',ADJUSTL(TRIM(PRESF_OUT))
    WRITE(*,*)'NUMERO DE IMPRESSOES:',MCMCPRESF
    WRITE(*,*)'ARQUIVO COM POSICOES DE MONITORAMENTO..:',ADJUSTL(TRIM(PRESF_IN))
    IF(IFLAG_PRESF==1) CALL COND_POINTS(PRESF_IN,NCONDPR,PCONDPR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TEMPO=(tt-tzero)/REAL(NVEL)
    DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
    IF(MCMCCONC.GE.NVEL)THEN
       MCMCCONC = NVEL
       DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
    ELSE
       DO WHILE(DMOD(DTPRT_CONC,TEMPO).GT.TOL)
          MCMCCONC = MCMCCONC+1
          DTPRT_CONC = (tt-tzero)/REAL(MCMCCONC)
       END DO
    END IF
    TPRT_CONC = DTPRT_CONC
    write(*,*)'NUMERO de IMPRESSOES CONC',MCMCCONC
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_CONC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
    IF(MCMCPRODF.GE.NVEL)THEN
       MCMCPRODF= NVEL
       DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
    ELSE
       DO WHILE(DMOD(DTPRT_PRODF,TEMPO).GT.TOL)
          MCMCPRODF= MCMCPRODF+1
          DTPRT_PRODF= (tt-tzero)/REAL(MCMCPRODF)
       END DO
    END IF
    TPRT_PRODF=DTPRT_PRODF  
    write(*,*)'NUMERO de IMPRESSOES PROD',MCMCPRODF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_PRODF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
    IF(MCMCDISF.GE.NVEL)THEN
       MCMCDISF= NVEL
       DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
    ELSE
       DO WHILE(DMOD(DTPRT_DISF,TEMPO).GT.TOL)
          MCMCDISF= MCMCDISF+1
          DTPRT_DISF= (tt-tzero)/REAL(MCMCDISF)
       END DO
    END IF
    TPRT_DISF = DTPRT_DISF
    write(*,*)'NUMERO de IMPRESSOES DISP',MCMCDISF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_DISF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
    IF(MCMCPRESF.GE.NVEL)THEN
       MCMCPRESF= NVEL
       DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
    ELSE
       DO WHILE(DMOD(DTPRT_PRESF,TEMPO).GT.TOL)
          MCMCPRESF= MCMCPRESF+1
          DTPRT_PRESF= (tt-tzero)/REAL(MCMCPRESF)
       END DO
    END IF
    TPRT_PRESF=DTPRT_PRESF  
    write(*,*)'NUMERO de IMPRESSOES PRES',MCMCPRESF
    write(*,*)'Delta t IMPRESSAO        ',DTPRT_PRESF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    CALL MONITORPOSITION(PCONDC, NCONDC, ELEM_CONDC)
    CALL MONITORPOSITION(PCONDP, NCONDP, ELEM_CONDP_Hidro)
    CALL MONITORPOSITION(PCONDD, NCONDD, ELEM_CONDD)
    CALL MONITORPOSITION(PCONDPR,NCONDPR,ELEM_CONDPR)
!
    contains
                             
      subroutine readRandMCMC(keyword_name,IFLAG, ArqIN, ArqOUT, MCMC, ierr)
!
        use mInputReader,    only: file_lines, findKeyword      
        
        implicit none
        
        integer keyword_line, ierr
        character(len=50) keyword_name
        integer :: IFLAG, MCMC
        character(len=128) :: ArqIN, ArqOUT

        ierr=0
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            ierr=1
            IFLAG=0
            return
        end if
        read(file_lines(keyword_line:),'(i5)') IFLAG
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),'(a)') ArqOut
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),*) MCMC
        keyword_line = keyword_line + 1
        if(keyword_name.ne.'mcmc_producao') &
           read(file_lines(keyword_line:),'(a)') ArqIn

    end subroutine 
    
  END SUBROUTINE INITMCMC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    SUBROUTINE COND_POINTS(file_in,N,PCOND)
!
      implicit none
!
      integer, intent(out)           :: N
      integer                        :: idata,i,j,ISTAT
      character(len=128), intent(in) :: FILE_IN
      real(8)                        :: TOL
      REAL(8), intent(out), DIMENSION(2,100) :: PCOND
!
      TOL=1.e-6
!
      idata=168
      OPEN(UNIT=IDATA,FILE=FILE_IN,STATUS='UNKNOWN',&
           ACTION='READ',IOSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_IN
         STOP
      ELSE
         WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FILE_IN))
      END IF
!
      read(idata,*)N
      write(*,101)N
!
      do i=1,N
      read(idata,*)(pcond(j,i),j=1,2)
      write(*,201)i,pcond(1,i),pcond(2,i)
      enddo
!
      do i=1,N
         do j=1,2
            pcond(j,i)=pcond(j,i)+TOL
         enddo
      enddo
      CLOSE(IDATA)
!
100   format(I5)
101   format('Numero de Pontos de Controle:',i5,/)
200   format(3F12.7)
201   format('Coordenada do Ponto',i5,' (',F10.4,',',F10.4,')')
!
    end SUBROUTINE COND_POINTS
!
!=======================================================================
!
    SUBROUTINE MONITORMASSATOTAL(U,POROS,T,iifmass_out)
      use mMalha,            only: numelReserv,numel
      use mPropGeoFisica,    only: hy,hx
      use mPropGeoFisica,    only: rhow
      use mPropGeoFisica,    only: GEOFORM
!
      IMPLICIT NONE
!
      INTEGER :: N,INFILE,NEL
      REAL(8) :: T,VOL
      REAL(8), DIMENSION(numelReserv) :: U,POROS
      INTEGER :: ISTAT,I,ZERO=0D0
      CHARACTER(len=128) :: NAME,FILE_IN
      CHARACTER(len=4)   :: EXT
      CHARACTER(len=5)   :: C
      CHARACTER(len=128) :: iifmass_out
!
      INFILE = 170
      WRITE(C,113)ZERO
      C=ADJUSTL(C)
      FILE_IN = IIFMASS_OUT
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM('_')//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
!
      OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
           ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',IIFMASS_OUT
         STOP
      END IF
!
      VOL=0.0D0
      DO NEL=1,numel
         IF(GEOFORM(NEL).EQ.'RESERVATORIO')THEN
            VOL = VOL+U(NEL)*POROS(NEL)
         END IF
      END DO
      write(*,*)'##### MONITORANDO MASSA TOTAL DE AGUA #####'
      WRITE(INFILE,400)T,VOL*RHOW*HX*HY
!
      CLOSE(INFILE)
!
113   FORMAT(I5)
400   FORMAT(100(e15.8,2x))
    END SUBROUTINE MONITORMASSATOTAL
    
    SUBROUTINE MONITORPOSITION(PCOND,N,ELEM_COND)
      use mMalha,         only: x,NUMEL,conecNodaisElem
      IMPLICIT NONE
!
      INTEGER :: N,J,I,NO,NEL
      REAL(8), DIMENSION(2,*) :: PCOND
      REAL(8) :: X1,X2,Y1,Y2
      INTEGER, DIMENSION(100)   :: ELEM_COND
!
      IF(N.GT.400)THEN
         WRITE(*,*)'NUMERO EXCESSIVO DE PONTOS'
         WRITE(*,*)'MODIFICAR EM MODULE MMCMC '
         STOP
      END IF
!
      WRITE(*,*)'########################################'
      WRITE(*,*)'### MONITORAMENTO ######################'
      WRITE(*,*)'### NUMERO DE PONTOS:',N
      WRITE(*,*)'########################################'
      DO 200 I=1,N
         DO NEL=1,NUMEL
            X1=x(1,conecNodaisElem(1,NEL))
            X2=x(1,conecNodaisElem(2,NEL))
            Y1=x(2,conecNodaisElem(1,NEL))
            Y2=x(2,conecNodaisElem(4,NEL))
            IF(PCOND(1,I).LT.X2)THEN
               IF(PCOND(1,I).GE.X1)THEN
                  IF(PCOND(2,I).LT.Y2)THEN
                     IF(PCOND(2,I).GE.Y1)THEN
                        ELEM_COND(I)=NEL
                        GO TO 200
                     END IF
                  END IF
               END IF
            END IF
         END DO
200      CONTINUE
!
       END SUBROUTINE MONITORPOSITION
!
!=======================================================================
!
  SUBROUTINE MONITOR(FNOME,PCOND,N,ELEM_COND,U,T)
    IMPLICIT NONE
!
    INTEGER :: N,INFILE
    REAL(8) :: T
    REAL(8), DIMENSION(*) :: U
    REAL(8), DIMENSION(2,100) :: PCOND
    CHARACTER(LEN=128) :: FNOME
    INTEGER :: ISTAT,I,ZERO=0D0
    INTEGER, DIMENSION(100)   :: ELEM_COND
    CHARACTER(len=128) :: NAME,FILE_IN
    CHARACTER(len=4)   :: EXT
    CHARACTER(len=5)   :: C
!
    INFILE = 163
    WRITE(C,113)ZERO
    C=ADJUSTL(C)
    FILE_IN = FNOME
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
!
    OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
         ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
       STOP
    END IF
!
    write(*,*)'##### MONITORANDO #####'
    WRITE(INFILE,400)T,(U(ELEM_COND(I)),I=1,N)
!
    CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(100(e15.8,2x))
  END SUBROUTINE MONITOR
!
!=======================================================================
!
  SUBROUTINE MONITORV(FNOME,PCOND,N,ELEM_COND,U,T)
    use mGeomecanica, only : ndofD
    USE mMalha,       only: x,numnp,conecNodaisElem
    IMPLICIT NONE
!
    INTEGER :: N,INFILE,NEL,J
    REAL(8) :: T
!    REAL(8), DIMENSION(ndofD,numnp) :: U
    REAL(8), DIMENSION(ndofD,*) :: U
    REAL(8), DIMENSION(2,*) :: PCOND
    CHARACTER(LEN=128) :: FNOME
    INTEGER :: ISTAT,I,ZERO=0D0
    INTEGER, DIMENSION(100)   :: ELEM_COND
    CHARACTER(len=128) :: NAME,FILE_IN
    CHARACTER(len=4)   :: EXT
    CHARACTER(len=5)   :: C
    REAL(8) :: X1,X2,Y1,Y2
!
    INFILE = 164
    WRITE(C,113)ZERO
    C=ADJUSTL(C)
    FILE_IN = FNOME
    EXT='.dat'
    NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
    NAME=ADJUSTL(TRIM(NAME))
!
    OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
         ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
       STOP
    END IF
!
    write(*,*)'##### MONITORANDO VET. ',ADJUSTL(TRIM(NAME)),' #####',N
    WRITE(INFILE,400)T,(U(2,conecNodaisElem(3,ELEM_COND(NEL))),NEL=1,N)
!
    CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(100(e15.8,2x))
  END SUBROUTINE MONITORV
!
!=======================================================================
!
      SUBROUTINE MONITORPROD(FNOME,PCOND,N,ELEM_COND,VEL,U,T)
      use mMalha,            only: numLadosReserv,numel
      use mMalha,            only: conecLadaisElem
      use mHidrodinamicart,  only: ndofV
      use mPropGeoFisica,    only: hy,xlw,xlo,XLT
!
      IMPLICIT NONE
!
      INTEGER :: N,INFILE,NEL,LADO
      REAL(8) :: T,PROD,OIL_MOBI,ST
      REAL(8), DIMENSION(numel) :: U
      REAL(8), DIMENSION(ndofV,numLadosReserv) :: VEL
      REAL(8), DIMENSION(2,100) :: PCOND
      CHARACTER(LEN=128) :: FNOME
      INTEGER :: ISTAT,I,ZERO=0D0
      INTEGER, DIMENSION(100)   :: ELEM_COND
      CHARACTER(len=128) :: NAME,FILE_IN
      CHARACTER(len=4)   :: EXT
      CHARACTER(len=5)   :: C
!
      INFILE = 169
      WRITE(C,113)ZERO
      C=ADJUSTL(C)
      FILE_IN = FNOME
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
!
      OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
           ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
         STOP
      END IF
!
      PROD=0.0D0
      DO NEL=1,N
         ST   = U(ELEM_COND(NEL))
         LADO = conecLadaisElem(2,ELEM_COND(NEL))
         OIL_MOBI = XLO(ST)/XLT(ST)
         PROD = PROD + VEL(1,LADO)*OIL_MOBI
      END DO
      write(*,*)'##### MONITORANDO PRODUCAO #####'
      WRITE(INFILE,400)T,PROD*HY!*dtprt_prodF
!
      CLOSE(INFILE)
!
113   FORMAT(I5)
400   FORMAT(100(e15.8,2x))
    END SUBROUTINE MONITORPROD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     subroutine imprimirCondicoesIniciaisMCMC(phi, satElem)
      use mMalha,            only: numelReserv
      use mLeituraEscritaSimHidroGeoMec,   only: iflag_mass, ifmass_out
!
      implicit none
!
      real*8, intent(in) :: phi(numelReserv), satElem(numelReserv)
      real*8  :: ZERO=0.D0
!
!.... imprime a condicao inicial da massa de agua
!
       if(iflag_mass==1)then
          call MONITORMASSATOTAL(satElem,phi,zero,ifmass_out)
       end if
      
    END SUBROUTINE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
      subroutine imprimirSolucaoNoTempoMCMC(sat,DIS,PORE,pressaoElem,velocLadal,NCONDP,PCONDP,tempo)
       use mHidrodinamicart,  only: ndofP, ndofV
       use mGeomecanica,      only: ndofD
       use mMalha,            only: nsd, numel, numelReserv, numLadosReserv, numLadosElem, conecLadaisElem
       use mMalha,            only: nen, numnp
       use mLeituraEscritaSimHidroGeoMec,   only: IFLAG_MASS,DTPRT_MASS, ifmass_out, TPRT_MASS
! 
      implicit none
!
      real*8, intent(in)  :: pressaoElem(ndofP,numelReserv), velocLadal(ndofV,numLadosReserv)
      real*8, intent(in)  :: sat(numelReserv),DIS(ndofD,numnp),PORE(numelReserv)
      integer, intent(in) :: NCONDP
      REAL(8), intent(in) :: PCONDP(2,400), tempo
!
      real*8        :: TOL=1e-6
! 
!
!.... imprime a condicao inicial da massa de agua
!
       IF(IFLAG_MASS==1)THEN
          IF(ABS(tempo-TPRT_MASS).LE.TOL)THEN
             TPRT_MASS = TPRT_MASS+DTPRT_MASS
             call MONITORMASSATOTAL(sat,PORE,tempo,ifmass_out)
          END IF
       end if
!
      IF(IFLAG_CONC==1)THEN
         IF(ABS(tempo-TPRT_CONC).LE.TOL)THEN
	   TPRT_CONC = TPRT_CONC+DTPRT_CONC
	   CALL MONITOR(CONC_OUT,PCONDC,NCONDC,ELEM_CONDC,sat,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_PRODF==1)THEN
         IF(ABS(tempo-TPRT_PRODF).LE.TOL)THEN
	   TPRT_PRODF = TPRT_PRODF+DTPRT_PRODF
	   CALL MONITORPROD(PRODF_OUT,PCONDP,NCONDP,ELEM_CONDP_Hidro,velocLadal,sat,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_PRESF==1)THEN
         IF(ABS(tempo-TPRT_PRESF).LE.TOL)THEN
	   TPRT_PRESF = TPRT_PRESF+DTPRT_PRESF
	   CALL MONITOR(PRESF_OUT,PCONDPR,NCONDPR,ELEM_CONDPR,pressaoElem,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_DISF==1)THEN
         IF(ABS(tempo-TPRT_DISF).LE.TOL)THEN
	   TPRT_DISF = TPRT_DISF+DTPRT_DISF
	   CALL MONITORV(DISF_OUT,PCONDD,NCONDD,ELEM_CONDD,DIS,tempo)
         ENDIF
      ENDIF

   END SUBROUTINE
    
END MODULE leituraEscritaMCMC
