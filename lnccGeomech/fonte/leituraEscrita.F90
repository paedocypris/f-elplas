!==============================================================================
!         programa de elementos finitos  
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
! 
!         + implementacoes de Abimael Loula
!
!         + novo projeto: modular e fortran 90, por
!         Eduardo Garcia,        bidu@lncc.br 
!         Tuane Lopes,           tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!==============================================================================
      module mLeituraEscrita
!
      implicit none
!
      integer :: iin,iecho,icoords,iconects,iconectsL
      integer :: ignuplot
      integer :: nprint
      integer :: iflag_tipoPrint
!
      contains
!
!==========================================================================
!
      subroutine abrirArquivosInformacoesMalha()
!  
!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
      character(len=20) :: nomeIecho
!
      iecho         = 12
      icoords       = 13
      iconects      = 14
      iconectsL     = 15
!
      nomeIecho='echo.dat'
      open(unit=iecho , file= nomeIecho)
!

   end subroutine 
!
!========================================================================
!
      subroutine fecharArquivosBase()
!
      close(iecho )
      close(icoords )
      close(ignuplot)
      close(iconectsL)
      close(iconects )
       
      end subroutine fecharArquivosBase
!      
! =======================================================================
!   
      subroutine echo
!
      implicit none
!
!.... program to echo input data
!
      character*4 ia(20)
      integer :: iech, i
!
!     cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return
!
      write(iecho,2000) iech
      backspace iin
!
      do 100 i=1,100000
         read(iin,3000,end=200) ia
         if (mod(i,50).eq.1) write(iecho,4000)
         write(iecho,5000) ia
  100 continue
!
  200 continue
      rewind iin
      read(iin,1000) iech
!
      return
!
 500  format('programa de elementos finitos em fortran 90 baseado em:',// &
      'The Finite Element Method, Hughes, T. J. R., (2003)'//)
 1000 format(16i10)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
     ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
     '    eq. 0, no echo of input data                        ',   /5x,&
     '    eq. 1, echo input data                              ',   ///)
 3000 format(20a4)
 4000 format(' ',8('123456789*'),//)
 5000 format(' ',20a4)

      end subroutine echo
!
!**** new *******************************************************************
!
      subroutine printp(a,idiag,neq,nsq,*)
      use mGlobaisEscalares
!
!.... program to print array d after Crout factorization 
!        a = u(transpose) * d * u
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: neq, nsq, iecho
      real*8 :: a(*)
      integer :: idiag(*)
!
      integer :: n, i
!
      do 100 n=1,neq
         if (mod(n,50).eq.1) write(iecho,1000) nsq
         write(iecho,1000)
         i = idiag(n)
         write(iecho,2000) n,a(i)
  100 continue
!
      return 1
!
 1000 format('1',' array d of factorization',/&
     ' a = u(transpose) * d * u ',                                //5x,&
     ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end subroutine
!
!**** new *******************************************************************
!
      subroutine prntel(mat,conectElem,nen,numel,tipo)
      implicit none
!
!.... program to print data for element with "nen" nodes
!
!        note: presently the label formats are limited to
!              elements with one to nine nodes
!
      integer :: nen, numel
      integer :: mat(*),conectElem(nen,*)
      integer :: tipo
!
      integer n, i
!
      if(tipo==1) then
      write(iconects,*) "# Conectividades nodais"
      do n=1,numel
        write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if

      if(tipo==2) then
      write(iconectsL,*) "# Conectividades ladais"
      do  n=1,numel
        write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if
!
      return
!
 2000 format(1x,i10,9(2x,i10))
 3000 format(1x,i10,7(2x,i10))
      end subroutine
!
!**** new *******************************************************************
!
      subroutine printResultado(dva, ndof, numnp, inicio, fim, icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof, numnp, inicio, fim, icode
      real*8  :: dva(ndof,numnp)
!
      integer :: n, i
!
      write(icode,*) "# Solucao"
      do 100 n=inicio,fim
         write(icode,2000) n,(dva(i,n),i=1,ndof)
         !write(*,*) n,(dva(i,n),i=1,ndof)
  100 continue
!
      return
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine
!
!**** new *******************************************************************
!
      subroutine prtgnup(name,x,dva,nsd,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: nsd, ndof, numnp, icode
      character*4 name(11)
      real*8 :: x(nsd,*),dva(ndof,*)
!
      integer :: n, j, i
!
      write(icode,*) name
      do 100 n=1,numnp
         write(icode,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
  100 continue
!
      return
!
 2000 format(6(1pe13.6,2x))
      end subroutine
      

!
!=======================================================================
!     
      subroutine prt(nsd,numel,t0,u,iunit)
!      
      use mMalha, only: xc
!      
      implicit none
!     
!     imprime campos escalares para o gnuplot ou para o matlab
!
      integer                   :: nsd,numel
      real(8), dimension(*)     :: u
      real(8)                   :: t0
      integer :: iunit
!     
      integer :: nel
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
!     
      do nel=1,numel
         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!     
!=======================================================================
!
      subroutine prtvB(nsd,numel,t0,velocLadal,ndofV,  conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: xc
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV ,numLadosElem
      real(8), dimension(ndofV,*) :: velocLadal
      real(8)                   :: t0
      integer                   :: conecLadaisElem(numLadosElem,numel)
!
      integer :: nel
      real(8) :: vc(nsd)
      real*8 :: mediaCentro
!
      integer :: iunit

      vc=0.0
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        vc(1) = (velocLadal(1,conecLadaisElem(2,nel))+velocLadal(1,conecLadaisElem(4,nel)))/2.0
        vc(2) = (velocLadal(1,conecLadaisElem(1,nel))+velocLadal(1,conecLadaisElem(3,nel)))/2.0
        if(nsd==3) vc(3) = (velocLadal(1,conecLadaisElem(5,nel))+velocLadal(1,conecLadaisElem(6,nel)))/2.0
        mediaCentro=sum(vc)/nsd
        write(iunit,"(6(f25.15,2x))")xc(1:nsd,nel),mediaCentro
!
      end do
!
      write(iunit,*)
!
      end subroutine
!     
!=======================================================================
!     
      subroutine prtv(nsd,numel,ndofV,numLadosReserv,t0,u,iunit)
!      
      use mMalha, only: xc
!
      implicit none
!     
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV, numLadosReserv
      real(8), dimension(ndofV,numLadosReserv) :: u
      real(8)                   :: t0
!
      integer :: nel
!     
      integer :: iunit
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!     
      do nel=1,numel
         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(1,nel),u(2,nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!
!=======================================================================
!
      subroutine prtvB3D(x,nen,nsd,numel,conecNodaisElem,t0,velocLadal,ndofV, numLados, conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: local 
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      real*8 :: x(nsd,*)
      integer                   :: nen,numel,nsd, ndofV ,numLados,numLadosElem
      real(8), dimension(ndofV,numLados) :: velocLadal
      real(8)                   :: t0
      integer                   :: conecNodaisElem(nen,numel),conecLadaisElem(numLadosElem,numel)
!
      integer  :: nel
      real*8   :: xl(nsd,nen)
      real(8)  :: xg,yg,zg,vcx,vcy,vcz
      integer :: lado1,lado2,lado3,lado4,lado5,lado6
!
      integer :: iunit
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
        xg = sum(xl(1,1:nen))/nen
        yg = sum(xl(2,1:nen))/nen
        zg = sum(xl(3,1:nen))/nen

        lado1 = conecLadaisElem(1,nel);  lado2 = conecLadaisElem(2,nel);
        lado3 = conecLadaisElem(3,nel);  lado4 = conecLadaisElem(4,nel);
        lado5 = conecLadaisElem(5,nel);  lado6 = conecLadaisElem(6,nel);
        vcx = (velocLadal(1,lado2)+velocLadal(1,lado4))/2.0
        vcy = (velocLadal(1,lado1)+velocLadal(1,lado3))/2.0
        vcz = (velocLadal(1,lado5)+velocLadal(1,lado6))/2.0
        write(iunit,"(6(f25.15,2x))") xg,yg,zg,vcx,vcy,vcz

      end do
!
      write(iunit,*)
!
      end subroutine
!
!
     END MODULE
!
!******** ************ ************ ************ *********** **********
!
