    !=============================================================================
    !         programa de elementos finitos
    !         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
    !
    !         + implementacoes de Abimael Loula
    !
    !         + novo projeto: modular e fortran 90, por
    !         Eduardo Garcia,        bidu@lncc.br
    !         Tuane Lopes,           tuane@lncc.br
    !
    !=============================================================================
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
    !
    !**** NEW **** MODIFIED FOR IRREGULAR MESH  *****************
    !
    subroutine printf(f,ndof,numnp,nlv,iecho)
    !
    !.... program to print prescribed force and boundary condition data
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer ndof, numnp, nlv, iecho
    real*8 :: f(ndof,numnp,*)
    !
    logical lzero
    integer :: nn, n, i
    !
    nn = 0
    !
    do 100 n=1,numnp
        call ztest(f(1,n,nlv),ndof,lzero)
        if (.not.lzero) then
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
            write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
        endif
100 continue
    !
    return
    !
1000 format('1',&
        ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
        '  b o u n d a r y   c o n d i t i o n s'//5x,&
        ' load vector number = ',i10///5x,&
        ' node no.',6(13x,'dof',i1,:)/)
2000 format(6x,i10,10x,6(1pe15.8,2x))
    end subroutine
    !
    !**** new ********************************************************************
    !
    subroutine printd(name,dva,ndof,numnp,icode)
    !
    !.... program to print kinematic data
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer :: ndof,numnp, icode
    character (LEN=*) ::  name
    real*8 dva(ndof,*)
    !
    logical lzero
    integer nn, n, i
    !
    nn = 0
    !
    do 100 n=1,numnp
        !       call ztest(dva(1,n),ndof,lzero)
        !       if (.not.lzero) then
        !          nn = nn + 1
        !          if (mod(nn,50).eq.1) &
        !            write(icode,1000) name,(i,i=1,ndof)
        write(icode,2000) n,(dva(i,n),i=1,ndof)
        !       endif
100 continue
    !
    return
    !
1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
2000 format(1x,i10,2x,6(1pe30.10,2x))
    end subroutine

    !
    !=============================================================================
    !
    subroutine ztest(a,n,lzero)
    !
    !.... program to determine if an array contains only zero entries
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer, intent(in)    :: n
    real*8, intent(in)     :: a(n)
    logical, intent(inout) :: lzero
    !
    integer :: i
    !
    lzero = .true.
    !
    do 100 i=1,n
        if (a(i).ne.0.0d0) then
            lzero = .false.
            return
        endif
100 continue
    !
    end
    !
    !**** new *******************************************************************
    !
    subroutine timing(tempoDeParede)
    implicit none
    !
    !.... program to determine elapsed cpu time
    !
    !
    real*8, intent(inout) :: tempoDeParede
    !
    character(LEN=8)  :: date
    character(LEN=10) :: minhaHora
    character(LEN=5)  :: zone
    integer,dimension(8) :: values
    integer ::  horas, minutos, segundos, milesimosSeg
    !
    call date_and_time(date,minhaHora,zone,values);
    !
    horas=values(5);      minutos=values(6);
    segundos=values(7); milesimosSeg=values(8);
    !
    tempoDeParede = (60*horas+minutos)*60.0+segundos+milesimosSeg/1000.00
    !
    return
    end

    !
    !=======================================================================
    !
    subroutine clear(a,m)
    !
    !.... program to clear a floating-point array
    !
    implicit real*8 (a-h,o-z)
    !
    !.... remove above card for single-precision operation
    !
    INTEGER :: M
    REAL*8 :: a(*)
    integer :: i
    !
    do 100 i=1,m
        a(i) = 0.0d0
100 continue
    !
    return
    end subroutine
    !
    !======================================================================
    !
    subroutine dividirTrabalho   (posI_, posF_, nth_, id_, inicio_, fim_)
    implicit none
    !
    integer, intent(in)   :: posI_, posF_, nth_, id_
    integer, intent(out)  :: inicio_, fim_
    !
    integer :: tamanhoBloco, resto, elementos
    LOGICAL :: distribuirResto = .true.

    tamanhoBloco = int((posF_-posI_+1)/nth_)
    resto = mod((posF_-posI_+1),nth_)
    inicio_ = posI_ + id_*tamanhoBloco
    fim_    = inicio_ + tamanhoBloco - 1

    if(resto > 0) then
        if (distribuirResto) then
            if(id_<resto) then
                if(id_ /= 0) inicio_ = inicio_+ id_
                fim_    = fim_ + id_ + 1
            else
                inicio_= inicio_ + resto
                fim_   = fim_    + resto
            end if
        else
            if(id_ == nth_ - 1) fim_ = fim_ + resto
        end if
    end if
    elementos = fim_ - inicio_ + 1
    ! write(*,'(a,4(a,i0))')rotulo,', id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos
    !  write(*,'(4(a,i0))')' id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos

    end subroutine dividirTrabalho
    !
    !**** new **********************************************************************
    !
    subroutine gerarLabel(label,tempo)
    !
    implicit none

    character(LEN=21), intent(out) :: label
    real*8, intent(in) :: tempo
    !
    character(LEN=21) :: labelAux, num
    integer :: i

    write(num,'(f20.4)') tempo
    !     labelAux="t="//ADJUSTL(num)
    labelAux="t="//num
    do i = 1, 21
        if(labelAux(i:i) .ne. ' ') then
            label(i:i) = labelAux(i:i)
        else
            label(i:i) = '0'
        end if
    end do

    end subroutine

    !
    !**** new **********************************************************************
    !
    subroutine diag(idiag,neq,n)
    !
    implicit none
    !
    !.... program to compute diagonal addresses of left-hand-side matrix
    !
    integer :: neq, n
    integer :: idiag(neq)
    !
    integer :: i

    n = 1
    idiag(1) = 1
    if (neq.eq.1) return
    !
    do 100 i=2,neq
        idiag(i) = idiag(i) + idiag(i-1) + 1
100 continue
    n = idiag(neq)
    !
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine load(id,f,brhs,ndof,numnp,nlvect)
    !
    !.... program to accumulate nodal forces and transfer into
    !        right-hand-side vector
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: f(ndof,numnp,*),brhs(*)
    integer :: ndof, numnp, nlvect
    !
    integer :: nlv
    integer :: i, j, k
    !
    do 300 i=1,ndof
        !
        do 200 j=1,numnp
            k = id(i,j)
            if (k.gt.0) then
                !
                do 100 nlv=1,nlvect
                    brhs(k) = brhs(k) + f(i,j,nlv)
100             continue
                !
            endif
            !
200     continue
        !
300 continue
    !
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine dirichletConditions(id,d,f,ndof,numnp,nlvect)
    !
    !.... program to compute displacement boundary conditions
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer*4:: ndof, numnp, nlvect
    integer*4:: id(ndof,*)
    real*8  :: d(ndof,*),f(ndof,numnp,*)
    !
    integer*4:: i, j, k, lv
    real*8  :: val
    !
    !
    do 300 i=1,ndof
        !
        do 200 j=1,numnp
            !
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.d0
            do 100 lv=1,nlvect
                val = val + f(i,j,lv)
100         continue
            !
            d(i,j) = val
            !
200     continue
        !
300 continue

    end subroutine


    SUBROUTINE LOADTIME(id,f,brhs,ndof,numnp,nlvect,XTIME)
    !
    !.... program to accumulate nodal forces and transfer into
    !        right-hand-side vector
    !
    IMPLICIT NONE
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: f(ndof,numnp,*),brhs(*)
    integer :: ndof, numnp, nlvect
    !
    integer :: nlv
    integer :: i, j, k
    REAL*8  :: XTIME
    !
    DO 300 I=1,NDOF
        DO 200 J=1,NUMNP
            K = ID(I,J)
            IF (K.GT.0) THEN
                DO 100 NLV=1,NLVECT
                    BRHS(K) = BRHS(K) + F(I,J,NLV)*XTIME
100             CONTINUE
            ENDIF
200     CONTINUE
300 CONTINUE
    !
    RETURN
    END SUBROUTINE
    !
    !**** new **********************************************************************
    subroutine btod(id,d,brhs,ndof,numnp)
    !
    !.... program to perform transfer from r.h.s. to displacement array
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: i, j, k
    integer :: ndof, numnp
    integer, dimension(ndof,numnp) :: id
    real*8,  dimension(ndof,numnp) :: d
    real*8,  dimension(*)          :: brhs
    !
    do 200 i=1,ndof
        !
        do 100 j=1,numnp
            k = id(i,j)
            if (k.gt.0) then
                d(i,j) = brhs(k)
            end if
100     continue
        !
200 continue
    !
    return
    end subroutine


    !**** new **********************************************************************
    subroutine kdbc(eleffm,elresf,dl,nee)
    !
    !.... program to adjust load vector for prescribed displacement
    !     boundary condition
    use mGlobaisEscalares, only: zero
    !
    implicit real*8 (a-h,o-z)
    !
    !.... remove above card for single-precision operation
    !
    integer :: nee
    real*8  :: eleffm(nee,*),elresf(*),dl(*)
    !
    integer :: i,j
    real*8  :: val
    !
    do 200 j=1,nee
        !
        val=dl(j)
        if(val.eq.zero) go to 200
        !
        do 100 i=1,nee
            elresf(i)=elresf(i)-eleffm(i,j)*val
100     continue
        !
200 continue
    !
    return
    end subroutine

    subroutine kdbc2(eleffm,elresf,dl,nee,lmt,nel)
    !
    !.... program to adjust load vector for prescribed displacement
    !     boundary condition
    use mglobaisescalares, only: zero
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !

    real*8 :: eleffm(nee,*),elresf(*),dl(*),val
    integer :: nel, nee
    integer*4, pointer :: lmt(:,:,:)
    !
    integer :: i, j, l
    integer*4:: lm(nee)

    lm(:) =reshape(lmt(:,:,nel),(/nee/));
    !
    do 200 j=1,nee
        l=lm(j)
        !
        val=dl(j)
        !
        !
        if(l.gt.0) go to 200
        if(val.eq.zero) go to 200
        !
        !
        do 100 i=1,nee
            elresf(i)=elresf(i)-eleffm(i,j)*val
100     continue
        !
200 continue
    !
    return
    end subroutine

    !**** new **********************************************************************
    subroutine pivots(a,idiag,neq,nsq,iecho,*)
    !
    !.... program to determine the number of zero and negative terms in
    !        array d of factorization a = u(transpose) * d * u
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(*)
    integer :: idiag(*)
    integer :: neq, nsq, iecho
    !
    integer :: iz, in, n, i
    !
    iz = 0
    in = 0
    !
    do 100 n=1,neq
        i = idiag(n)
        if (a(i).eq.0.) iz = iz + 1
        if (a(i).lt.0.) in = in + 1
100 continue
    !
    write(iecho,1000) nsq,iz,in
    !
    return 1
    !
1000 format(' ',&
        ' zero and/or negative pivots encountered                ', ///5x,&
        ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
        ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
        ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)
    !
    end subroutine

    !**** new **********************************************************************
    subroutine ftod(id,d,f,ndof,numnp,nlvect)
    !
    !.... program to compute displacement boundary conditions
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: d(ndof,*),f(ndof,numnp,*)
    integer :: ndof, numnp, nlvect
    !
    integer :: i, j, k, lv
    real*8  :: val
    !
    !
    do 300 i=1,ndof
        !
        do 200 j=1,numnp
            !
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
            do 100 lv=1,nlvect
                val = val + f(i,j,lv)
100         continue
            !
            d(i,j) = val
            !
200     continue
        !
300 continue
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine ftodTIME(id,d,f,ndof,numnp,nlvect,XTIME)
    !
    !.... program to compute displacement boundary conditions
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: d(ndof,*),f(ndof,numnp,*)
    integer :: ndof, numnp, nlvect
    !
    integer :: i, j, k, lv
    real*8  :: val, XTIME
    !
    do 300 i=1,ndof
        do 200 j=1,numnp
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
            do 100 lv=1,nlvect
                val = val + f(i,j,lv)*XTIME
100         continue
            d(i,j) = val
200     continue
300 continue
    !
    return
    !
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine colht(idiag,lm,ned,nen,numel,neq)
    !
    !.... program to compute column heights in global left-hand-side matrix
    !
    implicit none
    integer :: idiag(*),lm(ned,nen,*)
    integer :: ned, nen, numel, neq
    !
    integer :: i, j, k
    integer :: m, min, num
    !
    do 500 k=1,numel
        min = neq
        !
        do 200 j=1,nen
            !
            do 100 i=1,ned
                num = lm(i,j,k)
                if (num.gt.0) min = min0(min,num)
100         continue
            !
200     continue
        !
        do 400 j=1,nen
            do 300 i=1,ned
                num = lm(i,j,k)
                if (num.gt.0) then
                    m = num - min
                    if (m.gt.idiag(num)) then
                        idiag(num) = m
                    endif
                endif
                !
300         continue

400     continue
        !
500 continue
    !
    return
    end subroutine

    !**** new **********************************************************************
    function coldot(a,b,n)
    !
    !.... program to compute the dot product of vectors stored column-wise
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(*),b(*)
    integer :: n
    !
    real*8  :: coldot
    integer :: i
    !

    coldot = 0.0d0
    !
    !       do 100 i=1,n
    !       coldot = coldot + a(i)*b(i)
    coldot=dot_product(a(1:n),b(1:n))

    !   100 continue

    !
    return
    end function
    !
    !**** new **********************************************************************
    !
    subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
    !
    !.... program to add rectangular matrices
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(ma,*),b(mb,*),c(mc,*)
    integer :: ma,mb,mc,m,n,iopt
    !
    integer :: i,j
    !
    go to (1000,2000,3000),iopt
    !
    !.... iopt = 1, add entire matrices
    !
1000 do 1200 j=1,n
        !
        do 1100 i=1,m
            c(i,j) = a(i,j) + b(i,j)
1100    continue
        !
1200 continue
    return
    !
    !.... iopt = 2, add lower triangular and diagonal elements
    !
2000 do 2200 j=1,n
        !
        do 2100 i=j,m
            c(i,j) = a(i,j) + b(i,j)
2100    continue
        !
2200 continue
    return
    !
    !.... iopt = 3, add upper triangular and diagonal elements
    !
3000 do 3200 j=1,n
        !
        do 3100 i=1,j
            c(i,j) = a(i,j) + b(i,j)
3100    continue
        !
3200 continue
    return
    !
    end subroutine

    !**** new **********************************************************************
    function rowdot(a,b,ma,mb,n)
    !
    !.... program to compute the dot product of vectors stored row-wise
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    real*8  :: a(ma,*),b(mb,*)
    integer :: ma, mb, n
    !
    real*8  :: rowdot
    integer :: i
    !
    rowdot = 0.0d00
    !
    do i=1,n
        rowdot = rowdot + a(1,i)*b(1,i)
    enddo
    !
    return
    end function

    !**** new **********************************************************************
    subroutine btdb(elstif,b,db,nee,nrowb,nstr)
    !
    !.... program to multiply b(transpose) * db taking account of symmetry
    !        and accumulate into element stiffness matrix
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
    integer :: nee,nrowb,nstr
    real*8, external :: coldot
    !
    integer :: i,j
    !
    do 200 j=1,nee
        !
        do 100 i=1,j
            elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)
100     continue
        !
200 continue
    !
    return
    end subroutine
    !
    !----------------------------------------------------------------------
    !
    subroutine solverGaussBanda(s,x,r,ns,lb)
    implicit none

    integer :: ns,lb
    real*8 :: s(ns,lb), r(ns), x(ns)
    !
    integer :: l1,l2,ll,n,i,j,k,l,ls
    real*8 :: c
    real*8 :: part, ti, tf, tfAnt, tempoProc, tempoProcTotEst, somaTempoProcTotEst
    integer :: cont

    write(*,*) "..... Resolvendo o Sistema Pressao 3D, solver interno"

    part=0.02
    if(ns > 50000)  part=0.001
    !       print* , ' nt(part*ns) =', int(part*ns)

    call timing(ti)
    tfAnt=ti
    cont = 0
    somaTempoProcTotEst = 0.0

    ls=(lb-1)/2
    l1=ls+1
    l2=ls+2
    do n=1,ns
        !         call estimativasDeDesempenho()
        i=n
        do l=2,l1
            i=i+1
            if(i.gt.ns) cycle ! inicia nova iteracao, proximo l
            ll=n-i+l1
            if(s(i,ll).ne.0.) then
                c=s(i,ll)/s(n,l1)
                j=ll-1
                do k=l1,lb
                    j=j+1
                    if(s(n,k).eq.0.) cycle ! inicia nova iteracao, proximo k
                    s(i,j)=s(i,j)-c*s(n,k)
                end do
                r(i)=r(i)-c*r(n)
            end if
        end do
    end do
    r(ns)=r(ns)/s(ns,l1)
    n=ns - 1
    do n = ns - 1, 1, -1
        !        if(mod(n,int(0.05*ns)) == 0) write(*,*)'2, n = ', n
        l=n
        do k=l2,lb
            l=l+1
            if(l.gt.ns) exit ! finaliza o laco de repeticao, k > lb
            r(n)=r(n)-s(n,k)*r(l)
        end do
        r(n)=r(n)/s(n,l1)
    end do
    x=r

    contains
    subroutine estimativasDeDesempenho()

    if(mod(n,int(part*ns)) == 0) then
        cont = cont + 1
        call timing(tf)
        tempoProc           = tf - tfAnt
        tempoProcTotEst     = ns * tempoProc/(part*ns)
        somaTempoProcTotEst = somaTempoProcTotEst + tempoProcTotEst
        write(*,*)' n = ', n
        write(*,*) int(part*ns), 'iteracoes, tempo =', tempoProc
        write(*,*)' tempo total     (inst.) estimado =', tempoProc *      ns  / (part*ns)
        write(*,*)' tempo decorrido (inst.) estimado =', tempoProc *       n  / (part*ns)
        write(*,*)' tempo restante  (inst.) estimado =', tempoProc * (ns - n) / (part*ns)
        !         write(*,*)' tempo total (medio)       estimado =', somaTempoProcTotEst/cont
        !         write(*,*)' tempo restante            estimado =', somaTempoProcTotEst/cont -  (tf - ti)
        !write(*,*)' tempo de estimado total =' , ns * tempoProc/(part*ns)
        !write(*,*)' tempo restante          =' , ns * tempoProc/(part*ns) - (tf - ti)
        tfAnt=tf
    end if
    end subroutine estimativasDeDesempenho
    end  subroutine

    !
    !=======================================================================
    !
    subroutine solverUMFPack(alhs, b, Ap, Ai, neq, nonzeros)

    implicit none
    !
    real*8  :: alhs(*), b(neq)
    integer :: Ap(*), Ai(*)
    integer, intent(in) :: neq, nonzeros
    !
    real*8, allocatable  :: x(:)
    !
    allocate(x(neq)); x=0.0d0

#ifdef withumfpack
    call solverUMFPackPPD(Ap, Ai, alhs, x, b, neq, nonzeros)
#endif

    b = x
    deallocate(x)

    end subroutine solverUMFPack
    !
    !=======================================================================
    !
    subroutine solverUMFPackPPD(Ap, Ai, Ax, x, b,  neq, nonzeros)

    implicit none

    real*8,   intent(inout) :: Ax(*)
    integer,  intent(inout) :: Ap(*), Ai(*)
    real*8,   intent(inout) :: x(*)
    real*8,   intent(in)    :: b(*)
    integer, intent(in) :: neq, nonzeros
    !
    real*8  :: control (20), info (90)
    integer :: j, p, symbolic, numeric, sys

#ifdef withumfpack
    !     ----------------------------------------------------------------
    !     convert from 1-based to 0-based
    !     ----------------------------------------------------------------

    do 60 j = 1, neq+1
        Ap (j) = Ap (j) - 1
60  continue
    do 70 p = 1, nonzeros
        Ai (p) = Ai (p) - 1
70  continue

    !     ----------------------------------------------------------------
    !     factor the matrix and save to a file
    !     ----------------------------------------------------------------

    !     set default parameters
    call umf4def (control)

    !     print control parameters.  set control (1) to 1 to print
    !     error messages only
    control (1) = 1
    call umf4pcon (control)

    !     pre-order and symbolic analysis
    call umf4sym (neq, neq, Ap, Ai, Ax, symbolic, control, info)

    !   check umf4sym error condition
    !   if (info (1) .lt. 0) then
    !       print *, 'Error occurred in umf4sym: ', info (1)
    !       stop
    !   endif

    !     numeric factorization
    call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)

    !     print statistics for the numeric factorization
    !     call umf4pinf (control, info) could also be done.
    !     print 90, info (1), info (66), &
    !         (info (41) * info (4)) / 2**20, &
    !         (info (42) * info (4)) / 2**20, &
    !         info (43), info (44), info (45)
90  format ('numeric factorization:',/, &
        '   status:  ', f5.0, /, &
        '   time:    ', e10.2, /, &
        '   actual numeric LU statistics:', /, &
        '   size of LU:    ', f10.2, ' (MB)', /, &
        '   memory needed: ', f10.2, ' (MB)', /, &
        '   flop count:    ', e10.2, / &
        '   nnz (L):       ', f10.0, / &
        '   nnz (U):       ', f10.0)

    !     check umf4num error condition
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4num: ', info (1)
        stop
    endif
    !     solve Ax=b, without iterative refinement
    sys = 0
    call umf4sol (sys, x, b, numeric, control, info)
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4sol: ', info (1)
        stop
    endif

    !     print the residual.  x (i) should be 1 + i/n
    !     call resid (neq, nonzeros, Ap, Ai, Ax, pressaoElem, brhsP, r)

    !     ----------------------------------------------------------------
    !     solve Ax=b, with iterative refinement
    !     ----------------------------------------------------------------

    sys = 0
    call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info)
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4solr: ', info (1)
        stop
    endif

    !     print the residual.  x (i) should be 7 - i/n
    !     call resid (neq, nonzeros, Ap, Ai, Ax, x, b, r)

    !     free the numeric factorization
    call umf4fnum (numeric)
    !     free the symbolic analysis
    call umf4fsym (symbolic)

    do j = 1, neq+1
        Ap (j) = Ap (j) + 1
    enddo
    do p = 1, nonzeros
        Ai (p) = Ai (p) + 1
    enddo
    !
#endif
    end subroutine solverUMFPackPPD
    !         LNCC/MCT
    !         Petropolis, 07.2013
    !
    !**** NEW **** MODIFIED FOR IRREGULAR MESH  *****************
    !
    subroutine printf(f,ndof,numnp,nlv,iecho)
    !
    !.... program to print prescribed force and boundary condition data
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer ndof, numnp, nlv, iecho
    real*8 :: f(ndof,numnp,*)
    !
    logical lzero
    integer :: nn, n, i
    !
    nn = 0
    !
    do 100 n=1,numnp
        call ztest(f(1,n,nlv),ndof,lzero)
        if (.not.lzero) then
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
            write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
        endif
100 continue
    !
    return
    !
1000 format('1',&
        ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
        '  b o u n d a r y   c o n d i t i o n s'//5x,&
        ' load vector number = ',i10///5x,&
        ' node no.',6(13x,'dof',i1,:)/)
2000 format(6x,i10,10x,6(1pe15.8,2x))
    end subroutine
    !
    !**** new ********************************************************************
    !
    subroutine printd(name,dva,ndof,numnp,icode)
    !
    !.... program to print kinematic data
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer :: ndof,numnp, icode
    character (LEN=*) ::  name
    real*8 dva(ndof,*)
    !
    logical lzero
    integer nn, n, i
    !
    nn = 0
    !
    do 100 n=1,numnp
        !       call ztest(dva(1,n),ndof,lzero)
        !       if (.not.lzero) then
        !          nn = nn + 1
        !          if (mod(nn,50).eq.1) &
        !            write(icode,1000) name,(i,i=1,ndof)
        write(icode,2000) n,(dva(i,n),i=1,ndof)
        !       endif
100 continue
    !
    return
    !
1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
2000 format(1x,i10,2x,6(1pe30.10,2x))
    end subroutine

    !
    !=============================================================================
    !
    subroutine ztest(a,n,lzero)
    !
    !.... program to determine if an array contains only zero entries
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    integer, intent(in)    :: n
    real*8, intent(in)     :: a(n)
    logical, intent(inout) :: lzero
    !
    integer :: i
    !
    lzero = .true.
    !
    do 100 i=1,n
        if (a(i).ne.0.0d0) then
            lzero = .false.
            return
        endif
100 continue
    !
    end
    !
    !**** new *******************************************************************
    !
    subroutine timing(tempoDeParede)
    implicit none
    !
    !.... program to determine elapsed cpu time
    !
    !
    real*8, intent(inout) :: tempoDeParede
    !
    character(LEN=8)  :: date
    character(LEN=10) :: minhaHora
    character(LEN=5)  :: zone
    integer,dimension(8) :: values
    integer ::  horas, minutos, segundos, milesimosSeg
    !
    call date_and_time(date,minhaHora,zone,values);
    !
    horas=values(5);      minutos=values(6);
    segundos=values(7); milesimosSeg=values(8);
    !
    tempoDeParede = (60*horas+minutos)*60.0+segundos+milesimosSeg/1000.00
    !
    return
    end

    !
    !=======================================================================
    !
    subroutine clear(a,m)
    !
    !.... program to clear a floating-point array
    !
    implicit real*8 (a-h,o-z)
    !
    !.... remove above card for single-precision operation
    !
    INTEGER :: M
    REAL*8 :: a(*)
    integer :: i
    !
    do 100 i=1,m
        a(i) = 0.0d0
100 continue
    !
    return
    end subroutine
    !
    !======================================================================
    !
    subroutine dividirTrabalho   (posI_, posF_, nth_, id_, inicio_, fim_)
    implicit none
    !
    integer, intent(in)   :: posI_, posF_, nth_, id_
    integer, intent(out)  :: inicio_, fim_
    !
    integer :: tamanhoBloco, resto, elementos
    LOGICAL :: distribuirResto = .true.

    tamanhoBloco = int((posF_-posI_+1)/nth_)
    resto = mod((posF_-posI_+1),nth_)
    inicio_ = posI_ + id_*tamanhoBloco
    fim_    = inicio_ + tamanhoBloco - 1

    if(resto > 0) then
        if (distribuirResto) then
            if(id_<resto) then
                if(id_ /= 0) inicio_ = inicio_+ id_
                fim_    = fim_ + id_ + 1
            else
                inicio_= inicio_ + resto
                fim_   = fim_    + resto
            end if
        else
            if(id_ == nth_ - 1) fim_ = fim_ + resto
        end if
    end if
    elementos = fim_ - inicio_ + 1
    ! write(*,'(a,4(a,i0))')rotulo,', id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos
    !  write(*,'(4(a,i0))')' id=', id_,  ' em dividirTrabalho, ' , inicio_,' ',  fim_, ' ', elementos

    end subroutine dividirTrabalho
    !
    !**** new **********************************************************************
    !
    subroutine gerarLabel(label,tempo)
    !
    implicit none

    character(LEN=21), intent(out) :: label
    real*8, intent(in) :: tempo
    !
    character(LEN=21) :: labelAux, num
    integer :: i

    write(num,'(f20.4)') tempo
    !     labelAux="t="//ADJUSTL(num)
    labelAux="t="//num
    do i = 1, 21
        if(labelAux(i:i) .ne. ' ') then
            label(i:i) = labelAux(i:i)
        else
            label(i:i) = '0'
        end if
    end do

    end subroutine

    !
    !**** new **********************************************************************
    !
    subroutine diag(idiag,neq,n)
    !
    implicit none
    !
    !.... program to compute diagonal addresses of left-hand-side matrix
    !
    integer :: neq, n
    integer :: idiag(neq)
    !
    integer :: i

    n = 1
    idiag(1) = 1
    if (neq.eq.1) return
    !
    do 100 i=2,neq
        idiag(i) = idiag(i) + idiag(i-1) + 1
100 continue
    n = idiag(neq)
    !
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine load(id,f,brhs,ndof,numnp,nlvect)
    !
    !.... program to accumulate nodal forces and transfer into
    !        right-hand-side vector
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: f(ndof,numnp,*),brhs(*)
    integer :: ndof, numnp, nlvect
    !
    integer :: nlv
    integer :: i, j, k
    !
    do 300 i=1,ndof
        !
        do 200 j=1,numnp
            k = id(i,j)
            if (k.gt.0) then
                !
                do 100 nlv=1,nlvect
                    brhs(k) = brhs(k) + f(i,j,nlv)
100             continue
                !
            endif
            !
200     continue
        !
300 continue
    !
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    SUBROUTINE LOADTIME(id,f,brhs,ndof,numnp,nlvect,XTIME)
    !
    !.... program to accumulate nodal forces and transfer into
    !        right-hand-side vector
    !
    IMPLICIT NONE
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: f(ndof,numnp,*),brhs(*)
    integer :: ndof, numnp, nlvect
    !
    integer :: nlv
    integer :: i, j, k
    REAL*8  :: XTIME
    !
    DO 300 I=1,NDOF
        DO 200 J=1,NUMNP
            K = ID(I,J)
            IF (K.GT.0) THEN
                DO 100 NLV=1,NLVECT
                    BRHS(K) = BRHS(K) + F(I,J,NLV)*XTIME
100             CONTINUE
            ENDIF
200     CONTINUE
300 CONTINUE
    !
    RETURN
    END SUBROUTINE
    !
    !**** new **********************************************************************
    subroutine btod(id,d,brhs,ndof,numnp)
    !
    !.... program to perform transfer from r.h.s. to displacement array
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: i, j, k
    integer :: ndof, numnp
    integer, dimension(ndof,numnp) :: id
    real*8,  dimension(ndof,numnp) :: d
    real*8,  dimension(*)          :: brhs
    !
    do 200 i=1,ndof
        !
        do 100 j=1,numnp
            k = id(i,j)
            if (k.gt.0) then
                d(i,j) = brhs(k)
            end if
100     continue
        !
200 continue
    !
    return
    end subroutine


    !**** new **********************************************************************
    subroutine kdbc(eleffm,elresf,dl,nee)
    !
    !.... program to adjust load vector for prescribed displacement
    !     boundary condition
    use mGlobaisEscalares, only: zero
    !
    implicit real*8 (a-h,o-z)
    !
    !.... remove above card for single-precision operation
    !
    integer :: nee
    real*8  :: eleffm(nee,*),elresf(*),dl(*)
    !
    integer :: i,j
    real*8  :: val
    !
    do 200 j=1,nee
        !
        val=dl(j)
        if(val.eq.zero) go to 200
        !
        do 100 i=1,nee
            elresf(i)=elresf(i)-eleffm(i,j)*val
100     continue
        !
200 continue
    !
    return
    end subroutine

    !**** new **********************************************************************
    subroutine pivots(a,idiag,neq,nsq,iecho,*)
    !
    !.... program to determine the number of zero and negative terms in
    !        array d of factorization a = u(transpose) * d * u
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(*)
    integer :: idiag(*)
    integer :: neq, nsq, iecho
    !
    integer :: iz, in, n, i
    !
    iz = 0
    in = 0
    !
    do 100 n=1,neq
        i = idiag(n)
        if (a(i).eq.0.) iz = iz + 1
        if (a(i).lt.0.) in = in + 1
100 continue
    !
    write(iecho,1000) nsq,iz,in
    !
    return 1
    !
1000 format(' ',&
        ' zero and/or negative pivots encountered                ', ///5x,&
        ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
        ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
        ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)
    !
    end subroutine

    !**** new **********************************************************************
    subroutine ftod(id,d,f,ndof,numnp,nlvect)
    !
    !.... program to compute displacement boundary conditions
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: d(ndof,*),f(ndof,numnp,*)
    integer :: ndof, numnp, nlvect
    !
    integer :: i, j, k, lv
    real*8  :: val
    !
    !
    do 300 i=1,ndof
        !
        do 200 j=1,numnp
            !
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
            do 100 lv=1,nlvect
                val = val + f(i,j,lv)
100         continue
            !
            d(i,j) = val
            !
200     continue
        !
300 continue
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine ftodTIME(id,d,f,ndof,numnp,nlvect,XTIME)
    !
    !.... program to compute displacement boundary conditions
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: id(ndof,*)
    real*8  :: d(ndof,*),f(ndof,numnp,*)
    integer :: ndof, numnp, nlvect
    !
    integer :: i, j, k, lv
    real*8  :: val, XTIME
    !
    do 300 i=1,ndof
        do 200 j=1,numnp
            k = id(i,j)
            if (k.gt.0) go to 200
            val = 0.0d0
            do 100 lv=1,nlvect
                val = val + f(i,j,lv)*XTIME
100         continue
            d(i,j) = val
200     continue
300 continue
    !
    return
    !
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine colht(idiag,lm,ned,nen,numel,neq)
    !
    !.... program to compute column heights in global left-hand-side matrix
    !
    implicit none
    integer :: idiag(*),lm(ned,nen,*)
    integer :: ned, nen, numel, neq
    !
    integer :: i, j, k
    integer :: m, min, num
    !
    do 500 k=1,numel
        min = neq
        !
        do 200 j=1,nen
            !
            do 100 i=1,ned
                num = lm(i,j,k)
                if (num.gt.0) min = min0(min,num)
100         continue
            !
200     continue
        !
        do 400 j=1,nen
            do 300 i=1,ned
                num = lm(i,j,k)
                if (num.gt.0) then
                    m = num - min
                    if (m.gt.idiag(num)) then
                        idiag(num) = m
                    endif
                endif
                !
300         continue

400     continue
        !
500 continue
    !
    return
    end subroutine

    !**** new **********************************************************************
    function coldot(a,b,n)
    !
    !.... program to compute the dot product of vectors stored column-wise
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(*),b(*)
    integer :: n
    !
    real*8  :: coldot
    integer :: i
    !

    coldot = 0.0d0
    !
    !       do 100 i=1,n
    !       coldot = coldot + a(i)*b(i)
    coldot=dot_product(a(1:n),b(1:n))

    !   100 continue

    !
    return
    end function
    !
    !**** new **********************************************************************
    !
    subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
    !
    !.... program to add rectangular matrices
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(ma,*),b(mb,*),c(mc,*)
    integer :: ma,mb,mc,m,n,iopt
    !
    integer :: i,j
    !
    go to (1000,2000,3000),iopt
    !
    !.... iopt = 1, add entire matrices
    !
1000 do 1200 j=1,n
        !
        do 1100 i=1,m
            c(i,j) = a(i,j) + b(i,j)
1100    continue
        !
1200 continue
    return
    !
    !.... iopt = 2, add lower triangular and diagonal elements
    !
2000 do 2200 j=1,n
        !
        do 2100 i=j,m
            c(i,j) = a(i,j) + b(i,j)
2100    continue
        !
2200 continue
    return
    !
    !.... iopt = 3, add upper triangular and diagonal elements
    !
3000 do 3200 j=1,n
        !
        do 3100 i=1,j
            c(i,j) = a(i,j) + b(i,j)
3100    continue
        !
3200 continue
    return
    !
    end subroutine

    !**** new **********************************************************************
    function rowdot(a,b,ma,mb,n)
    !
    !.... program to compute the dot product of vectors stored row-wise
    !
    implicit none
    !
    !.... remove above card for single precision operation
    !
    real*8  :: a(ma,*),b(mb,*)
    integer :: ma, mb, n
    !
    real*8  :: rowdot
    integer :: i
    !
    rowdot = 0.0d00
    !
    do i=1,n
        rowdot = rowdot + a(1,i)*b(1,i)
    enddo
    !
    return
    end function

    !**** new **********************************************************************
    subroutine btdb(elstif,b,db,nee,nrowb,nstr)
    !
    !.... program to multiply b(transpose) * db taking account of symmetry
    !        and accumulate into element stiffness matrix
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
    integer :: nee,nrowb,nstr
    real*8, external :: coldot
    !
    integer :: i,j
    !
    do 200 j=1,nee
        !
        do 100 i=1,j
            elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)
100     continue
        !
200 continue
    !
    return
    end subroutine

    !
    !----------------------------------------------------------------------
    !
    subroutine solverGaussBanda(s,x,r,ns,lb)
    implicit none

    integer :: ns,lb
    real*8 :: s(ns,lb), r(ns), x(ns)
    !
    integer :: l1,l2,ll,n,i,j,k,l,ls
    real*8 :: c
    real*8 :: part, ti, tf, tfAnt, tempoProc, tempoProcTotEst, somaTempoProcTotEst
    integer :: cont

    write(*,*) "..... Resolvendo o Sistema Pressao 3D, solver interno"

    part=0.02
    if(ns > 50000)  part=0.001
    !       print* , ' nt(part*ns) =', int(part*ns)

    call timing(ti)
    tfAnt=ti
    cont = 0
    somaTempoProcTotEst = 0.0

    ls=(lb-1)/2
    l1=ls+1
    l2=ls+2
    do n=1,ns
        !         call estimativasDeDesempenho()
        i=n
        do l=2,l1
            i=i+1
            if(i.gt.ns) cycle ! inicia nova iteracao, proximo l
            ll=n-i+l1
            if(s(i,ll).ne.0.) then
                c=s(i,ll)/s(n,l1)
                j=ll-1
                do k=l1,lb
                    j=j+1
                    if(s(n,k).eq.0.) cycle ! inicia nova iteracao, proximo k
                    s(i,j)=s(i,j)-c*s(n,k)
                end do
                r(i)=r(i)-c*r(n)
            end if
        end do
    end do
    r(ns)=r(ns)/s(ns,l1)
    n=ns - 1
    do n = ns - 1, 1, -1
        !        if(mod(n,int(0.05*ns)) == 0) write(*,*)'2, n = ', n
        l=n
        do k=l2,lb
            l=l+1
            if(l.gt.ns) exit ! finaliza o laco de repeticao, k > lb
            r(n)=r(n)-s(n,k)*r(l)
        end do
        r(n)=r(n)/s(n,l1)
    end do
    x=r

    contains
    subroutine estimativasDeDesempenho()

    if(mod(n,int(part*ns)) == 0) then
        cont = cont + 1
        call timing(tf)
        tempoProc           = tf - tfAnt
        tempoProcTotEst     = ns * tempoProc/(part*ns)
        somaTempoProcTotEst = somaTempoProcTotEst + tempoProcTotEst
        write(*,*)' n = ', n
        write(*,*) int(part*ns), 'iteracoes, tempo =', tempoProc
        write(*,*)' tempo total     (inst.) estimado =', tempoProc *      ns  / (part*ns)
        write(*,*)' tempo decorrido (inst.) estimado =', tempoProc *       n  / (part*ns)
        write(*,*)' tempo restante  (inst.) estimado =', tempoProc * (ns - n) / (part*ns)
        !         write(*,*)' tempo total (medio)       estimado =', somaTempoProcTotEst/cont
        !         write(*,*)' tempo restante            estimado =', somaTempoProcTotEst/cont -  (tf - ti)
        !write(*,*)' tempo de estimado total =' , ns * tempoProc/(part*ns)
        !write(*,*)' tempo restante          =' , ns * tempoProc/(part*ns) - (tf - ti)
        tfAnt=tf
    end if
    end subroutine estimativasDeDesempenho
    end  subroutine

    !
    !=======================================================================
    !
    subroutine solverUMFPack(alhs, b, Ap, Ai, neq, nonzeros)

    implicit none
    !
    real*8  :: alhs(*), b(neq)
    integer :: Ap(*), Ai(*)
    integer, intent(in) :: neq, nonzeros
    !
    real*8, allocatable  :: x(:)
    !
    allocate(x(neq)); x=0.0d0

#ifdef withumfpack
    call solverUMFPackPPD(Ap, Ai, alhs, x, b, neq, nonzeros)
#endif

    b = x
    deallocate(x)

    end subroutine solverUMFPack
    !
    !=======================================================================
    !
    subroutine solverUMFPackPPD(Ap, Ai, Ax, x, b,  neq, nonzeros)

    implicit none

    real*8,   intent(inout) :: Ax(*)
    integer,  intent(inout) :: Ap(*), Ai(*)
    real*8,   intent(inout) :: x(*)
    real*8,   intent(in)    :: b(*)
    integer, intent(in) :: neq, nonzeros
    !
    real*8  :: control (20), info (90)
    integer :: j, p, symbolic, numeric, sys

#ifdef withumfpack
    !     ----------------------------------------------------------------
    !     convert from 1-based to 0-based
    !     ----------------------------------------------------------------

    do 60 j = 1, neq+1
        Ap (j) = Ap (j) - 1
60  continue
    do 70 p = 1, nonzeros
        Ai (p) = Ai (p) - 1
70  continue

    !     ----------------------------------------------------------------
    !     factor the matrix and save to a file
    !     ----------------------------------------------------------------

    !     set default parameters
    call umf4def (control)

    !     print control parameters.  set control (1) to 1 to print
    !     error messages only
    control (1) = 1
    call umf4pcon (control)

    !     pre-order and symbolic analysis
    call umf4sym (neq, neq, Ap, Ai, Ax, symbolic, control, info)

    !   check umf4sym error condition
    !   if (info (1) .lt. 0) then
    !       print *, 'Error occurred in umf4sym: ', info (1)
    !       stop
    !   endif

    !     numeric factorization
    call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info)

    !     print statistics for the numeric factorization
    !     call umf4pinf (control, info) could also be done.
    !     print 90, info (1), info (66), &
    !         (info (41) * info (4)) / 2**20, &
    !         (info (42) * info (4)) / 2**20, &
    !         info (43), info (44), info (45)
90  format ('numeric factorization:',/, &
        '   status:  ', f5.0, /, &
        '   time:    ', e10.2, /, &
        '   actual numeric LU statistics:', /, &
        '   size of LU:    ', f10.2, ' (MB)', /, &
        '   memory needed: ', f10.2, ' (MB)', /, &
        '   flop count:    ', e10.2, / &
        '   nnz (L):       ', f10.0, / &
        '   nnz (U):       ', f10.0)

    !     check umf4num error condition
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4num: ', info (1)
        stop
    endif
    !     solve Ax=b, without iterative refinement
    sys = 0
    call umf4sol (sys, x, b, numeric, control, info)
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4sol: ', info (1)
        stop
    endif

    !     print the residual.  x (i) should be 1 + i/n
    !     call resid (neq, nonzeros, Ap, Ai, Ax, pressaoElem, brhsP, r)

    !     ----------------------------------------------------------------
    !     solve Ax=b, with iterative refinement
    !     ----------------------------------------------------------------

    sys = 0
    call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info)
    if (info (1) .lt. 0) then
        print *, 'Error occurred in umf4solr: ', info (1)
        stop
    endif

    !     print the residual.  x (i) should be 7 - i/n
    !     call resid (neq, nonzeros, Ap, Ai, Ax, x, b, r)

    !     free the numeric factorization
    call umf4fnum (numeric)
    !     free the symbolic analysis
    call umf4fsym (symbolic)

    do j = 1, neq+1
        Ap (j) = Ap (j) + 1
    enddo
    do p = 1, nonzeros
        Ai (p) = Ai (p) + 1
    enddo
    !
#endif
    end subroutine solverUMFPackPPD