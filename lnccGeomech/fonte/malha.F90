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
    !
    module mMalha
    !
    real*8, allocatable :: x(:,:), xc(:,:)
    integer, allocatable :: conecNodaisElem(:,:), conecLadaisElem(:,:)
    integer, allocatable :: listaDosElemsPorNo(:,:)
    integer, allocatable :: listaDosElemsPorFace(:,:)
    integer :: nsd, numel, numnp, numnpD, nen, numLadosELem
    integer :: numnpReserv, numelReserv, numLadosReserv
    !
    INTEGER :: IrregMesh, Dirichlet, Neumann, I4SeaLoad
    INTEGER :: IDome, IGEOFORM
    REAL(8) :: POSTLINE, SALTLINE, RTOPLINE, RBTTLINE, RIFTLINE
    REAL(8) :: LEFTLINE, RGHTLINE, FRNTLINE, BACKLINE, BASELINE
    !
    public :: criarListaVizinhos, formlm
    public :: genel, genel1, genfl, genfl1
    public :: gensh, gensh1, gensh2, gensh3, igen
    PUBLIC :: FNCPROCESS, FNCINIT, FNCSOLID
    !
    contains

    !
    !=======================================================================
    !
    subroutine criarListaVizinhos(nen,numnp,numel,conecElem,listaDosElemsPorNo)
    implicit none
    !
    !     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
    !
    integer nen,numnp,numel
    integer, dimension(nen,numel)  :: conecElem
    integer, dimension(nen,numnp)  :: listaDosElemsPorNo
    integer :: no,nel,l

    do nel=1, numel
        do l=1, nen
            no=conecElem(l,nel)
            listaDosElemsPorNo(l,no) = nel
        end do
    end do

    end subroutine
    !
    !**** new **********************************************************************
    subroutine formlm (id,conecElem,lm,ndof,ned,nen,numel)
    !
    !.... program to form lm array
    !
    integer :: ndof,ned,nen,numel
    integer :: id(ndof,*),conecElem(nen,*),lm(ned,nen,*)
    !
    integer :: i,j,k,node
    !
    do 300 k=1,numel
        !
        do 200 j=1,nen
            node=conecElem(j,k)
            !
            do 100 i=1,ndof

                lm(i,j,k) = id(i,node)

100         continue
            !
200     continue
        !
300 continue

    !
    return
    end subroutine

    !******************************************************************************
    subroutine genel(conecElem,mat,nen, iin)
    !
    !.... program to read and generate element node and material numbers
    !
    !         conecElem(nen,numel) = element node numbers
    !         mat(numel)     = element material numbers
    !         nen            = number of element nodes (le.27)
    !         n              = element number
    !         ng             = generation parameter
    !         nel(i)         = number of elements in direction i
    !         incel(i)       = element number increment for direction i
    !         inc(i)         = node number increment for direction i
    !
    integer :: nen, iin
    integer :: conecElem(nen,*),mat(*),itemp(27)
    !
    integer :: m,ng,i
    integer :: n,nel(3),incel(3),inc(3)
    !
100 continue
    read(iin,1000) n,m,(itemp(i),i=1,nen),ng

    if (n.eq.0) return
    !       call imove(conecElem(1,n),itemp,nen)
    conecElem(1:nen,n)=itemp(1:nen)
    mat(n)=m
    if (ng.ne.0) then
        !
        !....... generate data
        !
        read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)
        !           write(*,1000) (nel(i),incel(i),inc(i),i=1,3)
        call genel1(conecElem,mat,nen,n,nel,incel,inc)
    endif
    go to 100
    !
1000 format(16i10,10x,14i10)
    !
    end subroutine

    !:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine imove(ia,ib,n)
    !
    !.... program to move an integer array
    !
    integer :: ia(1),ib(*)
    integer :: n
    !
    integer :: i
    !
    do 100 i=1,n
        ia(i)=ib(i)
100 continue

    return
    end subroutine
    !:::: new ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine move(a,b,n)
    !
    !.... program to move a floating-point array
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8 :: a(*),b(*)
    integer :: n
    !
    integer i
    !
    do 100 i=1,n
        a(i) = b(i)
100 continue
    !
    return
    end subroutine
    !
    !******************************************************************************
    !
    subroutine genel1(conecElem,mat,nen, n,nel,incel,inc)
    !
    !.... program to generate element node and material numbers
    !
    integer*4 ::  nen, conecElem(nen,*),mat(*)
    integer*4:: n,nel(3),incel(3),inc(3)

    integer*4:: i,j,k,ii,jj,kk,ie,le,je,ke
    !
    !.... set defaults
    !
    call geneld
    !
    !.... generation algorithm
    !
    ie = n
    je = n
    ke = n
    !
    ii = nel(1)
    jj = nel(2)
    kk = nel(3)
    !
    do 300 k=1,kk
        !
        do 200 j=1,jj
            !
            do 100 i=1,ii
                !
                if (i.ne.ii) then
                    le = ie
                    ie = le + incel(1)
                    call geneli(conecElem(1,ie),conecElem(1,le),inc(1),nen)
                    mat(ie) = mat(le)
                endif
100         continue
            !
            if (j.ne.jj) then
                le = je
                je = le + incel(2)
                call geneli(conecElem(1,je),conecElem(1,le),inc(2),nen)
                mat(je) = mat(le)
                ie = je
            endif
200     continue
        !
        if (k.ne.kk) then
            le = ke
            ke = le + incel(3)
            call geneli(conecElem(1,ke),conecElem(1,le),inc(3),nen)
            mat(ke) = mat(le)
            ie = ke
            je=ke
        endif
300 continue
    !
    return
    contains
    !******************************************************************************
    subroutine geneld
    !
    !.... program to set defaults for element node
    !        and material number generation
    !
    if (nel(1).eq.0) nel(1) = 1
    if (nel(2).eq.0) nel(2) = 1
    if (nel(3).eq.0) nel(3) = 1
    !
    if (incel(1).eq.0) incel(1) = 1
    if (incel(2).eq.0) incel(2) = nel(1)
    if (incel(3).eq.0) incel(3) = nel(1)*nel(2)
    !
    if (inc(1).eq.0) inc(1) = 1
    if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)
    if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)
    !
    return
    end subroutine
    !******************************************************************************
    subroutine geneli(conecElem2,conecElem1,inc,nen)
    !
    !.... program to increment element node numbers
    !
    integer :: inc, nen
    integer:: conecElem1(*),conecElem2(*)
    !
    integer :: i
    !
    do 100 i=1,nen
        if (conecElem1(i).eq.0) then
            conecElem2(i) = 0
        else
            conecElem2(i) = conecElem1(i) + inc
        endif
100 continue
    !
    return
    end  subroutine
    end subroutine genel1
    !******************************************************************************
    subroutine genfl(a,nra,iin)
    use mGlobaisEscalares
    !
    !.... program to read and generate floating-point nodal data
    !
    !         a       = input array
    !         nra     = number of rows in a (le.6)
    !         n       = node number
    !         numgp   = number of generation points
    !         ninc(i) = number of increments for direction i
    !         inc(i)  = increment for direction i
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: nra, iin
    real*8  :: a(nra,*)
    !
    real*8  :: temp(6,20)
    integer :: n,numgp,ninc(3),inc(3)
    integer :: i, j, m, mgen
    !
100 continue
    read(iin,1000) n,numgp,(temp(i,1),i=1,nra)
    !         write(*,1000) n,numgp,(temp(i,1),i=1,nra)

    if (n.eq.0) return
    call move(a(1,n),temp,nra)
    !     a(1:nra,n) = temp(1:nra,1)

    if (numgp.ne.0) then
        do 200 j=2,numgp
            !
            read(iin,1000) m,mgen,(temp(i,j),i=1,nra)
            !          write(*,1000) m,mgen,(temp(i,j),i=1,nra)

            if (mgen.ne.0) call move(temp(1,j),a(1,m),nra)
            !          if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m)
            !
200     continue
        read(iin,2000) (ninc(i),inc(i),i=1,3)
        !          write(*,2000) (ninc(i),inc(i),i=1,3)

        call genfl1(a,nra, temp, n, numgp, ninc, inc)
    endif
    go to 100
    !
    ! 1000 format(2i10,6f10.0)
1000 format(2i10,6f10.0)
2000 format(16i10)
    !
    end  subroutine
    !******************************************************************************
    subroutine genfl1(a,nra, temp, n, numgp, ninc, inc)
    use mGlobaisEscalares
    !
    !.... program to generate floating-point nodal data
    !        via isoparametric interpolation
    !
    !         iopt = 1, generation along a line
    !              = 2, generation over a surface
    !              = 3, generation within a volume
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: nra
    real*8  :: a(nra,*)
    real*8  :: temp(6,20)
    integer :: n,numgp,ninc(3),inc(3)
    !
    real*8  :: sh(20,1), dr, ds, dt, r, s, t
    integer :: iopt, ni, nj, nk, ii, jj, kk, i, j, k
    !
    iopt = 3
    if (ninc(3).eq.0) iopt = 2
    if (ninc(2).eq.0) iopt = 1
    !
    dr = 0.0
    ds = 0.0
    dt = 0.0
    !
    if (ninc(1).ne.0) dr = two/ninc(1)
    if (ninc(2).ne.0) ds = two/ninc(2)
    if (ninc(3).ne.0) dt = two/ninc(3)
    !
    ii = ninc(1)+1
    jj = ninc(2)+1
    kk = ninc(3)+1
    !
    ni = n
    nj = n
    nk = n
    !
    t = -one
    do 300 k=1,kk
        !
        s = -one
        do 200 j=1,jj
            !
            r = -one
            do 100 i=1,ii
                !
                call gensh(r,s,t,sh,numgp,iopt)

                !       resultado=matmul(temp,sh)
                !       a(1:nra,ni)=resultado(1:nra,1)
                call multab(temp,sh,a(1,ni),6,20,nra,numgp,nra,1,1)

                ni = ni + inc(1)
                r = r + dr
100         continue
            !
            nj = nj + inc(2)
            ni = nj
            s = s + ds
200     continue
        !
        nk = nk + inc(3)
        ni = nk
        t = t + dt
300 continue
    !
    return
    end  subroutine
    !**** new **********************************************************************
    subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
    !
    !.... program to multiply two matrices
    !
    !        l = range of dot-product index
    !        m = number of active rows in c
    !        n = number of active columns in c
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(ma,*),b(mb,*),c(mc,*)
    integer :: ma,mb,mc,l,m,n,iopt

    real*8, external :: coldot, rowdot
    !
    integer :: i,j
    !
    go to (1000,2000,3000,4000),iopt
    !
    !.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
    !
1000 do 1200 i=1,m
        !
        do 1100 j=1,n
            c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
1100    continue
        !
1200 continue
    return
    !                                            t
    !.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
    !
2000 do 2200 i=1,m
        !
        do 2100 j=1,n
            c(i,j) = coldot(a(1,i),b(1,j),l)
2100    continue
        !
2200 continue
    return
    !                                                t
    !.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
    !
3000 do 3200 i=1,m
        !
        do 3100 j=1,n
            c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
3100    continue
        !
3200 continue
    return
    !                                            t    t
    !.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
    !
4000 do 4200 i=1,m
        !
        do 4100 j=1,n
            c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
4100    continue
        !
4200 continue

    !
    return
    end subroutine
    !
    !**** NEW ******************************************************************
    !
    FUNCTION FNCPROCESS(TASK)
    !
    !.... PROGRAM TO COMPUTE COEFICIENT FOR TERZAGHIS PROBLEM
    !
    IMPLICIT NONE
    !
    REAL(8)        :: FNCPROCESS
    CHARACTER*25   :: FLAG, TASK
    !
    FLAG = 'TERZAGHI'
    !
    IF (TRIM(FLAG).EQ.TRIM(TASK)) THEN
        FNCPROCESS = 1.0D0
    ELSE
        FNCPROCESS = 3.0D0
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW ******************************************************************
    !
    FUNCTION FNCINIT(TASK)
    !
    !.... PROGRAM TO SETUP INITIALIZATION OF STRESS AND DISPLACEMENTS
    !
    IMPLICIT NONE
    !
    LOGICAL        :: FNCINIT, L1, L2
    CHARACTER*25   :: FLAG1, FLAG2, TASK
    !
    L1 = .FALSE.
    L2 = .FALSE.
    FLAG1 = 'TERZAGHI'
    FLAG2 = 'MANDEL'
    !
    IF (TRIM(FLAG1).EQ.TRIM(TASK)) L1=.TRUE.
    IF (TRIM(FLAG2).EQ.TRIM(TASK)) L2=.TRUE.
    IF (L1.OR.L2) THEN
        FNCINIT = .FALSE.
    ELSE
        FNCINIT = .TRUE.
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !**** NEW ******************************************************************
    !
    FUNCTION FNCSOLID(TASK)
    !
    !.... PROGRAM TO SETUP SOLID MODEL EXAMPLES
    !
    IMPLICIT NONE
    !
    LOGICAL        :: FNCSOLID
    CHARACTER*25   :: FLAG1, TASK
    !
    FLAG1 = 'SOLID_ONLY'
    
    IF (TRIM(FLAG1).EQ.TRIM(TASK)) THEN
        FNCSOLID = .TRUE.
    ELSE
        FNCSOLID = .FALSE.
    ENDIF
    !
    RETURN
    !
    END FUNCTION
    !
    !**** new **********************************************************************
    !
    function rcdot(a,b,ma,n)
    !
    !.... program to compute the dot product of a vector stored row-wise
    !        with a vector stored column-wise
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    real*8  :: a(ma,*),b(*)
    integer :: ma, n
    !
    real*8  :: rcdot
    integer :: i
    !
    rcdot = 0.0
    !
    do 100 i=1,n
        rcdot = rcdot + a(1,i)*b(i)
100 continue
    !
    return
    end function
    !
    !******************************************************************************
    !
    subroutine gensh(r,s,t,sh,numgp,iopt)
    !
    !.... program to call shape function routines
    !        for isoparametric generation
    !
    implicit none
    !
    !.... modify above card for single-precision operation
    !
    real*8  :: r, s, t, sh(*)
    !
    integer :: numgp, iopt
    !
    go to (100,200,300),iopt
    !
100 call gensh1(r,sh,numgp)
    return
    !
200 call gensh2(r,s,sh,numgp)
    return
    !
300 call gensh3(r,s,t,sh,numgp)
    return
    !
    end subroutine
    !******************************************************************************
    subroutine gensh1(r,sh,n)
    use mGlobaisEscalares
    !
    !.... program to compute 1d shape functions
    !        for isoparametric generation
    !
    implicit none
    !
    !.... modify above card(s) for single-precision operation
    !
    real*8  :: r, sh(*)
    integer :: n
    !
    sh(2) = pt5*r
    sh(1) = pt5 - sh(2)
    sh(2) = pt5 + sh(2)
    if (n.eq.3) then
        sh(3) = one - r*r
        sh(1) = sh(1) - pt5*sh(3)
        sh(2) = sh(2) - pt5*sh(3)
    endif
    !
    return
    end subroutine
    !******************************************************************************
    subroutine gensh2(r,s,sh,n)
    use mGlobaisEscalares
    !
    !.... program to compute 2d shape functions
    !        for isoparametric generation
    !
    implicit none
    !
    !.... modify above card for single-precision operation
    !
    real*8  :: r, s, sh(*)
    integer :: n
    !
    real*8  :: r1, r2, r3, s1, s2, s3
    !
    r2 = pt5*r
    r1 = pt5 - r2
    r2 = pt5 + r2
    s2 = pt5*s
    s1 = pt5 - s2
    s2 = pt5 + s2
    sh(1) = r1*s1
    sh(2) = r2*s1
    sh(3) = r2*s2
    sh(4) = r1*s2
    if (n.eq.4) return
    !
    r3 = one - r*r
    s3 = one - s*s
    sh(5) = r3*s1
    sh(6) = s3*r2
    sh(7) = r3*s2
    sh(8) = s3*r1
    sh(1) = sh(1) - pt5*(sh(5) + sh(8))
    sh(2) = sh(2) - pt5*(sh(6) + sh(5))
    sh(3) = sh(3) - pt5*(sh(7) + sh(6))
    sh(4) = sh(4) - pt5*(sh(8) + sh(7))
    !
    return
    end subroutine
    !******************************************************************************
    subroutine gensh3(r,s,t,sh,n)
    use mGlobaisEscalares
    !
    !.... program to compute 3d shape functions
    !        for isoparametric generation
    !
    implicit none
    !
    !.... modify above card for single-precision operation
    !
    real*8  :: r, s, t, sh(*)
    integer :: n
    !
    real*8  :: r1, r2, r3, rs1, rs2, rs3, rs4
    real*8  :: s1, s2, s3, t1, t2, t3
    !
    r2 = pt5*r
    r1 = pt5 - r2
    r2 = pt5 + r2
    s2 = pt5*s
    s1 = pt5 - s2
    s2 = pt5 + s2
    t2 = pt5*t
    t1 = pt5 - t2
    t2 = pt5 + t2
    !
    rs1 = r1*s1
    rs2 = r2*s1
    rs3 = r2*s2
    rs4 = r1*s2
    sh(1) = rs1*t1
    sh(2) = rs2*t1
    sh(3) = rs3*t1
    sh(4) = rs4*t1
    sh(5) = rs1*t2
    sh(6) = rs2*t2
    sh(7) = rs3*t2
    sh(8) = rs4*t2
    if (n.eq.8) return
    !
    r3 = one - r*r
    s3 = one - s*s
    t3 = one - t*t
    sh(17) = t3*rs1
    sh(18) = t3*rs2
    sh(19) = t3*rs3
    sh(20) = t3*rs4
    rs1 = r3*s1
    rs2 = s3*r2
    rs3 = r3*s2
    rs4 = s3*r1
    sh( 9) = rs1*t1
    sh(10) = rs2*t1
    sh(11) = rs3*t1
    sh(12) = rs4*t1
    sh(13) = rs1*t2
    sh(14) = rs2*t2
    sh(15) = rs3*t2
    sh(16) = rs4*t2
    !
    sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
    sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
    sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
    sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
    sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
    sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
    sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
    sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
    !
    return
    end subroutine
    !
    !**** new **********************************************************************
    !
    subroutine igen(ia,m,iin)
    use mGlobaisEscalares
    !
    !.... program to read and generate integer nodal data
    !
    !        ia = input array
    !         m = number of rows in ia
    !         n = node number
    !        ne = end node in generation sequence
    !        ng = generation increment
    !
    integer :: m, ia(m,*), iin
    !
    integer :: ib(m)
    integer :: n, ne, ng
    integer :: i
    !
100 continue
    read(iin,1000) n,ne,ng,(ib(i),i=1,m)

    if (n.eq.0) return

    if (ng.eq.0) then
        ne = n
        ng = 1
    else
        ne = ne - mod(ne-n,ng)
    endif
    !
    do 200 i=n,ne,ng
        !     call imove(ia(1,i),ib,m)
        ia(:,i)=ib
200 continue
    !
    go to 100
    !
1000 format(16i10)
    end subroutine
    !
    !**** new **********************************************************************
    subroutine local(conectElem,x,xl,nen,nrowx,nrowxl)
    !
    !.... program to localize a global array
    !
    !        note: it is assumed nrowxl.le.nrowx
    !
    implicit none
    !
    !.... remove above card for single-precision operation
    !
    integer :: conectElem(*)
    integer :: nrowx, nrowxl, nen
    double precision :: x(nrowx,*),xl(nrowxl,*)
    !
    integer :: i, j, node
    !
    do 200 j=1,nen
        node = conectElem(j)
        !
        do 100 i=1,nrowxl
            xl(i,j)= x(i,node)
100     continue
        !
200 continue
    !
    return
    end subroutine


    end module
