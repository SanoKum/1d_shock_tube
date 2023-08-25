! 1D shock tube using ROE Solver + MUSCL + 4th order Runge Kutta
! coded by Kumpei Sano. 2022 05 05

module Roe_library
    contains

    subroutine calcE( Q, gam, E)
        implicit none

        real(8)             :: Q(3) , gam
        real(8),intent(out) :: E(3)
        real(8)             :: rho , u , energy , p

        rho      = Q(1)
        u        = Q(2)/rho
        energy   = Q(3)

        p = (gam-1.0)*(energy - (rho*u**2)*0.5)

        E(1) = Q(2)
        E(2) = rho*u**2 + p
        E(3) = (energy+p)*u

    end subroutine calcE

    subroutine calcR( u, H, gam, R )
        implicit none

        real(8)             :: u , gam
        real(8)             :: H   , c
        real(8),intent(out) :: R(3,3)

        c = dsqrt((gam-1.0)*(H-0.5*u**2))

        R(1,1) = 1.0   ; R(1,2) = 1.0      ; R(1,3) = 1.0
        R(2,1) = u-c   ; R(2,2) = u        ; R(2,3) = u+c
        R(3,1) = H-c*u ; R(3,2) = 0.5*u**2 ; R(3,3) = H+c*u

    end subroutine calcR

    subroutine calcR_inv(u, H, gam, R_inv )
        implicit none

        real(8)             :: u , gam
        real(8)             :: H   ,  c
        real(8)             :: b1  ,  b2
        real(8),intent(out) :: R_inv(3,3)

        c = dsqrt((gam-1.0)*(H-0.5*u**2))

        b1 = 0.5*u*u*(gam - 1.0)/(c*c)
        b2 = (gam-1.d0)/(c*c)

        R_inv(1,1) = 0.5*(b1 +u/c) ; R_inv(1,2) = -0.5d0*(1.d0/c +b2*u) ; R_inv(1,3) = 0.5*b2
        R_inv(2,1) = 1.d0-b1       ; R_inv(2,2) = b2*u                  ; R_inv(2,3) = -b2
        R_inv(3,1) = 0.5*(b1 -u/c) ; R_inv(3,2) = 0.5d0*(1.d0/c-b2*u)   ; R_inv(3,3) = 0.5*b2

    end subroutine calcR_inv

    subroutine calcLambda( u_ave, H_ave, Correction, uL, HL, uR, HR, gam, Lm )
        implicit none

        integer             :: Correction ! for Expantion
        real(8)             :: u_ave , gam
        real(8)             :: H_ave   , c_ave
        real(8)             :: uL , HL , cL
        real(8)             :: uR , HR , cR
        real(8),intent(out) :: Lm(3,3)

        c_ave = dsqrt((gam-1.0)*(H_ave-0.5*u_ave**2))
        cL    = dsqrt((gam-1.0)*(HL   -0.5*uL**2))
        cR    = dsqrt((gam-1.0)*(HR   -0.5*uR**2))

        Lm(1:3,1:3) = 0.d0

        if ((Correction == 1) .and. (uL-cL < 0.0) .and. (0.0 < uR-cR)) then 
            Lm(1,1) = dabs(u_ave - c_ave) + 0.5*((uR-cR) - (uL-cL))
        else
            Lm(1,1) = dabs(u_ave - c_ave)
        end if

        if ((Correction == 1) .and. (uL < 0.0) .and. (0.0 < uR)) then 
            Lm(2,2) = dabs(u_ave) + 0.5*(uR - uL)
        else
            Lm(2,2) = dabs(u_ave)
        end if

        if ((Correction == 1) .and. (uL+cL < 0.0) .and. (0.0 < uR+cR)) then 
            Lm(3,3) = dabs(u_ave + c_ave) + 0.5*((uR+cR) - (uL+cL))
        else
            Lm(3,3) = dabs(u_ave + c_ave)
        end if

    end subroutine calcLambda

    subroutine calcRoeAverage(rL,rR,uL,uR,HL,HR , r2,u2,H2)
        implicit none

        real(8),intent(in)  :: rL , rR
        real(8),intent(in)  :: uL , uR
        real(8),intent(in)  :: HL , HR
        real(8),intent(out) :: r2 , u2 , H2

        r2 = dsqrt(rL*rR)
        u2 = (dsqrt(rL)*uL + dsqrt(rR)*uR)/(dsqrt(rL)+dsqrt(rR))
        H2 = (dsqrt(rL)*HL + dsqrt(rR)*HR)/(dsqrt(rL)+dsqrt(rR))

    end subroutine calcRoeAverage

    function calcH( rho, u, e, gam)
        implicit none

        real(8) :: rho , u , e , gam
        real(8) :: p 
        real(8) :: calcH

        p     = (gam-1.0)*(e-(rho*u**2)*0.5)
        calcH = (e + p)/rho

    end function calcH

    function calcC( rho, u, e, gam )
        implicit none

        real(8) :: rho , u , e , gam
        real(8) :: p   , c2, H
        real(8) :: calcC

        p = (gam-1.0)*(e-(rho*u**2)*0.5)
        H = (e + p)/rho
        c2= (gam-1.0)*(H-0.5*u**2)
        if (c2 .le.0.0) then
            calcC = 0.0
        else
            calcC = dsqrt(c2)
        end if
    
    end function calcC

end module

program main
    use Roe_library
    implicit none

    ! basic parameter
    integer :: nCell
    integer :: nSolPt
    integer :: DT_method
    integer :: exp_correction
    logical :: MUSCL
    real(8) :: DT , DTN
    real(8) :: End_time
    integer :: writePeriod
    real(8) :: time
    integer :: iSTEP , LastSTEP 
    real(8) :: xmax , xmin , DX
    real(8) :: CFL
    real(8) :: kappa
    integer :: left  =0
    integer :: center=1
    integer :: right =2


    ! Runge-Kutta 
    integer :: iRK
    integer :: nRungeKutta ! runge kutta order
    real(8),allocatable :: RK_K(:,:,:) !iCell,iSol,k,equation

    ! temp value
    integer :: iCell , iSol
    real(8) :: R      (3,3)
    real(8) :: R_inv  (3,3)
    real(8) :: Lm     (3,3)
    real(8) :: A      (3,3)
    real(8) :: uL , uR , uRoe
    real(8) :: rL , rR , rRoe
    real(8) :: HL , HR , HRoe
    real(8) :: eL , eR
    real(8) :: Q_L(3) , Q_R(3)
    real(8) :: E_L(3) , E_R(3)
    real(8) :: ctemp
    real(8) :: del_pls(3), del_mns(3)
    real(8) :: epsilon
    real(8) :: limiter(3)
    character*20 :: cSTEP
    character*40 :: fname

    ! mesh and variables
    real(8),allocatable   :: x (:,:) ! Cell id, 0:left plane, 1:cell center, 2:right plane

    real(8),allocatable   :: u    (:,:)
    real(8),allocatable   :: uN   (:,:)
    real(8),allocatable   :: rho  (:,:)
    real(8),allocatable   :: rhoN (:,:)
    real(8),allocatable   :: prs  (:,:)
    real(8),allocatable   :: prsN (:,:)
    real(8),allocatable   :: teng (:,:)
    real(8),allocatable   :: tengN(:,:)

    real(8),allocatable   :: int_eng (:,:)

    real(8),allocatable   :: Q      (:,:,:) , QN(:,:,:)
    real(8),allocatable   :: E      (:,:,:) , E2(:,:,:)

    real(8),allocatable   :: temp(:)
    real(8),allocatable   :: limiter_func(:,:)
    real(8),allocatable   :: limiter_func2(:,:)

    ! shock parameter
    real(8),parameter     :: gam     = 1.4d0
    real(8),parameter     :: rhoL    = 1.d0  
    real(8),parameter     :: prsL    = 1.d0 
    real(8),parameter     :: velL    = 0.d0
    real(8),parameter     :: rhoR    = 0.125d0
    real(8),parameter     :: prsR    = 0.1d0
    real(8),parameter     :: velR    = 0.d0

    ! -----------------------
    ! *** initial setting ***
    ! -----------------------
    nCell  = 100

    DT_method = 0     ! 0:constant DT, 1:CFL
    DT        = 2.d-3*100/nCell
    CFL       = 0.4d0 ! Courant Number
    !DT        = 0.0001
    End_time  = 0.2d0

    exp_correction = 0 ! 0: off, 1: on

    LastSTEP = 2000

    writePeriod= nCell/100
    xmax   = 1.d0 ; xmin = 0.d0
    nSolPt = 1

    nRungeKutta = 4
    !nRungeKutta = 4

    
    MUSCL   = .TRUE. ! MUSCL interpolation flag
    !MUSCL   = .FALSE.
        kappa   = 1.d0/3.d0 ! 1/3: 3rd order MUSCL, -1: 2nd upwind
        !kappa   = -1.d0     ! 1/3: 3rd order MUSCL, -1: 2nd upwind
        epsilon = 1e-20

    allocate( x      (nCell,0:2) )
    allocate( u      (nCell,0:2) )
    allocate( uN     (nCell,0:2) )

    allocate( rho    (nCell,0:2) )
    allocate( rhoN   (nCell,0:2) )
    allocate( prs    (nCell,0:2) )
    allocate( prsN   (nCell,0:2) )
    allocate( teng   (nCell,0:2) )
    allocate( tengN  (nCell,0:2) )
    allocate( int_eng(nCell,0:2) )

    allocate( Q      (nCell,0:2,3) )
    allocate( QN     (nCell,0:2,3) )
    allocate( E      (nCell,0:2,3) )
    allocate( E2     (nCell,0:2,3) )

    allocate( RK_K   (nCell,nRungeKutta,3))

    allocate( temp   (nCell) )
    allocate( limiter_func (nCell,3) )
    allocate( limiter_func2(nCell,3) )

    ! ----------------------------
    ! *** set cell coordinates ***
    ! ----------------------------
    DX = (xmax - xmin)/nCell
    do iCell=1,nCell
        x(iCell,0) = (iCell-1)*DX           ! left plane
        x(iCell,1) = (iCell-1)*DX + DX/2.d0 ! cell center
        x(iCell,2) = (iCell  )*DX           ! right plane
    end do

    ! ------------------------------------------------
    ! *** set initial value for shock tube problem ***
    ! ------------------------------------------------
    open(10,file="initVal.dat",status="replace")
    write(10,*) "# x , rho, u , teng , prs"
    do iCell=1,nCell
        if (iCell.le.nCell/2) then
            u   (iCell,center)= velL
            rho (iCell,center)= rhoL
            prs (iCell,center)= prsL
            teng(iCell,center)= prsL/(gam-1.0) + rhoL*velL**2*0.5
        else 
            u   (iCell,center)= velR
            rho (iCell,center)= rhoR
            prs (iCell,center)= prsR
            teng(iCell,center)= prsR/(gam-1.0) + rhoR*velR**2*0.5
        end if
        Q(iCell,center,1) = rho (iCell,center)
        Q(iCell,center,2) = rho (iCell,center)*u(iCell,center)
        Q(iCell,center,3) = teng(iCell,center)
        write(10,*) x(iCell,center),rho(iCell,center),u(iCell,center),teng(iCell,center),prs(iCell,center)
    end do
    close(10)


! -------------------------
! *** Start Calculation ***
! -------------------------
    time = 0.d0

    ! ----------------------------
    ! *** output initial value ***
    ! ----------------------------
    write(cSTEP,'(f8.6)') time
    fname = "results_time_"//trim(adjustl(cSTEP))//".dat"
    open(10,file=fname,status="replace")
    write(10,*) "# x   rho   u   teng   prs   int-energy"
    do iCell=1,nCell
        rho    (iCell,center) = Q(iCell,center,1)
        u      (iCell,center) = Q(iCell,center,2)/rho(iCell,center)
        teng   (iCell,center) = Q(iCell,center,3)
        prs    (iCell,center) = (gam-1.0)*(teng(iCell,center)-0.5*rho(iCell,center)*u(iCell,center)**2)
        int_eng(iCell,center) = prs(iCell,center)/(gam-1.0)/rho(iCell,center)

        write(10,*) x(iCell,center),rho(iCell,center),u(iCell,center),teng(iCell,center), &
                    prs(iCell,center),int_eng(iCell,center)
    end do
    close(10)

    do iSTEP=1,LastSTEP
        ! -------------------------------
        ! *** Set DT by CFL condition ***
        ! -------------------------------
        if (DT_method .eq. 1) then
            DTN = DT
            do iCell=1,nCell
                ctemp       = calcC( Q(iCell,center,1), Q(iCell,center,2)/Q(iCell,center,1), Q(iCell,center,3), gam )
                temp(iCell) = abs(Q(iCell,center,2)/Q(iCell,center,1)) + ctemp
            end do
            DT = CFL*DX/maxval(temp)
            !if (iSTEP.ne.1) then
            !    DT = max(DT, DTN*0.90)
            !    DT = min(DT, DTN*1.10)
            !end if
        end if

        time = time + DT

        write(*,*) "STEP = ", iSTEP, "DT = ", DT, "time = ", time

        do iRK=1,nRungeKutta
        ! --------------------------------
        ! *** Set Q and E at interface ***
        ! --------------------------------
            do iCell=1,nCell

                if (iRK.eq.1) then ! first step for Runge-Kutta
                    QN(iCell,left:right,1:3) = Q(iCell,left:right,1:3)
                end if

                call calcE(Q(iCell,center,1:3),gam,E(iCell,center,1:3))

                if (iCell == 1 .or. iCell == nCell .or. MUSCL .eqv. .FALSE.) then ! no MUSCL
                    E(iCell,left,:)  = E(iCell,center,:)
                    Q(iCell,left,:)  = Q(iCell,center,:)

                    Q(iCell,right,:) = Q(iCell,center,:)
                    E(iCell,right,:) = E(iCell,center,:)

                else ! MUSCL
                    del_pls = Q(iCell+1,center,:) - Q(iCell  ,center,:)
                    del_mns = Q(iCell  ,center,:) - Q(iCell-1,center,:)
                    limiter = (2.d0*del_pls*del_mns + epsilon)/(del_pls**2 + del_mns**2 + epsilon)
                    limiter_func(iCell,1:3) = limiter(1:3)

                    limiter_func2(iCell,1) = (2.d0*del_pls(1)*del_mns(1) + epsilon)/(del_pls(1)**2 + del_mns(1)**2 + epsilon)
                    limiter_func2(iCell,2) = (2.d0*del_pls(2)*del_mns(2) + epsilon)/(del_pls(2)**2 + del_mns(2)**2 + epsilon)
                    limiter_func2(iCell,3) = (2.d0*del_pls(3)*del_mns(3) + epsilon)/(del_pls(3)**2 + del_mns(3)**2 + epsilon)
                    
                    ! ------------
                    ! *** left ***
                    ! ------------
                    Q(iCell,left,:) = Q(iCell,center,:) -0.25d0*limiter*((1.d0-kappa*limiter)*del_pls &
                                    + (1.d0+kappa*limiter)*del_mns)
                    call calcE(Q(iCell,left,:),gam,E(iCell,left,:))

                    ! -------------
                    ! *** right ***
                    ! -------------
                    Q(iCell,right,:) = Q(iCell,center,:) +0.25d0*limiter*((1.d0-kappa*limiter)*del_mns &
                                     + (1.d0+kappa*limiter)*del_pls)
                    call calcE(Q(iCell,right,:),gam,E(iCell,right,:))

                end if
            end do

        ! ---------------------------------------------------
        ! *** set common flux at interface using Roe flux ***
        ! ---------------------------------------------------
            do iCell=1,nCell
                ! ----------------------
                ! *** left interface ***
                ! ----------------------
                if (iCell .eq. 1) then ! left boundary
                    E2(iCell, left, 1:3) = E(iCell, center, 1:3)
                else
                    Q_R(1:3) = Q(iCell  , left , 1:3)
                    E_R(1:3) = E(iCell  , left , 1:3)
                    rR       = Q(iCell  , left , 1)
                    uR       = Q(iCell  , left , 2)/rR
                    eR       = Q(iCell  , left , 3)
                    HR       = calcH(rR,uR,eR,gam)

                    Q_L(1:3) = Q(iCell-1, right, 1:3)
                    E_L(1:3) = E(iCell-1, right, 1:3)
                    rL       = Q(iCell-1, right, 1)
                    uL       = Q(iCell-1, right, 2)/rL
                    eL       = Q(iCell-1, right, 3)
                    HL       = calcH(rL,uL,eL,gam)

                    call calcRoeAverage(rL,rR,uL,uR,HL,HR , rRoe,uRoe,HRoe)
                    call calcR         (uRoe, HRoe, gam, R )
                    call calcR_inv     (uRoe, HRoe, gam, R_inv )
                    call calcLambda    (uRoe, HRoe, exp_correction, uL, HL, uR, HR, gam, Lm )

                    A = matmul(R,Lm)
                    A = matmul(A,R_inv)
                    E2(iCell, left, 1:3) = 0.5*(E_L + E_R) -0.5*matmul(A,(Q_R-Q_L))
                end if

                ! right interface
                if (iCell .eq. nCell) then ! right boundary
                    E2(iCell, 2, 1:3) = E(iCell,1,1:3)
                else
                    Q_R(1:3) = Q(iCell+1, left, 1:3)
                    E_R(1:3) = E(iCell+1, left, 1:3)
                    rR       = Q(iCell+1, left, 1)
                    uR       = Q(iCell+1, left, 2)/rR
                    eR       = Q(iCell+1, left, 3)
                    HR       = calcH(rR,uR,eR,gam)

                    Q_L(1:3) = Q(iCell, right, 1:3)
                    E_L(1:3) = E(iCell, right, 1:3)
                    rL       = Q(iCell, right, 1)
                    uL       = Q(iCell, right, 2)/rL
                    eL       = Q(iCell, right, 3)
                    HL       = calcH(rL,uL,eL,gam)

                    call calcRoeAverage(rL,rR,uL,uR,HL,HR , rRoe,uRoe,HRoe)
                    call calcR         (uRoe, HRoe, gam, R )
                    call calcR_inv     (uRoe, HRoe, gam, R_inv )
                    call calcLambda    (uRoe, HRoe, exp_correction, uL, HL, uR, HR, gam, Lm )

                    A = matmul(R,Lm)
                    A = matmul(A,R_inv)
                    E2(iCell, right, 1:3) = 0.5*(E_L + E_R) -0.5*matmul(A,(Q_R-Q_L))

                end if

            end do

            ! ------------------------
            ! *** Time Integration ***
            ! ------------------------

            if (nRungeKutta .eq. 1) then ! 1st order explicit Euler
                do iCell=1,nCell
                    Q(iCell, center, 1:3) = QN(iCell, center, 1:3) &
                                     - DT/DX*(E2(iCell, right, 1:3) - E2(iCell, left, 1:3))
                end do

            else if (nRungeKutta .eq. 4) then ! 4th order Runge kutta
                do iCell=1,nCell
                    RK_K(iCell, iRK, 1:3) = -(E2(iCell, right, 1:3) - E2(iCell, left, 1:3))/DX
        
                    if      (iRK.eq.1) then
                        Q(iCell, center, 1:3) = QN(iCell, center, 1:3) &
                                                + DT*RK_K(iCell, iRK, 1:3)*0.5
                    else if (iRK.eq.2) then
                        Q(iCell, center, 1:3) = QN(iCell, center, 1:3) &
                                                + DT*RK_K (iCell, iRK, 1:3)*0.5
                    else if (iRK.eq.3) then
                        Q(iCell, center, 1:3) = QN(iCell, center, 1:3) &
                                                + DT*RK_K (iCell, iRK, 1:3)
                    else if (iRK.eq.4) then
                        Q(iCell, center, 1:3) = QN(iCell, center, 1:3)            &
                                              + DT*(   RK_K(iCell,1,1:3) &
                                                    +2*RK_K(iCell,2,1:3) &
                                                    +2*RK_K(iCell,3,1:3) &
                                                    +  RK_K(iCell,4,1:3) )/6
                    else 
                        write(*,*) "Error: something wrong"
                        stop
                    end if
        
                end do
            end if

        end do

        ! --------------------
        ! *** output value ***
        ! --------------------
        if (iSTEP-iSTEP/writePeriod*writePeriod .eq. 0) then
            write(cSTEP,'(f8.6)') time
            fname = "results_time_"//trim(adjustl(cSTEP))//".dat"
            open(10,file=fname,status="replace")
            write(10,*) "# x   rho   u   teng   prs   int-energy"
            do iCell=1,nCell
                rho    (iCell,center) = Q(iCell,center,1)
                u      (iCell,center) = Q(iCell,center,2)/rho(iCell,center)
                teng   (iCell,center) = Q(iCell,center,3)
                prs    (iCell,center) = (gam-1.0)*(teng(iCell,center)-0.5*rho(iCell,center)*u(iCell,center)**2)
                int_eng(iCell,center) = prs(iCell,center)/(gam-1.0)/rho(iCell,center)

                write(10,*) x(iCell,center),rho(iCell,center),u(iCell,center),teng(iCell,center), &
                            prs(iCell,center),int_eng(iCell,center)
            end do
            close(10)
        end if

        if (time .ge. End_time) then
            exit
        end if

    end do

    ! ------------------------
    ! *** output end value ***
    ! ------------------------
    write(cSTEP,'(f8.6)') time
    fname = "endVal.dat"
    open(10,file=fname,status="replace")
    write(10,*) "# x   rho   u   teng   prs   int-energy limiter"
    do iCell=1,nCell
        rho    (iCell,center) = Q(iCell,center,1)
        u      (iCell,center) = Q(iCell,center,2)/rho(iCell,center)
        teng   (iCell,center) = Q(iCell,center,3)
        prs    (iCell,center) = (gam-1.0)*(teng(iCell,center)-0.5*rho(iCell,center)*u(iCell,center)**2)
        int_eng(iCell,center) = prs(iCell,center)/(gam-1.0)/rho(iCell,center)

        write(10,*) x(iCell,center),rho(iCell,center),u(iCell,center),teng(iCell,center),prs(iCell,center), &
                    int_eng(iCell,center),limiter_func(iCell,1),limiter_func(iCell,2),limiter_func(iCell,3),&
                    limiter_func2(iCell,1),limiter_func2(iCell,2),limiter_func2(iCell,3)
    end do
    close(10)

end program main
