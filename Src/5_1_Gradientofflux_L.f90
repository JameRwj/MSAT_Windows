    SUBROUTINE Gradientofflux_L(i,j,Face,nx,ny,L,faiL)
    USE BasicData,ONLY:p2,rho,u,v,p,zero,two

    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: L,faiL
    REAL(p2),DIMENSION(4) :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus
    REAL(p2),DIMENSION(4) :: flux1_minus,flux2_minus,flux3_minus,flux4_minus
    REAL(p2),DIMENSION(4) :: UL,UR

    REAL(p2) :: faiL1,faiL2,faiL3,faiL4,faiR1,faiR2,faiR3,faiR4
    REAL(p2) :: UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4
    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    delta = 1.0E-7_p2

    !===========重构===========
    CALL Reconstruct_L(i,j,rho,Face,faiL1,UL1)
    CALL Reconstruct_L(i,j,u,  Face,faiL2,UL2)
    CALL Reconstruct_L(i,j,v,  Face,faiL3,UL3)
    CALL Reconstruct_L(i,j,p,  Face,faiL4,UL4)

    CALL Reconstruct_R(i,j,rho,Face,faiR1,UR1)
    CALL Reconstruct_R(i,j,u,  Face,faiR2,UR2)
    CALL Reconstruct_R(i,j,v,  Face,faiR3,UR3)
    CALL Reconstruct_R(i,j,p,  Face,faiR4,UR4)

    deltaL = delta
    deltaR = zero

    !===========左侧加减扰动计算通量===========
    !第一项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_plus)

    !第一项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_minus)

    !第二项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_plus)

    !第二项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_minus)

    !第三项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_plus)

    !第三项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_minus)

    !第四项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_plus)

    !第四项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_minus)

    !=======求每一个通量梯度=======

    L(:,1) = (flux1_plus-flux1_minus)/two/delta
    L(:,2) = (flux2_plus-flux2_minus)/two/delta
    L(:,3) = (flux3_plus-flux3_minus)/two/delta
    L(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                L(m,n) = -L(m,n)
            END DO
        END DO
    END IF

    faiL(1,1) = faiL1
    faiL(2,2) = faiL2
    faiL(3,3) = faiL3
    faiL(4,4) = faiL4

    END SUBROUTINE Gradientofflux_L