    PROGRAM Main
    USE BasicData
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D
    
    IMPLICIT NONE
    
    OPEN(103,FILE="Grid.dat")
    READ(103,*) id,jd
    CLOSE(103)

    
    ALLOCATE(u(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1),v(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1),&
			 rho(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1),p(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
    ALLOCATE(x(1:id,1:jd),y(1:id,1:jd))  
    ALLOCATE(sxx(1:id,1:jd-1),sxy(1:id,1:jd-1),syx(1:id-1,1:jd),syy(1:id-1,1:jd),vol(1:id-1,1:jd-1))
    ALLOCATE(nxx(1:id,1:jd-1),nxy(1:id,1:jd-1),nyx(1:id-1,1:jd),nyy(1:id-1,1:jd))
	ALLOCATE(sxMod(1:id,1:jd-1),syMod(1:id-1,1:jd))
    
    CALL ReadGrid()
    
    CALL ReadSettings()
    
    CALL Information_write()
    
    CALL Initialization()
    
    CALL CalculateEigen()
    
    CALL Output()
    
    DEALLOCATE(rho,u,v,p)
    DEALLOCATE(x,y,sxx,sxy,syx,syy,vol,nxx,nxy,nyx,nyy,sxMod,syMod)
    DEALLOCATE(WR,WI,VectorL,VectorR)

    WRITE(*,*)
    WRITE(*,"(A80)")"The computation is finished, corresponding files are outputted in Results folder."
    WRITE(*,*)
    PAUSE

    END PROGRAM Main