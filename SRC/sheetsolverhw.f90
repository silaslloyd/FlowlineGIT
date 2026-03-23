!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/*****************************************************************************/
! *
! * A prototype solver for advection-diffusion-reaction equation,
! * This equation is generic and intended for education purposes
! * but may also serve as a starting point for more complex solvers.
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE SheetSolverhw( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables (i.e. for a given element)
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: i, j, m, n, nb, nd, t, active, dimsheet   ! n = number of nodes (given element), n = number of degrees of freedom (given element), nd = number of DOFs on the boundary (given element), t = index for..., active = number of active element
  INTEGER :: iter, maxiter  ! nonlinear iteration number, max number of nonlinear iterations
  LOGICAL :: Found, Newton, BulkUpdate, RHSUpdate
  TYPE(Mesh_t), POINTER :: Mesh 
  TYPE(ValueList_t), POINTER :: BodyForce, Material, SolverParams 
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: hwSol, qwSol ! These all 'point' to memory with different names elsewhere
  REAL(KIND=dp), POINTER :: hw(:), qw(:)
  REAL(KIND=dp), ALLOCATABLE :: hwOld(:)  ! copy of values to retain when pointer is overwritten by new iteration values
  INTEGER, POINTER :: hwPerm(:), qwPerm(:)   ! Used to match up node number with solution value at that node (not obvious due to how Elmer stores values)
  INTEGER, ALLOCATABLE :: hwOldPerm(:)  ! copy of values to retain when pointer is overwritten by new iteration values
  TYPE(Nodes_t) :: ElementNodes
  INTEGER :: dim, qwNDOFs, k, rankA, rankM
  REAL(KIND=dp), ALLOCATABLE :: nodalhw(:), dhwdx(:,:), gradPhi0(:,:), dBasisdx(:,:), nodalhwOld(:), dhwdxOld(:,:), nodalqw(:,:)
  REAL(KIND=dp), ALLOCATABLE :: DensityWater(:), LatentHeat(:), Phi0(:), HydraulicConductivity(:), EffectivePressure(:)
  REAL(KIND=dp), ALLOCATABLE :: dEffectivePressuredx(:), ddEffectivePressuredx(:), dHydraulicConductivitydx(:)
  REAL(KIND=dp), ALLOCATABLE :: q0(:,:), qh(:,:), QQh(:)
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:,:), FORCE(:,:)

SAVE MASS, STIFF, LOAD, FORCE

!PointerToSolver => Solver    ! https://fortran-lang.discourse.group/t/understanding-fortran-pointers/1142

SolverParams => GetSolverParams()    ! Access information (keywords) from the relevant solver section in the sif

! Details of mesh - IS THIS NEEDED HERE?
  Mesh => Solver % Mesh
  dim = Mesh % MeshDim     ! 1, 2 or 3 D
  m = Mesh % NumberOfNodes      ! Number of nodes in the mesh
  n = Mesh % MaxElementNodes    ! Maximum number of nodes that there could be in a single element

! Point to water sheet thickness solution: the values and how those values match up with the nodes
  hwSol => Solver % Variable
  hwPerm  => hwSol % Perm    ! Identifier for which of the values (hw) is the value for a given node
  hw => hwSol % Values   ! Values of the water sheet thickness
  hwOld = hw  ! a copy of the hw values from the previous iteration so I can use them after the FEM system solve (during which the pointer values will be updated)
  hwOldPerm = hwPerm

! Point to water flux solution
  qwSol => VariableGet( Mesh % Variables, 'Water Flux' )
  IF ( ASSOCIATED( qwSol ) ) THEN
    qwPerm => qwSol % Perm
    qwNDOFs = qwSol % DOFs  ! number of vector components of the water flux (i.e. 3) 
    qw => qwSol % Values  ! all components of
  ELSE
  CALL FATAL("sheetsolverhw", "Pointer to Water Flux not associated")
  END IF

  CALL DefaultStart()
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL Info('SheetSolverhw','Sheet solver iteration: '//I2S(iter))

    Newton = GetNewtonActive()  ! logical that is 1 if Newton iterations are to be used (Picard otherwise)

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    
    ! Build matrix arising from the bulk elements
    ! ------------------
    Active = GetNOFActive()   ! number of active elements (i.e. not passive)

    DO t=1,Active   ! for each active element..
      Element => GetActiveElement(t)   ! information about that element
      n  = GetElementNOFNodes()        ! number of nodes
      nd = GetElementNOFDOFs()         ! number of degrees of freedom (DOF)
      nb = GetElementNOFBDOFs()        ! number of DOFs that are on the boundary
      CALL LocalMatrix(  Element, n, nd+nb, hw(hwPerm(Element % NodeIndexes(1:n))) )
      dimsheet =  Element % TYPE % DIMENSION
     ! WRITE(*,*) 'DIMENSION OF SHEET', dimsheet
    END DO

    !CALL DefaultFinishBulkAssembly()
    ! If the second and third entries are TRUE logicals, the matrices and force vector will be saved
    BulkUpdate = .TRUE.   ! mass and stiffness matrices are saved if true
    RHSUpdate = .TRUE.    ! force vector is saved if true
    CALL DefaultFinishBulkAssembly(Solver,BulkUpdate,RHSUpdate)
    
    rankM = RANK(MASS)
    rankA = RANK(STIFF)

    CALL Info('SheetSolverhw','After Bulk Assembly, rank of Matrix: '//I2S(rankM)//', rank of A: '//I2S(rankA)//'& 
       , nd: '//I2S(nd)//', n: '//I2S(n)//', total nodes: '//I2S(m),Level=1)   !
    
    ! NO BC APPLIED! We want the natural BC (no normal flux) at the coldtemp boundary.
    ! At the grounding line, we want effective pressure = 0, but this is applied as a 
    ! Dirichlet BC in the sif (i.e. strongly) so is not imposed here.
    !! Build the matrix arising from the boundary elements
    !! ----------------------
    !Active = GetNOFBoundaryElements()   ! number of elements that have at least one node on the domain boundary
    !DO t=1,Active
    !  Element => GetBoundaryElement(t)
    !  IF(ActiveBoundaryElement()) THEN
    !    n  = GetElementNOFNodes()
    !    nd = GetElementNOFDOFs()
    !    CALL LocalMatrixBC(  Element, n, nd )
    !  END IF
    !END DO
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()
    IF( DefaultConverged() ) EXIT    

  END DO

!---------------------------------------------------
! Calculate water flux 
!---------------------------------------------------
! Do this using whichever form was used for the linearisation (i.e. Newton or Picard)
DO t=1, Solver % NumberOfActiveElements   ! for each active element...
  Element => GetActiveElement(t,Solver)
  n = GetElementNOFNodes(Element)   ! number of nodes in element
  CALL GetElementNodes( ElementNodes )   ! ElementNodes = nodal coordinates (x,y,z) for nodes in the element
  CALL GetParameters(Element, Material, n, DensityWater, LatentHeat, Phi0, EffectivePressure, HydraulicConductivity, &
    dEffectivePressuredx, ddEffectivePressuredx, dHydraulicConductivitydx) ! Get parameter values (at nodes) from sif file

  DO i=1, n   ! ... for each node within the element
    j = Element % NodeIndexes(i)  ! Node index of the node within given element

    ! The different components of the flux when linearised:
    !-------------------------------------------------------
    nodalhw(i) = hw(hwPerm(j))   ! hw at current iteration (just calculated)
    dhwdx(i,1:2) = nodalhw(i) * dBasisdx(i,1:2)   ! grad(hw) at current iteration (just calculated)
    gradPhi0(i,1:2) = Phi0(i) * dBasisdx(i,1:2)  ! where N = Phi0 - Phi
    nodalhwOld = hwOld(hwOldPerm(j))   ! hw at previous iteration
    dhwdxOld(i,1:2) = nodalhwOld(i) * dBasisdx(i,1:2)   ! grad(hw) at previous iteration

    ! Coefficient of hw in flux linearisation
    IF (Newton) THEN
      qh(i,1:2) = dHydraulicConductivitydx(i)*dEffectivePressuredx(i)*(gradPhi0(i,1:2) - &
        dEffectivePressuredx(i)*dhwdx(i,1:2)) - HydraulicConductivity(i)*ddEffectivePressuredx(i)*dhwdx(i,1:2)
    ELSE
      qh(i,1:2) = 0
    END IF

    ! Coefficient of grad(hw) in flux linearisation
    IF (Newton) THEN
      QQh(i) = -HydraulicConductivity(i)*dEffectivePressuredx(i)
    ELSE
      QQh(i) = -HydraulicConductivity(i)*dEffectivePressuredx(i)
    END IF

    ! Order 1 term in flux linearisation
    IF (Newton) THEN
      q0(i,1:2) = HydraulicConductivity(i)*gradPhi0(i,1:2) - qh(i,1:2)
    ELSE
      q0(i,1:2) = HydraulicConductivity(i)*gradPhi0(i,1:2)
    END IF

  ! CALCULATE WATER FLUX (NEW ITERATION) FROM THESE COMPONENTS AND hw, gradhw
  ! AT NEW AND PREVIOUS ITERARION.
  
  ! The pointer qw contains all the components of the Water Flux for all the nodes (global)
  ! stacked together, i.e qx1, qy1, qz1, qx2, qy2, qz2, ...
  ! To reassmeble 
    DO k = 1,2   ! dim and qwNDOFs should be equal
      nodalqw(i,k) = q0(i,k) + qh(i,k)*nodalhwOld(i) + QQh(i)*dhwdxOld(i,k)  ! value of each component (k) of qw at each node (i) within the element
      qw((qwPerm(j)-1)*qwNDOFs+k) = nodalqw(i,k)  ! value of  each global node (qwPerm(j)) 
    END DO

  END DO

END DO 

CALL DefaultFinish()
 
CONTAINS

! Get parametrs from sif file
! -----------------------------------------------------------------------------
  SUBROUTINE GetParameters(Element, Material, n, DensityWater, LatentHeat, Phi0, EffectivePressure, HydraulicConductivity, &
    dEffectivePressuredx, ddEffectivePressuredx, dHydraulicConductivitydx)
! -----------------------------------------------------------------------------
  REAL(KIND=dp) :: DensityWater(n), LatentHeat(n), Phi0(n), EffectivePressure(n), HydraulicConductivity(n), &
    dEffectivePressuredx(n), ddEffectivePressuredx(n), dHydraulicConductivitydx(n)
  INTEGER :: n   ! number of nodes
  LOGICAL :: Found = .FALSE.
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Nodes_t) :: Nodes

  CALL GetElementNodes( Nodes )

  Material => GetMaterial()

  ! Get parameter values from the Material section of the sif file, using the keywords used there
  DensityWater(1:n) = ListGetReal( Material, 'Density Water', n, Element % NodeIndexes, Found)
  LatentHeat(1:n) = ListGetReal( Material, 'Latent Heat Capacity', n, Element % NodeIndexes, Found)
  Phi0(1:n) = ListGetReal( Material, 'Phi0', n, Element % NodeIndexes, Found)  ! where Effective Pressure = Phi0 - Hydraulic Conductivity
  EffectivePressure(1:n) = ListGetReal( Material, 'Effective Pressure', n, Element % NodeIndexes, Found)
  HydraulicConductivity(1:n) = ListGetReal( Material, 'Hydraulic Conductivity', n, Element % NodeIndexes, Found)
  dEffectivePressuredx(1:n) = ListGetReal( Material, 'Effective Pressure First Derivative', n, Element % NodeIndexes, Found)
  ddEffectivePressuredx(1:n) = ListGetReal( Material, 'Effective Pressure Second Derivative', n, Element % NodeIndexes, Found)
  dHydraulicConductivitydx(1:n) = ListGetReal( Material, 'Hydraulic Conductivity First Derivative', n, Element % NodeIndexes, Found) ! as a function of N, not hw

!-------------------------------------------------------------------------------
  END SUBROUTINE GetParameters
!-------------------------------------------------------------------------------



! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, nodalhw )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: DensityWater(n), LatentHeat(n), nodalhw(n), Phi0(n), dhwdx(n,3)
    REAL(KIND=dp) :: gradPhi0(n,3), HydraulicConductivity(n), EffectivePressure(n)
    REAL(KIND=dp) :: dEffectivePressuredx(n), ddEffectivePressuredx(n), dHydraulicConductivitydx(n)
    REAL(KIND=dp) :: DensityWaterAtIP, LatentHeatAtIP, hwAtIP, LoadAtIP, Weight
    REAL(KIND=dp) :: HydraulicConductivityAtIP, EffectivePressureAtIP, dEffectivePressuredxAtIP
    REAL(KIND=dp) :: ddEffectivePressuredxAtIP, dHydraulicConductivitydxAtIP, QQhAtIp
    REAL(KIND=dp) :: q0AtIP(dim), dhwdxAtIP(3), gradPhi0AtIP(3), qhAtIP(3)
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim, rankA, rankM
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Water sheet thickness (from previous iteration, I hope) at integration point
    hwAtIP = SUM( nodalhw(1:n) * Basis(1:n) )   ! hw
    dhwdxAtIP = MATMUL( nodalhw(1:n) , dBasisdx(1:n,1:2))   ! grad(hw)
    !WRITE(*,*) 'DHDX', dhwdxAtIP
    ! Get the basal melt source (from the heat solver) as the source of the hydrology solver
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
       Load(1:n) = GetReal( BodyForce,'Melt Source', Found )   ! Supply of water from melting ice
    
    ! Get material property parameter values (at nodes for this element) from sif file
    CALL GetParameters(Element, Material, n, DensityWater, LatentHeat, Phi0, EffectivePressure, HydraulicConductivity, &
      dEffectivePressuredx, ddEffectivePressuredx, dHydraulicConductivitydx)

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )   ! integration points
    IF( Element % ElementIndex == 1 ) THEN
      CALL Info('SheetSolver','Integration points in 1st element: '//I2S(IP % n),Level=8)
    END IF

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      
      ! Parameters at integration point:
      !---------------------------------
      DensityWaterAtIP = SUM( DensityWater(1:n) * Basis(1:n))
      LatentHeatAtIP = SUM( LatentHeat(1:n) * Basis(1:n))

      ! The different components of the flux when linearised:
      !------------------------------------------------------- 
      gradPhi0AtIP = MATMUL( Phi0(1:n) , dBasisdx(1:n,1:2) )  ! where N = Phi0 - Phi
      HydraulicConductivityAtIP = SUM( HydraulicConductivity(1:n) * Basis(1:n) )
      EffectivePressureAtIP = SUM( EffectivePressure(1:n) * Basis(1:n) )
      dEffectivePressuredxAtIP = SUM( dEffectivePressuredx(1:n) * Basis(1:n) )
      
      ! Higher derivative needed for Newton iteration
      IF (Newton) THEN
        dHydraulicConductivitydxAtIP = SUM( dHydraulicConductivitydx(1:n) * Basis(1:n) )
        ddEffectivePressuredxAtIP = SUM( ddEffectivePressuredx(1:n) * Basis(1:n) )
      END IF

      ! Coefficient of hw in flux linearisation
      IF (Newton) THEN
        !WRITE(*,*) 'NEWTON'
        qhAtIP = dHydraulicConductivitydxAtIP*dEffectivePressuredxAtIP*(gradPhi0AtIP - &
          dEffectivePressuredx*dhwdxAtIP) - HydraulicConductivityAtIP*ddEffectivePressuredxAtIP*dhwdxAtIP
      ELSE
        !WRITE(*,*) 'NOT NEWTON'
        !qhAtIP = 0
        qhAtIP = dHydraulicConductivitydxAtIP*dEffectivePressuredxAtIP*(gradPhi0AtIP - &
          dEffectivePressuredx*dhwdxAtIP) - HydraulicConductivityAtIP*ddEffectivePressuredxAtIP*dhwdxAtIP
          !WRITE(*,*) dHydraulicConductivitydxAtIP*dEffectivePressuredxAtIP*(gradPhi0AtIP)
          !WRITE(*,*) '-', dHydraulicConductivitydxAtIP
          !WRITE(*,*) '-', dEffectivePressuredxAtIP
          !WRITE(*,*) '-', dhwdxAtIP

          !WRITE(*,*) dEffectivePressuredx*dhwdxAtIP
          !WRITE(*,*) '-',dEffectivePressuredx
          !WRITE(*,*) '-',dhwdxAtIP

          !WRITE(*,*) HydraulicConductivityAtIP*ddEffectivePressuredxAtIP*dhwdxAtIP
          !WRITE(*,*) '-',HydraulicConductivityAtIP
          !WRITE(*,*) '-',ddEffectivePressuredxAtIP
          !WRITE(*,*) '-',dhwdxAtIP
          !WRITE(*,*) qhAtIP

      END IF

      ! Coefficient of grad(hw) in flux linearisation
      IF (Newton) THEN
        QQhAtIP = -HydraulicConductivityAtIP*dEffectivePressuredxAtIP
      ELSE
        QQhAtIP = -HydraulicConductivityAtIP*dEffectivePressuredxAtIP
      END IF

      ! Order 1 term in flux linearisation
      IF (Newton) THEN
        q0AtIP = HydraulicConductivityAtIP*gradPhi0AtIP - qhAtIP
      ELSE
        q0AtIP = HydraulicConductivityAtIP*gradPhi0AtIP
      END IF

      Weight = IP % s(t) * DetJ

      DO q=1,nd
        ! Melt source
        ! ------------------------------
        FORCE(q) = FORCE(q) + Weight * LoadAtIP * Basis(q)
        
        ! Flux contribution: q0 term (g0,grad(theta))
        FORCE(q) = FORCE(q) + Weight * DensityWaterAtIP * &
          LatentHeatAtIP * SUM( dBasisdx(q,1:2) * q0AtIP(dim) )
        DO p=1,nd
          ! Flux contribution: hw term (hw*qh,grad(theta))
          ! -----------------------------------
          STIFF(p,q) = STIFF(p,q) - Weight * &
            DensityWaterAtIP * LatentHeatAtIP * &
            SUM(qhAtIP(1:2)*dBasisdx(q,1:2)) * Basis(p)
          !WRITE(*,*) "--------------", p,q
 
          !WRITE(*,*) 'SUM', SUM(qhAtIP(1:2)*dBasisdx(q,1:2))
          !WRITE(*,*) 'qhAtIP',qhAtIP(1:2)
          !WRITE(*,*) 'dBasisdx', dBasisdx(q,1:2)

      
          !WRITE(*,*) 'SIFFPQ',STIFF(p,q)

          ! Flux contribution: grad(hw) term (QQh*grad(hw),grad(theta))
          ! -----------------------------------
          STIFF(p,q) = STIFF(p,q) - Weight * DensityWaterAtIP * &
            LatentHeatAtIP * QQhatIP*SUM(dBasisdx(q,1:2) * dBasisdx(p,1:2))

          ! Time derivative (rho*Lw*dhw/dt,theta):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * DensityWaterAtIP * LatentHeatAtIP * &
            Basis(q) * Basis(p)
        END DO
      END DO
    END DO

    ! Check ranks of siffness and mass matrices are full
    rankA = RANK(STIFF)
    rankM = RANK(MASS)

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    !CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !CALL Info('SheetSolverhw','Rank of Matrix: '//I2S(rankM)//', rank of A: '//I2S(rankA)//', nd: '//I2S(nd)//', n: '//I2S(n),Level=1)

    !WRITE(*,*) 'Mass Matrix size = ', SIZE(MASS,1), ' x ', SIZE(MASS,2)
    !DO i = 1, SIZE(MASS,1)
    !    DO j = 1, SIZE(MASS,2)
    !        WRITE(*,*) 'M', i, j, ':', MASS(i,j)
    !    END DO
    !END DO

   !WRITE(*,*) 'Stiffness Matrix size = ', SIZE(STIFF,1), ' x ', SIZE(STIFF,2)
   ! DO i = 1, SIZE(STIFF,1)
   !     DO j = 1, SIZE(STIFF,2)
   !         WRITE(*,*) 'M', i, j, ':', STIFF(i,j)
   !     END DO
   ! END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
!  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
!    INTEGER :: n, nd
!    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
!    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
!    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
!    LOGICAL :: Stat,Found
!    INTEGER :: i,t,p,q,dim
!    TYPE(GaussIntegrationPoints_t) :: IP

!    TYPE(ValueList_t), POINTER :: BC

!    TYPE(Nodes_t) :: Nodes
!    SAVE Nodes
!------------------------------------------------------------------------------
!    BC => GetBC()
!    IF (.NOT.ASSOCIATED(BC) ) RETURN

!    dim = CoordinateSystemDimension()

!    CALL GetElementNodes( Nodes )
!    STIFF = 0._dp
!    FORCE = 0._dp
!    LOAD = 0._dp

!    Flux(1:n)  = GetReal( BC,'field flux', Found )
!    Coeff(1:n) = GetReal( BC,'robin coefficient', Found )
!    Ext_t(1:n) = GetReal( BC,'external field', Found )

        
    ! Numerical integration:
    !-----------------------
!    IP = GaussPoints( Element )
!    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
!      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
!              IP % W(t), detJ, Basis, dBasisdx )

!      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
!      F = SUM(Basis(1:n)*flux(1:n))

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
!      C = SUM(Basis(1:n)*coeff(1:n))
!      Ext = SUM(Basis(1:n)*ext_t(1:n))

!      DO p=1,nd
!        DO q=1,nd
!          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
!        END DO
!      END DO

!      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
!      TotLen = TotLen + Weight 
!    END DO
!    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
!  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SheetSolverhw
!------------------------------------------------------------------------------