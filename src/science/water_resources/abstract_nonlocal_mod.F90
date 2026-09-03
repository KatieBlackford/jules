!******************************COPYRIGHT**************************************
! (c) UK Centre for Ecology & Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE abstract_nonlocal_mod

!------------------------------------------------------------------------------
! Description:
!   Calculate abstraction of water from non-local water sources.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in HYDROLOGY
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!------------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE  !  private scope by default
PUBLIC abstract_nonlocal

!------------------------------------------------------------------------------
! Module parameters.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  water_min = 1.0e-10
    ! A minimum amount of water (either a demand or available water) at or
    ! below which some calculations are not performed (kg). This is introduced
    ! primarily so as to avoid very small values in denominators, but this also
    ! prevents resources being spent calculating tiny fluxes that are
    ! physically insignificant. Ignoring these values can mean that (for
    ! example) a small demand is not met when there is water available, but the
    ! amounts are physically insignificant - e.g. a demand of 1.0e-10 kg (for
    ! the gridbox) which is neglected once per day for 100 years means 4e-6kg
    ! of demand is ignored in total.

CONTAINS

!##############################################################################

SUBROUTINE abstract_nonlocal( global_land_pts, priority_order,                &
                              nonlocal_network, demand_nl,                    &
                              demand_unmet, sw_abstracted,                    &
                              nonlocal_abstracted, sw_avail)

USE jules_water_resources_mod, ONLY:                                          &
  nwater_use, l_prioritise, n_sw_source, n_nonlocal_max 

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Abstract water demand from non-local surface water sources 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  global_land_pts
    ! Number of land points in the full model grid.

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  priority_order(global_land_pts,nwater_use),                                 &
    ! Priorities of water demands at each gridpoint, in order of decreasing
    ! priority. Values are the index in multi-sector arrays.
  nonlocal_network(global_land_pts,n_nonlocal_max)
    ! List of gridbox numbers that nonlocal abstractions can draw from.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  demand_nl(global_land_pts,nwater_use)
    ! Demand for water accumulated over the water resource timestep to be met
    ! from non-local surface water (kg).
!------------------------------------------------------------------------------
! Array arguments with INTENT(IN OUT).
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  demand_unmet(global_land_pts,nwater_use),                                   &
    ! Unmet demands for water (kg).
  sw_abstracted(global_land_pts,n_sw_source),                                 &
    ! Net abstraction of water from surface water (kg).
  nonlocal_abstracted(global_land_pts,n_sw_source),                           &
    ! Net abstraction of water from non-local surface water (kg).
  sw_avail(global_land_pts,n_sw_source)
    ! Surface water that is available for abstraction (kg).

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ABSTRACT_NONLOCAL'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                    &
  j, k, l, p, s, i
    ! Loop counters and indices.

REAL(KIND=real_jlslsm) ::                                                     &
  tot_sw_demand,                                                              &
    ! Total mass of water to be abstracted from surface water (kg).
  sw_abs,                                                                     &
    ! Amount of surface water extracted (kg).
  tot_sw_avail
    ! Total mass of available surface water (kg).

!------------------------------------------------------------------------------
! Local array variables.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  sw_avail_start(global_land_pts,n_sw_source),                                &
    ! Surface water that is available for abstraction, at the start of this
    ! call (kg).
  demand(nwater_use)
    ! Updating demand (kg).

!------------------------------------------------------------------------------
!end of header

! Save the initial amounts of available water.
sw_avail_start(:,:) = sw_avail(:,:)

DO l = 1, global_land_pts

  demand(:) = demand_nl(l,:)
  
  ! Loop through nonlocal network.
  DO p = 1, n_nonlocal_max
    i = nonlocal_network(l,p)
    IF ( i == imdi ) CYCLE

    ! Calculate the total demand for abstraction from surface water.
    tot_sw_demand = SUM( demand(:) )
    ! Calculate the total available surface water.
    tot_sw_avail  = SUM( sw_avail(i,:) )

    ! If there is negligible demand or negligible available water, move to
    ! the next gridbox in the network. (This also avoids a small value in a
    ! denominator below.)
    IF ( tot_sw_demand > water_min .AND. tot_sw_avail > water_min ) THEN

      IF ( l_prioritise ) THEN
        ! Demands are prioritised.
        DO j = 1, nwater_use
          k = priority_order(l,j)

          ! Loop over sources
          IF ( demand(k) < water_min ) CYCLE

          DO s = 1, n_sw_source
            IF ( demand(k) < sw_avail(i,s) ) THEN
              ! Demand from this sector can be met in full from this 
              ! source.
              sw_abs    = demand(k)
              demand(k) = 0.0
            ELSE
              ! Demand from this sector cannot be met in full from this 
              ! source.
              sw_abs    = sw_avail(i,s)
              demand(k) = demand(k) - sw_abs
            END IF
            demand_unmet(l,k) = demand_unmet(l,k) - sw_abs
            sw_avail(i,s)     = sw_avail(i,s) - sw_abs
          END DO  !  sources

        END DO  !  j (water uses)

      ELSE

        ! Demands are not prioritised.
        IF ( tot_sw_demand <= tot_sw_avail ) THEN
          ! There is enough surface water for all demands.
          demand_unmet(l,:) = demand_unmet(l,:) - demand(:)
          demand(:)         = 0.0

          DO s = 1, n_sw_source
            IF ( tot_sw_demand < sw_avail(i,s) ) THEN
              ! The total demand can be met from this source.
              sw_abs        = tot_sw_demand
              tot_sw_demand = 0.0
            ELSE
              ! The total demand cannot be met from this source. 
              ! Abstract the remaining available water. We will meet the
              ! rest of the demand from other sources.
              sw_abs        = sw_avail(i,s)
              tot_sw_demand = tot_sw_demand - sw_abs
            END IF
            sw_avail(i,s) = sw_avail(i,s) - sw_abs
          END DO  !  sources

        ELSE
          ! There is insufficient water to meet all demands.
          ! Meet a fraction of each demand
          demand_unmet(l,:) = demand_unmet(l,:) - demand(:) *                 &
                              tot_sw_avail / tot_sw_demand
          demand(:)  = demand(:) * (1.0 - tot_sw_avail / tot_sw_demand)
          ! All sources of water have been exhausted
          sw_avail(i,:) = 0.0

        END IF  !  tot_sw_demand v. tot_sw_avail

      END IF  !  l_prioritise

    END IF  !  tot_sw_demand and tot_sw_avail > water_min

  END DO  !  nonlocal_network loop

END DO  !  global_land_pts loop

! Calculate the abstraction from each surface water source
sw_abstracted(:,:) = sw_abstracted(:,:) + sw_avail_start(:,:) - sw_avail(:,:)
nonlocal_abstracted(:,:) = sw_avail_start(:,:) - sw_avail(:,:)

RETURN
END SUBROUTINE abstract_nonlocal

!##############################################################################

END MODULE abstract_nonlocal_mod

