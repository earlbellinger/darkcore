! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras
      
            use star_lib
            use star_def
            use const_def
            use math_lib
            use utils_lib, only: mesa_error
            
            implicit none
            
            ! these routines are called by the standard run_star check_model
            contains
            
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               
               ! this is the place to set any procedure pointers you want to change
               ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
               
               
               ! the extras functions in this file will not be called
               ! unless you set their function pointers as done below.
               ! otherwise we use a null_ version which does nothing (except warn).
               
               s% extras_startup => extras_startup
               s% extras_start_step => extras_start_step
               s% extras_check_model => extras_check_model
               s% extras_finish_step => extras_finish_step
               s% extras_after_evolve => extras_after_evolve
               s% how_many_extra_history_columns => how_many_extra_history_columns
               s% data_for_extra_history_columns => data_for_extra_history_columns
               s% how_many_extra_profile_columns => how_many_extra_profile_columns
               s% data_for_extra_profile_columns => data_for_extra_profile_columns  
               
               s% how_many_extra_history_header_items => how_many_extra_history_header_items
               s% data_for_extra_history_header_items => data_for_extra_history_header_items
               s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
               s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

               s% other_cgrav => other_cgrav
               
            end subroutine extras_controls
            
            subroutine do1_relax_R_center(s, new_Rcenter, ierr)
               ! adjust all lnR's to keep same density for each cell as 1st guess for next model
               type (star_info), pointer :: s
               real(dp), intent(in) :: new_Rcenter ! cm
               integer, intent(out) :: ierr
               real(dp) :: dm, rho, dr3, rp13
               integer :: k
               ierr = 0
               s% R_center = new_Rcenter
               ! adjust lnR's
               rp13 = s% R_center*s% R_center*s% R_center
               do k = s% nz, 1, -1
                  dm = s% dm(k)
                  rho = s% rho(k)
                  dr3 = dm/(rho*four_thirds_pi) ! dm/rho is cell volume
                  s% xh(s% i_lnR,k) = log(rp13 + dr3)*one_third
                  rp13 = rp13 + dr3
               end do
            end subroutine do1_relax_R_center
            
            subroutine dark_core(id, s, ierr)
                integer, intent(in) :: id
                type (star_info), pointer :: s
                integer, intent(out) :: ierr
                
                real(dp) :: G, c_s
                real(dp) :: M_core, R_B
                real(dp) :: core_avg_rho, core_avg_eps, new_core_mass
                
                G       = s% cgrav(s% nz)   ! gravitational constant (cm^3 / g s^2)
                c_s     = s% csound(s% nz)  ! speed of sound         (cm / s)
                
                M_core = s% xtra(1)                 ! core mass (g)
                R_B = s% x_ctrl(1) * 2 * G * M_core / pow(c_s, 2)  ! Bondi radius (cm)
                new_core_mass = M_core / Msun       ! new core mass (Msun)
                core_avg_eps = 0                    ! dark core 
                core_avg_rho = 1 / (4 / 3 * pi) * M_core / pow(R_B, 3) ! average core density (g / cm^3)
                
                s% xtra(2) = R_B
                write(*, *) "Bondi radius:", R_B
                
                call star_relax_core( &
                    id, new_core_mass, s% job% dlg_core_mass_per_step, &
                    s% job% relax_core_years_for_dt, &
                    core_avg_rho, core_avg_eps, ierr)
                
            end subroutine dark_core
            
            
            subroutine extras_startup(id, restart, ierr)
                integer, intent(in) :: id
                logical, intent(in) :: restart
                integer, intent(out) :: ierr
                type (star_info), pointer :: s
                
                ierr = 0
                call star_ptr(id, s, ierr)
                if (ierr /= 0) return
                
                if (s% x_logical_ctrl(1)) then
                    s% xtra(1) = s% job% new_core_mass * Msun
                    call dark_core(id, s, ierr) 
                end if
            end subroutine extras_startup
            
            
            integer function extras_start_step(id)
                integer, intent(in) :: id
                integer :: ierr
                type (star_info), pointer :: s
      
                ierr = 0
                call star_ptr(id, s, ierr)
                if (ierr /= 0) return
                extras_start_step = 0
            end function extras_start_step
            
            
            ! returns either keep_going, retry, or terminate.
            integer function extras_check_model(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               extras_check_model = keep_going
               if (.false. .and. s% star_mass_h1 < 0.35d0) then
                  ! stop when star hydrogen mass drops to specified level
                  extras_check_model = terminate
                  write(*, *) 'have reached desired hydrogen mass'
                  return
               end if
               ! if you want to check multiple conditions, it can be useful
               ! to set a different termination code depending on which
               ! condition was triggered.  MESA provides 9 customizeable
               ! termination codes, named t_xtra1 .. t_xtra9.  You can
               ! customize the messages that will be printed upon exit by
               ! setting the corresponding termination_code_str value.
               ! termination_code_str(t_xtra1) = 'my termination condition'
      
               ! by default, indicate where (in the code) MESA terminated
               if (extras_check_model == terminate) s% termination_code = t_extras_check_model
            end function extras_check_model
      
      
            integer function how_many_extra_history_columns(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               how_many_extra_history_columns = 2
            end function how_many_extra_history_columns
            
            
            subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
               integer, intent(in) :: id, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
               type (star_info), pointer :: s
               
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               
               names(1) = "M_core"
               names(2) = "R_B"
               vals(1:2) = 0
               if (s% x_logical_ctrl(1)) then 
                  vals(1) = s% xtra(1) / Msun   ! dark core mass in Msun
                  vals(2) = s% xtra(2) / Rsun   ! R_B / Rsun 
               end if
               
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               
               ! note: do NOT add the extras names to history_columns.list
               ! the history_columns.list is only for the built-in history column options.
               ! it must not include the new column names you are adding here.
               
            end subroutine data_for_extra_history_columns
      
      
            integer function how_many_extra_profile_columns(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               how_many_extra_profile_columns = 0
            end function how_many_extra_profile_columns
            
            
            subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
               use star_def, only: star_info, maxlen_profile_column_name
               use const_def, only: dp
               integer, intent(in) :: id, n, nz
               character (len=maxlen_profile_column_name) :: names(n)
               real(dp) :: vals(nz,n)
               integer, intent(out) :: ierr
               type (star_info), pointer :: s
               integer :: k, op_err, net_lwork
               logical :: okay
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
          end subroutine data_for_extra_profile_columns
      
      
            integer function how_many_extra_history_header_items(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               how_many_extra_history_header_items = 0
            end function how_many_extra_history_header_items
      
      
            subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
               integer, intent(in) :: id, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               type(star_info), pointer :: s
               integer, intent(out) :: ierr
               ierr = 0
               call star_ptr(id,s,ierr)
               if(ierr/=0) return
      
               ! here is an example for adding an extra history header item
               ! also set how_many_extra_history_header_items
               ! names(1) = 'mixing_length_alpha'
               ! vals(1) = s% mixing_length_alpha
      
            end subroutine data_for_extra_history_header_items
      
      
            integer function how_many_extra_profile_header_items(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               how_many_extra_profile_header_items = 0
            end function how_many_extra_profile_header_items
      
      
            subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
               integer, intent(in) :: id, n
               character (len=maxlen_profile_column_name) :: names(n)
               real(dp) :: vals(n)
               type(star_info), pointer :: s
               integer, intent(out) :: ierr
               ierr = 0
               call star_ptr(id,s,ierr)
               if(ierr/=0) return
      
               ! here is an example for adding an extra profile header item
               ! also set how_many_extra_profile_header_items
               ! names(1) = 'mixing_length_alpha'
               ! vals(1) = s% mixing_length_alpha
      
            end subroutine data_for_extra_profile_header_items
      
      
            ! returns either keep_going or terminate.
            ! note: cannot request retry; extras_check_model can do that.
            integer function extras_finish_step(id)
               integer, intent(in) :: id
               integer :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               extras_finish_step = keep_going
      
               ! to save a profile, 
                  ! s% need_to_save_profiles_now = .true.
               ! to update the star log,
                  ! s% need_to_update_history_now = .true.
      
               ! see extras_check_model for information about custom termination codes
               ! by default, indicate where (in the code) MESA terminated
               if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
            end function extras_finish_step
            
            
            subroutine extras_after_evolve(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
               type (star_info), pointer :: s
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
            end subroutine extras_after_evolve

            subroutine other_cgrav(id, ierr)
               use star_def
               integer, intent(in) :: id
               integer, intent(out) :: ierr
               type (star_info), pointer :: s
               integer :: k
               ierr = 0
               call star_ptr(id, s, ierr)
               if (ierr /= 0) return
               s% cgrav(:) = 6.674484d-8
            end subroutine other_cgrav
      
            end module run_star_extras
