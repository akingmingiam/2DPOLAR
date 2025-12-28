program polar_sim
    use mesh_types,        only : PolarMesh, TimeControlParameters
    use flow_fields,       only : PrimitiveVariables, ConservedVariables
    use fluid_properties,  only : StiffenedGas
    use flux_types,        only : FluxFields, allocate_flux_fields
    use mesh_ops,          only : mesh_setup
    use initial_ops,       only : initialize_flow_fields, initialize_sedov
    use output_writer,     only : write_primitive_tecplot, print_time_info, write_run_info
    use time_step_control, only : compute_global_timestep
    use flux_ops,          only : compute_all_fluxes
    use field_update_ops,  only : update_conserved_field
    use field_conversion,  only : update_primitive_from_conserved
    use boundary_update,   only : update_boundary_primitive_conserved

    implicit none

    type(PolarMesh)             :: mesh
    type(TimeControlParameters) :: time
    type(PrimitiveVariables)    :: prim
    type(ConservedVariables)    :: cons
    type(FluxFields)            :: flux
    type(StiffenedGas)          :: gas


    !-----------------------------------------------------------
    ! Construct mesh + allocate flow field arrays + initialized
    ! allocate flux arrays + Create result folder + Output initial fields
    !-----------------------------------------------------------
    call mesh_setup(mesh, time, prim, cons)
    call initialize_flow_fields(mesh, prim, cons, gas)
    call allocate_flux_fields(flux, mesh)
    call execute_command_line("mkdir -p result")
    call write_run_info("result", mesh%params, time)
    call write_primitive_tecplot(time%t, time%step, time%output_step, "result/", mesh, prim)
    time%output_step = time%output_step + 1
    
    !--------------------------------------------------------------
    ! Main Time Loop
    !--------------------------------------------------------------
    do while (time%t < time%t_end .and. time%step < time%max_steps)

        time%dt = compute_global_timestep(mesh, prim, gas, time%CFL, time%dt_max, time%dt_min)
        time%step = time%step + 1
        call print_time_info(time%step, time%t, time%dt)

        call compute_all_fluxes(mesh, prim, gas, flux)

        call update_conserved_field(mesh, cons, prim, flux, time%dt)

        call update_primitive_from_conserved(prim, cons, gas, &
                                             mesh%params%i_lo_phys, mesh%params%i_hi_phys, &
                                             mesh%params%j_lo_phys, mesh%params%j_hi_phys)

        call update_boundary_primitive_conserved(mesh, prim, cons, gas)
        
        time%t = time%t + time%dt

        if ( (time%t >= time%time_outputs(time%output_step) .and. &
              time%output_step <= time%max_output_steps) ) then
            call write_primitive_tecplot(time%t, time%step, time%output_step, "result/", mesh, prim)
            time%output_step = time%output_step + 1

            if (time%output_step > time%max_output_steps) then
                stop
            end if
        end if

 
        ! if (mod(time%step, 20) == 0) then
        !     call write_primitive_tecplot(time%t, time%step, time%output_step, "result/", mesh, prim)
        !     time%output_step = time%output_step + 1
        ! end if

    end do
end program polar_sim

