!*****************************************************************************************
subroutine cal_ann_cent2(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_poisson):: poisson
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
    real(8):: dipole(3) 
    real(8):: cte=0.529177210d0
    real(8):: cte_2=1.d0/0.529177210d0
    character(6) :: filename
    call f_routine(id='cal_ann_cent2')
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%chi_i(1:atoms%nat))
        allocate(ann_arr%chi_1(1:atoms%nat))
        allocate(ann_arr%chi_2(1:atoms%nat))
        allocate(ann_arr%chi_d(1:atoms%nat))
        allocate(ann_arr%a(1:(2*atoms%nat+1)*(2*atoms%nat+1)))
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_1=0.d0
        ann_arr%chi_2=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_1=0.d0
        ann_arr%chi_2=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    endif
    allocate(ann_arr%gama(1:2*atoms%nat))
    ann_arr%gama(1:2*atoms%nat)=0.d0
    if(parini%iverbose>=2) call cpu_time(time1)
    call init_electrostatic_cent2(parini,atoms,ann_arr,ann_arr%a,poisson)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr_tmp)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fatpq(1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%stresspq(1:3,1:3,1:symfunc%linked_lists%maxbound_rad))
    endif
    if(parini%iverbose>=2) call cpu_time(time3)
    over_iat: do iat=1,atoms%nat
        i=atoms%itypat(iat)
        ng=ann_arr%ann(i)%nn(0)
        ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat)
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),out_ann)
            call cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
        elseif(trim(ann_arr%event)=='train') then
            call cal_architecture_der(ann_arr%ann(i),out_ann)
            ann_arr%chi_i(iat)=out_ann
            tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
            ann_arr%chi_1(iat)=ann_arr%ann(i)%chi1
            ann_arr%chi_2(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ann_arr%num(1),iat)=ann_arr%g_per_atom(1:ann_arr%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_gama_cent2(atoms,ann_arr) 
    call get_qat_from_chi_cent2(parini,ann_arr,atoms,poisson,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    !--------------------------------------------------------------------------
    if(parini%iverbose>=2 .and. trim(ann_arr%event)=='evalu') then
        write(filename,'(a4,I2.2)')'dpm_',ann_arr%istep_opt_ann
        open(unit=7820,file=trim(filename),access='append')
        !write(*,*) ann_arr%istep_opt_ann,atoms%ratp(1,1)
        dipole(1)=0.d0 ; dipole(2)=0.d0 ; dipole(3)=0.d0
        call update_ratp(atoms)
        write(7820,*) 'conf:'
        write(7820,*)'nat: ', atoms%nat
        do iat=1,atoms%nat
            dipole(1)=dipole(1)+(atoms%qat_1(iat)+atoms%qat_2(iat))*atoms%ratp(1,iat)
            dipole(2)=dipole(2)+(atoms%qat_1(iat)+atoms%qat_2(iat))*atoms%ratp(2,iat)
            dipole(3)=dipole(3)+(atoms%qat_1(iat)+atoms%qat_2(iat))*atoms%ratp(3,iat)
            write(7820,'(a3,3es16.8)') atoms%sat(iat),cte*atoms%ratp(1,iat),cte*atoms%ratp(2,iat),cte*atoms%ratp(3,iat)
        enddo
        write(7820,'(a11,3es16.8)')  'dpm(cent): ',dipole(1),dipole(2),dipole(3) 
        write(7820,'(a11,3es16.8)')  'dpm(org): ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3)
        write(7820,'(a11,3es16.8)')  'dpm(err): ',atoms%dpm(1)-dipole(1),atoms%dpm(2)-dipole(2),atoms%dpm(3)-dipole(3)
        write(7820,'(a11,3es16.8)')  'dpm(rmse):',&
            sqrt((atoms%dpm(1)-dipole(1))**2+(atoms%dpm(2)-dipole(2))**2+(atoms%dpm(3)-dipole(3))**2)
        !call yaml_mapping_open('cep',flow=.true.)
        !call yaml_map('dpx',dipole(1),fmt='(f10.3)')
        !call yaml_map('dpy',dipole(2),fmt='(f10.3)')
        !call yaml_map('dpz',dipole(3),fmt='(f10.3)')
        !call yaml_mapping_close()
        close(7820)
    endif
    !--------------------------------------------------------------------------
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call get_electrostatic_cent2(parini,atoms,ann_arr,epot_c,ann_arr%a,poisson)
!    call get_gama_force_cent2(atoms,ann_arr)
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        call yaml_mapping_open('Timing of CENT2')
        call yaml_map('initialize matrix',time2-time1)
        call yaml_map('calculation of symfunc',time3-time2)
        call yaml_map('neural network process',time4-time3)
        call yaml_map('linear equations solver',time5-time4)
        call yaml_map('force (SR term)',time6-time5)
        call yaml_map('energy (SR+LR), force (LR)',time7-time6)
        call yaml_map('total time',time7-time1)
        call yaml_mapping_close()
    endif !end of if for printing out timing.
    atoms%epot=epot_c
    if(trim(ann_arr%event)=='evalu') then
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        do iat=1,atoms%nat
            fx_es=atoms%fat(1,iat)-ann_arr%fat_chi(1,iat)
            fy_es=atoms%fat(2,iat)-ann_arr%fat_chi(2,iat)
            fz_es=atoms%fat(3,iat)-ann_arr%fat_chi(3,iat)
            tt1=tt1+fx_es**2+fy_es**2+fz_es**2
            tt2=tt2+ann_arr%fat_chi(1,iat)**2+ann_arr%fat_chi(2,iat)**2+ann_arr%fat_chi(3,iat)**2
            tt3=tt3+fx_es*ann_arr%fat_chi(1,iat)+fy_es*ann_arr%fat_chi(2,iat)+fz_es*ann_arr%fat_chi(3,iat)
        enddo
        tt1=sqrt(tt1)
        tt2=sqrt(tt2)
        ann_arr%fchi_angle=tt3/(tt1*tt2)
        ann_arr%fchi_norm=tt2/tt1
    endif
    call fini_electrostatic_cent2(parini,ann_arr,atoms,poisson)
    !call repulsive_potential_cent(parini,atoms,ann_arr)
    call getvol_alborz(atoms%cellvec,vol)
    call invertmat_alborz(atoms%cellvec,hinv)
    !The following line is inconsistent with the definition of stress tensor
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)*vol
    do i=1,3
    do j=1,3
        atoms%celldv(i,j)=vol*(atoms%stress(i,1)*hinv(j,1)+atoms%stress(i,2)*hinv(j,2)+atoms%stress(i,3)*hinv(j,3))
    enddo
    enddo

    deallocate(symfunc%linked_lists%prime_bound)
    deallocate(symfunc%linked_lists%bound_rad)
    deallocate(symfunc%linked_lists%bound_ang)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%chi_i)
        deallocate(ann_arr%chi_1)
        deallocate(ann_arr%chi_2)
        deallocate(ann_arr%chi_d)
        deallocate(ann_arr%a)
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
    endif
    deallocate(ann_arr%gama)
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine get_qat_from_chi_cent2(parini,ann_arr,atoms,poisson,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(inout):: a(2*atoms%nat+1,2*atoms%nat+1)
    !local variables
    integer:: iat
    real(8):: qtot
    character(200):: smsg
    if(trim(ann_arr%syslinsolver)=='direct') then
        if(trim(atoms%boundcond)=='free') then
            call get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,a)
        else
            smsg='ERROR: solving linear system of equations with direct approach is'
            smsg=trim(smsg)//' possible only for free BC, atoms%boundcond= '
            write(*,'(2a)') trim(smsg),trim(atoms%boundcond)
            stop
        endif
    elseif(trim(ann_arr%syslinsolver)=='apply_matrix') then
        if(trim(atoms%boundcond)=='free') then
            call get_qat_from_chi_iter(parini,ann_arr,atoms,a)
        else
            smsg='ERROR: solving linear system of equations with explicitly'
            smsg=trim(smsg)//' applying matrix is possible only for free BC, atoms%boundcond= '
            write(*,'(2a)') trim(smsg),trim(atoms%boundcond)
            stop
        endif
    elseif(trim(ann_arr%syslinsolver)=='operator') then
        call get_qat_from_chi_operator(parini,poisson,ann_arr,atoms)
    else
        stop 'ERROR: unknown syslinsolver'
    endif
    if(parini%iverbose>=2) then
        call yaml_map('charge_1 on atoms',atoms%qat_1(1:atoms%nat),fmt='(f10.5)')
        call yaml_map('charge_2 on atoms',atoms%qat_2(1:atoms%nat),fmt='(f10.5)')
    endif
end subroutine get_qat_from_chi_cent2
!*****************************************************************************************
subroutine get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(2*atoms%nat+1,2*atoms%nat+1)
    !local variables
    integer:: info !, iat
    character(len=40) :: format_string, str_nat
    associate(nat=>atoms%nat)
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%ipiv(1:2*nat+1))
        allocate(ann_arr%qq(1:2*nat+1))
    endif
    call DGETRF(2*nat+1,2*nat+1,a,2*nat+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    ann_arr%qq(1:nat)=-ann_arr%chi_1(1:nat)-ann_arr%gama(1:nat)
    ann_arr%qq(nat+1:2*nat)=-ann_arr%chi_2(1:nat)-ann_arr%gama(nat+1:2*nat)
    ann_arr%qq(2*nat+1)=atoms%qtot
    call DGETRS('N',2*nat+1,1,a,2*nat+1,ann_arr%ipiv,ann_arr%qq,2*nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat_1(1:nat)=ann_arr%qq(1:nat)
    atoms%qat_2(1:nat)=ann_arr%qq(nat+1:2*nat)
    if(parini%iverbose>=3) then
        write(str_nat,'(I2.2)') nat
        format_string = '(a8,'//trim(str_nat)//'f10.4)'
         write(57,trim(format_string)) 'qat_1',atoms%qat_1
         write(70,trim(format_string)) 'qat_2',atoms%qat_2
    endif
    call charge_analysis_cent2(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',ann_arr%qq(2*nat+1))
        !write(*,*) 'Lagrange ',ann_arr%qq(nat+1)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%ipiv)
        deallocate(ann_arr%qq)
    endif
    end associate
end subroutine get_qat_from_chi_dir_cent2
!*****************************************************************************************
subroutine init_electrostatic_cent2(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(2*atoms%nat+1,2*atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    real(8),allocatable :: gausswidth(:)
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    pi=4.d0*atan(1.d0)
    if(trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train') then
        if(trim(atoms%boundcond)=='free' .and. parini%free_bc_direct) then
            ann_arr%syslinsolver='direct'
        else
            ann_arr%syslinsolver=trim(parini%syslinsolver_ann)
        endif
    else
        ann_arr%syslinsolver=trim(parini%syslinsolver_ann)
    endif
    ann_arr%ener_ref=0.d0
    do iat=1,atoms%nat
        ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
    enddo
    if (.not. parini%ewald) then 
        poisson%alpha = max(maxval(ann_arr%ann(:)%gausswidth_1),maxval(ann_arr%ann(:)%gausswidth_2))
    else 
        if (parini%alpha_ewald<0.d0) then
            call getvol_alborz(atoms%cellvec,vol)
            c=2.2d0
            poisson%alpha = 1.d0/(c*sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
            write(*,*)"optimized alpha = ", poisson%alpha
        else
            poisson%alpha=parini%alpha_ewald
        endif
    end if
    if(trim(ann_arr%syslinsolver)=='direct' .or. trim(ann_arr%syslinsolver)=='apply_matrix') then
        if(trim(atoms%boundcond)/='free') then
            write(*,*) 'ERROR: syslinsolver=direct can be used only for free BC.'
        endif
        call get_amat_cent2(atoms,ann_arr,a)
    elseif(trim(ann_arr%syslinsolver)=='operator') then
        if(trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='slab') then
            allocate(gausswidth(atoms%nat))
            gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
            poisson%task_finit="alloc_rho:set_ngp"
            call init_hartree(parini,atoms,poisson,gausswidth)
            deallocate(gausswidth)
        else
            write(*,'(a)',advance='no') 'ERROR: currently syslinsolver=operator only '
            write(*,'(2a)') 'for BC=bulk/slab, but now BC=',trim(trim(atoms%boundcond))
            stop
        endif
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(ann_arr%syslinsolver)
        stop
    endif
    end associate
end subroutine init_electrostatic_cent2
!*****************************************************************************************
subroutine get_amat_cent2(atoms,ann_arr,a)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(2*atoms%nat+1,2*atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, pi
    real(8):: beta_iat_1, beta_iat_2, beta_jat_1, beta_jat_2
    real(8):: gama_1, gama_2, gama_3, gama_4
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        beta_iat_1=ann_arr%ann(atoms%itypat(iat))%gausswidth_1
        beta_iat_2=ann_arr%ann(atoms%itypat(iat))%gausswidth_2
        gama_1=1.d0/sqrt(beta_iat_1**2+beta_iat_1**2)
        gama_2=1.d0/sqrt(beta_iat_1**2+beta_iat_2**2)
        gama_3=1.d0/sqrt(beta_iat_2**2+beta_iat_1**2)
        gama_4=1.d0/sqrt(beta_iat_2**2+beta_iat_2**2)
        a(iat,iat)=gama_1*2.d0/sqrt(pi)+ann_arr%ann(atoms%itypat(iat))%hardness_1
        a(atoms%nat+iat,atoms%nat+iat)=gama_4*2.d0/sqrt(pi)+ann_arr%ann(atoms%itypat(iat))%hardness_2
        a(iat,atoms%nat+iat)=(gama_2+gama_3)/sqrt(pi)
        a(atoms%nat+iat,iat)=(gama_2+gama_3)/sqrt(pi)
        a(iat,2*atoms%nat+1)=1.d0
        a(atoms%nat+iat,2*atoms%nat+1)=1.d0
        a(2*atoms%nat+1,iat)=1.d0
        a(2*atoms%nat+1,atoms%nat+iat)=1.d0
        do jat=iat+1,atoms%nat
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat_1=ann_arr%ann(atoms%itypat(jat))%gausswidth_1
            beta_jat_2=ann_arr%ann(atoms%itypat(jat))%gausswidth_2
            gama_1=1.d0/sqrt(beta_iat_1**2+beta_jat_1**2)
            gama_2=1.d0/sqrt(beta_iat_1**2+beta_jat_2**2)
            gama_3=1.d0/sqrt(beta_iat_2**2+beta_jat_1**2)
            gama_4=1.d0/sqrt(beta_iat_2**2+beta_jat_2**2)
            a(iat,jat)=erf(gama_1*r)/r
            a(jat,iat)=a(iat,jat)
            a(atoms%nat+iat,atoms%nat+jat)=erf(gama_4*r)/r
            a(atoms%nat+jat,atoms%nat+iat)=a(atoms%nat+iat,atoms%nat+jat)
            a(iat,atoms%nat+jat)=erf(gama_2*r)/r
            a(jat,atoms%nat+iat)=a(iat,atoms%nat+jat)
            a(atoms%nat+iat,jat)=erf(gama_3*r)/r
            a(atoms%nat+jat,iat)=a(atoms%nat+iat,jat)
        enddo
    enddo
    a((2*atoms%nat)+1,(2*atoms%nat)+1)=0.d0
end subroutine get_amat_cent2
!*****************************************************************************************
subroutine fini_electrostatic_cent2(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    !use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    if(trim(ann_arr%syslinsolver)=='operator') then
        call fini_hartree(parini,atoms,poisson)
    endif
end subroutine fini_electrostatic_cent2
!*****************************************************************************************
subroutine get_electrostatic_cent2(parini,atoms,ann_arr,epot_c,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(2*atoms%nat+1,2*atoms%nat+1) ! ASK Dr. Ghasemi
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, tt3, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    tt1=0.d0
    tt2=0.d0
    tt3=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_1(iat)*atoms%qat_1(iat)+&
                ann_arr%chi_2(iat)*atoms%qat_2(iat)
        tt2=tt2+atoms%qat_1(iat)**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness_1+&
                atoms%qat_2(iat)**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness_2
        tt3=tt3+atoms%qat_1(iat)*ann_arr%gama(iat)+atoms%qat_2(iat)*ann_arr%gama(iat+atoms%nat) !FAKEPOT
    enddo
    call cal_electrostatic_ann_cent2(parini,atoms,ann_arr,a,poisson)
    epot_c=epot_es+tt1+tt2+tt3+ann_arr%ener_ref
    end associate
end subroutine get_electrostatic_cent2
!*****************************************************************************************
subroutine cal_electrostatic_ann_cent2(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(2*atoms%nat+1,2*atoms%nat+1)! ASK Dr. Ghasemi
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: tt2, tt3, tt4, ttf, pi
    real(8):: gama_1, gama_2, gama_3, gama_4 
    real(8):: beta_iat_1, beta_iat_2, beta_jat_1, beta_jat_2
    real(8):: dx, dy, dz, r, ehartree_t
    real(8), allocatable:: gausswidth(:)
    allocate(gausswidth(1:atoms%nat))
    if(trim(atoms%boundcond)=='free') then
        pi=4.d0*atan(1.d0)
        tt2=0.d0
        tt3=0.d0
        tt4=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            beta_iat_1=ann_arr%ann(atoms%itypat(iat))%gausswidth_1
            beta_iat_2=ann_arr%ann(atoms%itypat(iat))%gausswidth_2
            gama_1=1.d0/sqrt(beta_iat_1**2+beta_iat_1**2)
            gama_2=1.d0/sqrt(beta_iat_1**2+beta_iat_2**2)
            gama_3=1.d0/sqrt(beta_iat_2**2+beta_iat_1**2)
            gama_4=1.d0/sqrt(beta_iat_2**2+beta_iat_2**2)
            tt2=tt2+(atoms%qat_1(iat)**2*gama_1 + atoms%qat_2(iat)**2*gama_4)/sqrt(pi)
            tt4=tt4+(atoms%qat_1(iat)*atoms%qat_2(iat)*(gama_2+gama_3))/sqrt(pi)
            do jat=iat+1,atoms%nat
                dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
                dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
                dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                beta_jat_1=ann_arr%ann(atoms%itypat(jat))%gausswidth_1
                beta_jat_2=ann_arr%ann(atoms%itypat(jat))%gausswidth_2
                gama_1=1.d0/sqrt(beta_iat_1**2+beta_jat_1**2)
                gama_2=1.d0/sqrt(beta_iat_1**2+beta_jat_2**2)
                gama_3=1.d0/sqrt(beta_iat_2**2+beta_jat_1**2)
                gama_4=1.d0/sqrt(beta_iat_2**2+beta_jat_2**2)
                tt3=tt3+( atoms%qat_1(iat)*atoms%qat_1(jat)*erf(gama_1*r) + atoms%qat_1(iat)*atoms%qat_2(jat)*erf(gama_2*r) &
                    + atoms%qat_2(iat)*atoms%qat_1(jat)*erf(gama_3*r) + atoms%qat_2(iat)*atoms%qat_2(jat)*erf(gama_4*r))/r
                !tt3=tt3+atoms%qat(iat)*atoms%qat(jat)*erf(gama*r)/r
                ttf=(2.d0/sqrt(pi)*gama_1*exp(-(gama_1*r)**2)/r**2-erf(gama_1*r)/r**3)*atoms%qat_1(iat)*atoms%qat_1(jat)
                ttf=ttf+(2.d0/sqrt(pi)*gama_2*exp(-(gama_2*r)**2)/r**2-erf(gama_2*r)/r**3)*atoms%qat_1(iat)*atoms%qat_2(jat)
                ttf=ttf+(2.d0/sqrt(pi)*gama_3*exp(-(gama_3*r)**2)/r**2-erf(gama_3*r)/r**3)*atoms%qat_2(iat)*atoms%qat_1(jat)
                ttf=ttf+(2.d0/sqrt(pi)*gama_4*exp(-(gama_4*r)**2)/r**2-erf(gama_4*r)/r**3)*atoms%qat_2(iat)*atoms%qat_2(jat)
                atoms%fat(1,iat)=atoms%fat(1,iat)+ttf*dx
                atoms%fat(2,iat)=atoms%fat(2,iat)+ttf*dy
                atoms%fat(3,iat)=atoms%fat(3,iat)+ttf*dz
                atoms%fat(1,jat)=atoms%fat(1,jat)-ttf*dx
                atoms%fat(2,jat)=atoms%fat(2,jat)-ttf*dy
                atoms%fat(3,jat)=atoms%fat(3,jat)-ttf*dz
            enddo
        enddo
        ann_arr%epot_es=tt2+tt3+tt4
    elseif(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
        call get_hartree(parini,poisson,atoms,gausswidth,ehartree_t)
        poisson%gw(1:poisson%nat)=poisson%gw_ewald(1:poisson%nat)
        call get_hartree_force(parini,poisson,atoms)
    else
        stop 'ERROR: the requested BCs is not yet implemented.'
    endif
    deallocate(gausswidth)
end subroutine cal_electrostatic_ann_cent2
!*****************************************************************************************
subroutine charge_analysis_cent2(parini,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: iat, i 
    real(8):: q_1, q_2, c_1, c_2, q 
    real(8):: chi_1_min_per_conf(10), chi_1_max_per_conf(10)
    real(8):: chi_2_min_per_conf(10), chi_2_max_per_conf(10)
    chi_1_min_per_conf(1:10)= 1.d20
    chi_1_max_per_conf(1:10)=-1.d20
    chi_2_min_per_conf(1:10)= 1.d20
    chi_2_max_per_conf(1:10)=-1.d20
    do iat=1,atoms%nat
        q=atoms%qat_1(iat)+atoms%qat_2(iat)
        q_1=atoms%qat_1(iat)
        q_2=atoms%qat_2(iat)
        c_1=ann_arr%chi_1(iat)
        c_2=ann_arr%chi_2(iat)
        i=atoms%itypat(iat)
        ann_arr%natsum(i)=ann_arr%natsum(i)+1

        ann_arr%q_1_min(i)=min(q_1,ann_arr%q_1_min(i))
        ann_arr%q_1_max(i)=max(q_1,ann_arr%q_1_max(i))
        ann_arr%q_1_sum(i)=ann_arr%q_1_sum(i)+q_1
        ann_arr%q_2_min(i)=min(q_2,ann_arr%q_2_min(i))
        ann_arr%q_2_max(i)=max(q_2,ann_arr%q_2_max(i))
        ann_arr%q_2_sum(i)=ann_arr%q_2_sum(i)+q_2
        ann_arr%qmin(i)=min(q,ann_arr%qmin(i))
        ann_arr%qmax(i)=max(q,ann_arr%qmax(i))
        ann_arr%qsum(i)=ann_arr%qsum(i)+q

        ann_arr%chi_1_min(i)=min(c_1,ann_arr%chi_1_min(i))
        ann_arr%chi_1_max(i)=max(c_1,ann_arr%chi_1_max(i))
        ann_arr%chi_2_min(i)=min(c_2,ann_arr%chi_2_min(i))
        ann_arr%chi_2_max(i)=max(c_2,ann_arr%chi_2_max(i))
        chi_1_min_per_conf(i)=min(c_1,chi_1_min_per_conf(i))
        chi_1_max_per_conf(i)=max(c_1,chi_1_max_per_conf(i))
        chi_2_min_per_conf(i)=min(c_2,chi_2_min_per_conf(i))
        chi_2_max_per_conf(i)=max(c_2,chi_2_max_per_conf(i))
        ann_arr%chi_1_sum(i)=ann_arr%chi_1_sum(i)+c_1
        ann_arr%chi_2_sum(i)=ann_arr%chi_2_sum(i)+c_2
    enddo
    do i=1,ann_arr%nann
        ann_arr%chi_1_delta(i)=max(ann_arr%chi_1_delta(i),chi_1_max_per_conf(i)-chi_1_min_per_conf(i))
        ann_arr%chi_2_delta(i)=max(ann_arr%chi_2_delta(i),chi_2_max_per_conf(i)-chi_2_min_per_conf(i))
    enddo
end subroutine charge_analysis_cent2
!*****************************************************************************************
subroutine get_gama_cent2(atoms,ann_arr)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: iat
    real(8):: dx, dy, dz, rsq, pi
    real(8):: a1, a2, sigma, omega
    real(8):: cft1, cfb1, cft2, cfb2 
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        a1=ann_arr%ann(atoms%itypat(iat))%gausswidth_1
        a2=ann_arr%ann(atoms%itypat(iat))%gausswidth_2
        sigma=atoms%fakegw
        omega=atoms%fakecoeff
        if (iat==atoms%fakeindex) then
            ann_arr%gama(iat)=3.d0*omega*(a1**2)*(sigma**2)/(2.d0*((a1**2+sigma**2)**2.5d0)*(pi**1.5d0))
            ann_arr%gama(atoms%nat+iat)=3.d0*omega*(a2**2)*(sigma**2)/(2.d0*((a2**2+sigma**2)**2.5d0)*(pi**1.5d0))
        else
            dx=atoms%ratp(1,iat)-atoms%ratp(1,atoms%fakeindex)
            dy=atoms%ratp(2,iat)-atoms%ratp(2,atoms%fakeindex)
            dz=atoms%ratp(3,iat)-atoms%ratp(3,atoms%fakeindex)
            rsq=dx*dx+dy*dy+dz*dz
            cft1=omega*sigma*(3.d0*(a1**4)+3.d0*(a1**2)*(sigma**2)+2.d0*(sigma**2)*rsq)
            cfb1=2.d0*a1*sqrt(a1**(-2) + sigma**(-2))*((a1**2+sigma**2)**3)*(pi**1.5)
            cft2=omega*sigma*(3.d0*(a2**4)+3.d0*(a2**2)*(sigma**2)+2.d0*(sigma**2)*rsq)
            cfb2=2.d0*a2*sqrt(a2**(-2) + sigma**(-2))*((a2**2+sigma**2)**3)*(pi**1.5)
            ann_arr%gama(iat)=cft1*exp(-rsq/(a1**2+sigma**2))/cfb1
            ann_arr%gama(atoms%nat+iat)=cft2*exp(-rsq/(a2**2+sigma**2))/cfb2
        endif
    enddo
end subroutine get_gama_cent2
!*****************************************************************************************
subroutine get_gama_force_cent2(atoms,ann_arr)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: iat
    real(8):: dx, dy, dz, rsq, pi
    real(8):: a1, a2, sigma, omega
    real(8):: gama_d1, gama_d2
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    gama_d1=0.d0
    gama_d2=0.d0
    do iat=1,atoms%nat
        a1=ann_arr%ann(atoms%itypat(iat))%gausswidth_1
        a2=ann_arr%ann(atoms%itypat(iat))%gausswidth_2
        sigma=atoms%fakegw
        omega=atoms%fakecoeff
        if (iat==atoms%fakeindex) then
            dx=1.d0;dy=1.d0;dz=1.d0
            gama_d1=-((omega*sigma*(3*a1**4 + a1**2*sigma**2 - 2*sigma**4))&
                /(a1*Pi**1.5*Sqrt(a1**(-2) + sigma**(-2))*(a1**2 + sigma**2)**4))
            gama_d2=-((omega*sigma*(3*a2**4 + a2**2*sigma**2 - 2*sigma**4))&
                /(a2*Pi**1.5*Sqrt(a2**(-2) + sigma**(-2))*(a2**2 + sigma**2)**4))
        else
            dx=atoms%ratp(1,iat)-atoms%ratp(1,atoms%fakeindex)
            dy=atoms%ratp(2,iat)-atoms%ratp(2,atoms%fakeindex)
            dz=atoms%ratp(3,iat)-atoms%ratp(3,atoms%fakeindex)
            rsq=dx*dx+dy*dy+dz*dz
            gama_d1=-((omega*sigma*(3*a1**4 + a1**2*sigma**2 + 2*rsq*sigma**2 - 2*sigma**4))&
                     /(a1*Exp(rsq/(a1**2 + sigma**2))*Pi**1.5*Sqrt(a1**(-2) + sigma**(-2))&
                     *(a1**2 + sigma**2)**4))
            gama_d2=-((omega*sigma*(3*a2**4 + a2**2*sigma**2 + 2*rsq*sigma**2 - 2*sigma**4))&
                    /(a2*Exp(rsq/(a2**2 + sigma**2))*Pi**1.5*Sqrt(a2**(-2) + sigma**(-2))&
                    *(a2**2 + sigma**2)**4))
        endif
!        write(46,*) 'before' , iat , atoms%fat(1,iat), atoms%fat(2,iat), atoms%fat(3,iat)
        atoms%fat(1,iat)=atoms%fat(1,iat)+(gama_d1*atoms%qat_1(iat)+gama_d2*atoms%qat_2(iat))*dx
        atoms%fat(2,iat)=atoms%fat(2,iat)+(gama_d1*atoms%qat_1(iat)+gama_d2*atoms%qat_2(iat))*dy
        atoms%fat(3,iat)=atoms%fat(3,iat)+(gama_d1*atoms%qat_1(iat)+gama_d2*atoms%qat_2(iat))*dz
!        write(46,*) 'after' , iat , atoms%fat(1,iat), atoms%fat(2,iat), atoms%fat(3,iat)
    enddo
end subroutine get_gama_force_cent2
!*****************************************************************************************
