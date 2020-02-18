!*****************************************************************************************
module mod_symfunc
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    private
    public :: symfunc_deallocate
    !public:: 
    type, public:: typ_symfunc
        integer:: ng=-1
        integer:: nat=-1
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: y0d_bond(:,:)
        real(8), allocatable:: y0d(:,:,:)
        real(8), allocatable:: y0dr(:,:,:)
        type(typ_linked_lists):: linked_lists
    end type typ_symfunc
    type, public:: typ_symfunc_arr
        integer:: nconf=-1
        type(typ_symfunc), allocatable:: symfunc(:)
    end type typ_symfunc_arr
contains
!*****************************************************************************************
subroutine symfunc_deallocate(symfunc)
    use dynamic_memory
    implicit none
    type(typ_symfunc), intent(inout):: symfunc
    if(allocated(symfunc%y)) deallocate(symfunc%y)
    if(allocated(symfunc%y0d_bond)) deallocate(symfunc%y0d_bond)
    if(allocated(symfunc%y0d)) deallocate(symfunc%y0d)
    if(allocated(symfunc%y0dr)) deallocate(symfunc%y0dr)

end subroutine symfunc_deallocate
end module mod_symfunc
!*****************************************************************************************
