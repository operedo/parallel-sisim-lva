module boostdijkstra
    use iso_c_binding

    private
    public :: dijkstra

    include "boostdijkstra_cdef.f90"

contains ! Implementation of the functions. We just wrap the C function here.
    subroutine dijkstra(n,a,b,c,d,e,f,g)
        implicit none
        integer(C_INT), intent(in) :: n
        integer(C_INT), intent(in) :: a
        integer(C_INT), intent(in) :: b(:)
        integer(C_INT), intent(in) :: c(:)
        real(C_DOUBLE), intent(in) :: d(:)
        integer(C_INT), intent(in) :: e
        integer(C_INT), intent(in) :: f(:)
        real(C_DOUBLE), intent(in) :: g(:,:)

        call dijkstra_c(n,a,b,c,d,e,f,g)
    end subroutine
end module
