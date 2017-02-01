!------------------------------------------------------------------------------------------------------------------
! MODULE: utils
!
!> @author Christopher J. Gust
!
!DESCRIPTION: 
!> Contains subroutines for reading and writing to disk. 
!-----------------------------------------------------------------------------------------------------------------
module utils

implicit none 
public read_array,write_array,loadobservables,read_matrix,write_matrix,kronecker

contains

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: write_array
!> @author Luca Guierrieri

!> @brief  Writes an array to a file.
!> @param[in] xarray Array to write to disk.
!> @param[in] nelem Number of elements of array.
!> @param[in] filename Name of file with array data.
!
!-----------------------------------------------------------------------------------------------------------

subroutine write_array(Xarray,Nelem,filename)
       
implicit none
integer, intent(in) :: nelem
double precision, intent(in) :: xarray(nelem)
character(len=100), intent(in) :: filename

integer :: i

!Open File and write matrix into file
open (unit=4, file=filename, status='replace', action='write')
do i=1,nelem
     write(4,100) xarray(i)
end do

100 format (ES30.16)

close(unit=4)

end subroutine write_array

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: write_array
!> @author Luca Guierrieri

!> @brief  Writes an array to a file.
!> @param[in] xarray Array to write to disk.
!> @param[in] nelem Number of elements of array.
!> @param[in] filename Name of file with array data.
!
!-----------------------------------------------------------------------------------------------------------

SUBROUTINE read_array(xarray,Nelem,filename)

implicit none 
integer, intent(in) :: Nelem
character(len=100), intent(in) :: filename
double precision, intent(out) :: xarray(Nelem)

integer :: i

!Open file for reading
open (unit=4, file=filename, status='OLD', action='READ')
do i=1,nelem
     read(4,100) xarray(i)
end do

close(UNIT=4)

100 format (ES30.16)

end subroutine read_array

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: loadobservables
!> @author Matthew Smith
!> @date 2-29-12

!> @brief Reads a matrix from disk.
!> @param[in] ntimeperiods Columns of matrix, observables.
!> @param[in] nseries Rows of matrix, observables.
!> @param[in] filename Name of file read from disk.
!> @param[out] observables Matrix output read from disk.  
!
!-----------------------------------------------------------------------------------------------------------


subroutine loadobservables(ntimeperiods,nseries,observables,filename)
           
implicit none 

integer, intent(in) :: ntimeperiods
integer, intent(in) :: nseries
character(len=100), intent(in) :: filename
real(8), intent(out) :: observables(nseries,ntimeperiods)

!Local variables
integer :: t

open(unit=4,file=filename,status='OLD', action='READ')

do t = 1,ntimeperiods
    read(4,*) observables(:,t)
end do

close(unit=4)

end subroutine loadobservables

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: read_matrix
!> @author Christopher Gust
!> @date 7-9-15

!> @brief Reads an array from disk and stores as 2 dimensional matrix.
!> @param[in] inputfile Filename with data to be read from disk
!> @param[in] nrow Number of rows of the output matrix.
!> @param[in] ncol Number of columns of the output matrix.
!> @param[out] matrixout Output matrix.
!
!-----------------------------------------------------------------------------------------------------------
subroutine read_matrix(inputfile,nrow,ncol,matrixout)

implicit none 
integer, intent(in) :: nrow,ncol
character(len=100), intent(in) :: inputfile
double precision, intent(out) :: matrixout(nrow,ncol)

double precision :: inputvec(nrow*ncol)

call read_array(inputvec,nrow*ncol,inputfile)
matrixout = reshape(inputvec,(/nrow,ncol/)) 

end subroutine read_matrix

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: write_matrix
!> @author Christopher Gust
!> @date 7-9-15

!> @brief Reads a matrix to disk.
!> @param[in] outputfile Filename for data written to disk.
!> @param[in] nrow Number of rows of the input matrix.
!> @param[in] ncol Number of columns of the input matrix.
!> @param[out] matrixin Matrix to be written to disk.
!
!-----------------------------------------------------------------------------------------------------------
subroutine write_matrix(outputfile,nrow,ncol,matrixin)

implicit none 
integer, intent(in) :: nrow,ncol
character(len=100), intent(in) :: outputfile
double precision, intent(in) :: matrixin(nrow,ncol)

double precision :: outputvec(nrow*ncol)

outputvec = reshape(matrixin,(/nrow*ncol/))  
call write_array(outputvec,nrow*ncol,outputfile)

end subroutine write_matrix

!--------------------------------------------------------------------------------------------------------
! SUBROUTINE: kronecker
!> @author Christopher Gust and Edward Herbst 
!> @date 2-8-16

!> @brief Computes kronecker product: C = alpha*C + beta*A kron B.  Uses Herbst code amended slightly.
!> @param[in] ma Rows of A.
!> @param[in] na Cols of A.
!> @param[in] mb Rows of B.
!> @param[in] nb Cols of B.
!> @param[in] mc Rows of C.
!> @param[in] nc Cols of C.
!> @param[in] alpha Scalar.
!> @param[in] beta Scalar.  
!> @param[in] A Input matrix, maxna.
!> @param[in] B Input matrix, mbxnb.
!> @param     C Input\output matrix, mcxmc.
!-----------------------------------------------------------------------------------------------------------
subroutine kronecker(ma, na, mb, nb, mc, nc, alpha, beta, A, B, C)
 
implicit none
integer, intent(in) :: ma, na, mb, nb, mc, nc
double precision, intent(in) :: alpha, beta
double precision, dimension(ma,na), intent(in) :: A
double precision, dimension(mb,nb), intent(in) :: B
double precision, dimension(mc,nc), intent(inout) :: C

integer :: i, j, i1, j1 
 
do j = 1, na
   do i = 1, ma
      i1 = (i-1)*mb + 1
      j1 = (j-1)*nb + 1
      C(i1:i1+mb-1,j1:j1+nb-1) = alpha*C(i1:i1+mb-1,j1:j1+nb-1) + beta*A(i,j)*B
   end do
end do

end subroutine kronecker

end module utils
