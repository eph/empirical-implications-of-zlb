module class_RandomNumber

  use mkl_vsl
  use mkl_vsl_type

  include 'mpif.h'
  type RandomNumber
     
     integer :: seed = 1848

     integer :: brng = vsl_brng_mt19937

     integer :: methodu = VSL_METHOD_DUNIFORM_STD
     integer :: methodn = VSL_METHOD_DGAUSSIAN_BOXMULLER


     double precision :: normal_mean = 0.0d0
     double precision :: normal_std = 1.0d0

     double precision :: uniform_lb = 0.0d0
     double precision :: uniform_ub = 1.0d0 

     type(vsl_stream_state) :: stream 

     contains
       
       procedure :: norm_rvs 
       procedure :: uniform_rvs 

  end type RandomNumber

  interface RandomNumber
    module procedure new_RandomNumber
  end interface


  contains
    
    type(RandomNumber) function new_RandomNumber(seed) result(rn)

      integer, intent(in), optional :: seed

      integer :: errcode

      if (present(seed)) rn%seed = seed

      errcode = vslnewstream(rn%stream, rn%brng, rn%seed)

    end function new_RandomNumber

    function norm_rvs(rn, dim_a, dim_b, mu, sig) result(rvs)
      class(RandomNumber) :: rn
      integer, intent(in) :: dim_a, dim_b

      double precision :: rvs(dim_a, dim_b)

      double precision, intent(in), optional :: mu, sig
      double precision :: rvs_mu, rvs_sig
      
      integer :: errcode 

      rvs_mu = rn%normal_mean
      rvs_sig = rn%normal_std

      if (present(mu)) rvs_mu = mu
      if (present(sig)) rvs_sig = sig

      errcode = vdrnggaussian( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_mu, rvs_sig)

    end function norm_rvs

    function uniform_rvs(rn, dim_a, dim_b, lb, ub) result(rvs)
      class(RandomNumber) :: rn
      integer, intent(in) :: dim_a, dim_b

      double precision :: rvs(dim_a, dim_b)

      double precision, intent(in), optional :: lb, ub
      double precision :: rvs_lb, rvs_ub
      
      integer :: errcode 

      rvs_lb = rn%uniform_lb
      rvs_ub = rn%uniform_ub

      if (present(lb)) rvs_lb = lb
      if (present(ub)) rvs_ub = ub

      errcode = vdrnguniform( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_lb, rvs_ub)
    
    end function uniform_rvs
    

end module class_RandomNumber
