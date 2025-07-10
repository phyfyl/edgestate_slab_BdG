subroutine ham_qlayer2qlayerribbon_BdG(k,H00,H01)
     use para
     use Kronecker, only : KronProd     ! Kroneker product

     implicit none

     ! loop index
     integer :: iR, i1,i2,j1,j2

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic,zmin,delta_l

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib,nm

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

     complex(dp) :: ratio,sigmay(2,2),I_np(Np,Np),i00(Np*nslab*Num_wann,Np*nslab*Num_wann),I_mu(Np*nslab*num_wann,Np*nslab*num_wann)

     complex(Dp) :: Hijk(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann),Hij_k(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann)

     complex(Dp), intent(out) :: H00(2*Np*nslab*num_wann, 2*Np*nslab*num_wann),H01(2*Np*nslab*num_wann,2*Np*nslab*num_wann)

     complex(Dp) :: Hamk_temp(num_wann,num_wann),Hamk_d(nslab*num_wann,nslab*num_wann),Hamk_Delta(Np*nslab*num_wann,Np*nslab*num_wann), I_norb(num_wann/2,num_wann/2)
     nm=2*Np*nslab*num_wann
     call eye_mat(nm/2,I_mu)
     Hijk=0.0d0
     Hij_k=0.0d0
     sigmay(1,1)= (0.0d0, 0.0d0);   sigmay(1,2)= (0.0d0,-1.0d0)
     sigmay(2,1)= (0.0d0, 1.0d0);   sigmay(2,2)= (0.0d0, 0.0d0)

     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      !  call rotation(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ia= int(new_ia)
        inew_ib= int(new_ib)
        if (abs(new_ia).le.ijmax)then
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ic
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hijk(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hijk(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)

           kdotr=-k*new_ic
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)
           Hij_k(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij_k(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
        endif

     enddo
     H00=0d0
     H01=0d0
     ! i1,j1 row index
     do i1=1,Np
     do j1=1,nslab
     ! i2,j2 column index
     do i2=1,Np
     do j2=1,nslab
       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         H00((i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             (i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             (i2-1)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             (i2-1)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hijk(i2-i1,j2-j1,1:Num_wann,1:Num_wann)

         H00(nm/2+(i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             nm/2+(i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             nm/2+(i2-1)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             nm/2+(i2-1)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = -1.0d0*conjg(Hij_k(i2-i1,j2-j1,1:Num_wann,1:Num_wann))

       endif
     enddo
     enddo
     enddo
     enddo


     I_np=0d0
     do i1=1,np
        I_np(i1,i1)=1d0
     enddo

     I_norb=0d0
     Hamk_temp=0d0
     Hamk_d=0d0
     Hamk_Delta=0d0
     delta_l=100*Angstrom2atomic
     zmin=minval(Origin_cell%wannier_centers_cart(3, :))
        ! i2 row index
     do i2=1,nslab
       do j1=1,num_wann/2
          I_norb(j1,j1)=&
                Delta_BdG*eV2Hartree*exp(-(Origin_cell%wannier_centers_cart(3,j1)&
                +Origin_cell%Ruc(3)*(i2-1)-zmin)/delta_l)
       enddo
       Hamk_temp=kronProd(zi*sigmay,I_norb)
       Hamk_d((i2-1)*num_wann+1:i2*num_wann,(i2-1)*num_wann+1:i2*num_wann)=Hamk_temp
     enddo
     Hamk_Delta=kronProd(I_np,Hamk_d)


     ! i1,j1 row index
     do i1=1,Np
     do j1=1,nslab
     ! i2,j2 column index
     do i2=Np+1,Np*2
     do j2=1,nslab
       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         H01((i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             (i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             (i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             (i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hijk(i2-i1,j2-j1,1:Num_wann,1:Num_wann)

         H01(nm/2+(i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             nm/2+(i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             nm/2+(i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             nm/2+(i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = -1.0d0*conjg(Hij_k(i2-i1,j2-j1,1:Num_wann,1:Num_wann))

       endif
     enddo
     enddo
     enddo
     enddo

     H00(1:nm/2,nm/2+1:nm)=Hamk_Delta
     H00(nm/2+1:nm,1:nm/2)=transpose(conjg(Hamk_Delta))
     H00(1:nm/2,1:nm/2)=H00(1:nm/2,1:nm/2)-mu_BdG*eV2Hartree*I_mu
     H00(nm/2+1:nm,nm/2+1:nm)=H00(nm/2+1:nm,nm/2+1:nm)+mu_BdG*eV2Hartree*I_mu

     return
  end subroutine ham_qlayer2qlayerribbon_BdG