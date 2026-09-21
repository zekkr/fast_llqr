! Experimental multivariate Gaussian provider and deterministic Prim traversal.
subroutine ssqr_gaussian_tree(a,y,x,n,q,d,z,ne,h,tau,tol,maxit,threshold,min_keep,cache_flags, &
                              beta,hseq,diagnostics,visit_order,parent,edge_weight,mst_seconds, &
                              ierr,failed_eval)
  use weighted_qr_tree_core, only: fit_screen_tree
  use iso_fortran_env, only: real64,int64
  implicit none
  integer,intent(in)::n,q,d,ne,maxit,min_keep,cache_flags
  real(real64),intent(in)::a(n,q),y(n),x(n,d),z(ne,d),h(d),tau,tol,threshold(ne)
  real(real64),intent(out)::beta(ne,q),edge_weight(ne),mst_seconds
  integer,intent(out)::hseq(ne,q),diagnostics(ne,18),visit_order(ne),parent(ne),ierr,failed_eval
  real(real64)::scaled_z(ne,d),center(d),nearest2(ne),delta,dist2,tick_rate
  logical::in_tree(ne)
  integer::root,node,candidate,step,i,j,best_parent
  integer(int64)::tick0,tick1,clock_rate

  ierr=0; failed_eval=0; beta=0; hseq=0; diagnostics=0
  visit_order=0; parent=0; edge_weight=0; mst_seconds=0
  if(n<q .or. q/=d+1 .or. d<1 .or. ne<1 .or. any(h<=0)) then
    ierr=5; return
  end if
  do j=1,d
    scaled_z(:,j)=z(:,j)/h(j)
    center(j)=sum(scaled_z(:,j))/real(ne,real64)
  end do

  call system_clock(tick0,clock_rate)
  root=1;dist2=huge(1.0_real64)
  do i=1,ne
    delta=sum((scaled_z(i,:)-center)**2)
    if(delta<dist2) then
      dist2=delta;root=i
    end if
  end do
  in_tree=.false.;nearest2=huge(1.0_real64)
  in_tree(root)=.true.;visit_order(1)=root;parent(root)=0;edge_weight(root)=0
  do i=1,ne
    if(i==root) cycle
    nearest2(i)=sum((scaled_z(i,:)-scaled_z(root,:))**2)
    parent(i)=root
  end do
  do step=2,ne
    candidate=0;dist2=huge(1.0_real64)
    do i=1,ne
      if(.not.in_tree(i) .and. nearest2(i)<dist2) then
        dist2=nearest2(i);candidate=i
      end if
    end do
    if(candidate==0) then
      ierr=5; return
    end if
    node=candidate;in_tree(node)=.true.;visit_order(step)=node
    edge_weight(node)=sqrt(max(0.0_real64,nearest2(node)))
    do i=1,ne
      if(in_tree(i)) cycle
      delta=sum((scaled_z(i,:)-scaled_z(node,:))**2)
      if(delta<nearest2(i)) then
        nearest2(i)=delta;parent(i)=node
      end if
    end do
  end do
  call system_clock(tick1)
  if(clock_rate>0_int64) then
    tick_rate=real(clock_rate,real64)
    mst_seconds=real(tick1-tick0,real64)/tick_rate
  end if

  call fit_screen_tree(a,y,n,q,ne,visit_order,parent,tau,tol,maxit,threshold,min_keep,cache_flags,weights, &
                       beta,hseq,diagnostics,ierr,failed_eval)
contains
  subroutine weights(ev,nw,w,ids,na)
    integer,intent(in)::ev,nw
    real(real64),intent(inout)::w(nw)
    integer,intent(out)::ids(nw),na
    integer::row,col
    real(real64)::radius2,u
    na=0
    do row=1,nw
      radius2=0
      do col=1,d
        u=(x(row,col)-z(ev,col))/h(col)
        radius2=radius2+u*u
      end do
      w(row)=exp(-0.5_real64*radius2)
      if(w(row)>0) then
        na=na+1;ids(na)=row
      end if
    end do
  end subroutine weights
end subroutine ssqr_gaussian_tree
