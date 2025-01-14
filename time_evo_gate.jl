#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the time    #
# evolution Hanmiltonian.                   #
############################################# 

function Ham_mobile_gates(N::Int,sites,NBC1::Int,NBC2::Int,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)
 Jx = -J
 Jy = -J + Ky
 Jz = -J - Kz
 hzsites = [0 for n=1:N]
 hzsites[1:NBC1].=hz
 hzsites[NBC2:N].=-hz

 global gates = ITensor[]
  for j in 1:(N - 1)
    s1 = sites[j]
    s2 = sites[j + 1]
    hj_inner =
    hy*op("Sy", s1) +
    hx*op("Sx", s1) +
    hzsites[j]*op("Sz", s1) 
    Gj = exp(-im * tau/2 * hj_inner)
    push!(gates, Gj)
    
    hj_inter =
      Jz*op("Sz", s1) * op("Sz", s2) +
      Jy*op("Sy", s1) * op("Sy", s2) +
      Jx*op("Sx", s1) * op("Sx", s2) 

    Gj = exp(-im * tau/2 * hj_inter)
    push!(gates, Gj)

  end
  for j = N
    s1 = sites[j]
    hj_inner =
      hy*op("Sy", s1) +
      hx*op("Sx", s1) +
      hzsites[j]*op("Sz", s1) 
    Gj = exp(-im * tau/2 * hj_inner)
    push!(gates, Gj)
  end
 append!(gates, reverse(gates));

     ##### Hamiltonian ######
     global os_Ham= OpSum()
     for i =1:N-1
         global os_Ham += Jz,"Sz",i,"Sz",i+1
         global os_Ham += Jy,"Sy",i,"Sy",i+1
         global os_Ham += Jx,"Sx",i,"Sx",i+1
     end
     for i = 1:N
         global os_Ham += hy,"Sy",i
         global os_Ham += hx,"Sx",i
         global os_Ham += hzsites[i],"Sz",i
     end
     global H=MPO(os_Ham,sites)

 return H,gates
end

function Ham_time_dependent_gates(N::Int,sites,dhx::Float64,dhy::Float64)
  global os_Ham= OpSum()
  for i = 1:N
    global os_Ham += dhy,"Sy",i
    global os_Ham += dhx,"Sx",i
  end
  global H=MPO(os_Ham,sites)
  return H
end

function Ham_time_tot(N::Int,site,H_time::MPO,dhx::Float64,dhy::Float64,omega::Float64,t0::Float64)
  f0=map(ω -> (t -> 1), 0)
  f1cos=map(ω -> (t -> cos(ω * (t+t0))), omega)
  f1sin=map(ω -> (t -> sin(ω * (t+t0))), omega)
  
  f = (f0,f1cos,f1sin)
  H1cos=Ham_time_dependent_gates(N,sites,dhx,0.0)
  H1sin=Ham_time_dependent_gates(N,sites,0.0,dhy)
  H = (H_time , H1cos, H1sin)
  Ht=TimeDependentSum(f, H)

  return Ht
end