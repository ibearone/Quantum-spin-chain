#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the time    #
# evolution Hanmiltonian.                   #
############################################# 

############### Define evo_gates for TEBD simulation ##############################
function evo_gates_TEBD_LT_1(N::Int,sites,NBC1::Int,NBC2::Int,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)
 Jx = -J
 Jy = -J + Ky
 Jz = -J - Kz
 hzsites = [0 for n=1:N]
 hzsites[1:NBC1].=hz
 hzsites[NBC2:N].=-hz

 global evo_gates = ITensor[]
  for j in 1:(N - 1)
    s1 = sites[j]
    s2 = sites[j + 1]
    hj_inner =
    hy*op("Sy", s1) +
    hx*op("Sx", s1) +
    hzsites[j]*op("Sz", s1) 
    Gj = exp(-im * tau/2 * hj_inner)
    push!(evo_gates, Gj)
    
    hj_inter =
      Jz*op("Sz", s1) * op("Sz", s2) +
      Jy*op("Sy", s1) * op("Sy", s2) +
      Jx*op("Sx", s1) * op("Sx", s2) 

    Gj = exp(-im * tau/2 * hj_inter)
    push!(evo_gates, Gj)

  end
  for j = N
    s1 = sites[j]
    hj_inner =
      hy*op("Sy", s1) +
      hx*op("Sx", s1) +
      hzsites[j]*op("Sz", s1) 
    Gj = exp(-im * tau/2 * hj_inner)
    push!(evo_gates, Gj)
  end
  append!(evo_gates, reverse(evo_gates));


 return evo_gates
end

function evo_gates_TEBD_LT_1_Ht(N::Int,sites,NBC1::Int,NBC2::Int,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,dhx::Float64,dhy::Float64,omega::Float64,t::Float64,tau::Float64)
  Jx = -J
  Jy = -J + Ky
  Jz = -J - Kz
  hzsites = [0 for n=1:N]
  hzsites[1:NBC1].=hz
  hzsites[NBC2:N].=-hz
 
  global evo_gates = ITensor[]
   for j in 1:(N - 1)
     s1 = sites[j]
     s2 = sites[j + 1]
     hj_inner =
     (hy+dhy*sin(omega * t))*op("Sy", s1) +
     (hx+dhx*cos(omega * t))*op("Sx", s1) +
     hzsites[j]*op("Sz", s1) 
     Gj = exp(-im * tau/2 * hj_inner)
     push!(evo_gates, Gj)
     
     hj_inter =
       Jz*op("Sz", s1) * op("Sz", s2) +
       Jy*op("Sy", s1) * op("Sy", s2) +
       Jx*op("Sx", s1) * op("Sx", s2) 
 
     Gj = exp(-im * tau/2 * hj_inter)
     push!(evo_gates, Gj)
 
   end
   for j = N
     s1 = sites[j]
     hj_inner =
     (hy+dhy*cos(omega * t))*op("Sy", s1) +
     (hx+dhx*sin(omega * t))*op("Sx", s1) +
       hzsites[j]*op("Sz", s1) 
     Gj = exp(-im * tau/2 * hj_inner)
     push!(evo_gates, Gj)
   end
   append!(evo_gates, reverse(evo_gates));
 
 
  return evo_gates
end



function evo_gates_TEBD_LT_5_Ht(N::Int,sites,hx::Float64,hy::Float64,hz::Float64,dhx::Float64,dhy::Float64,omega::Float64,t::Float64,tau::Float64)

  global evo_gates = ITensor[]
  for j in 1:N
    s_i = sites[j]
    hj_inner =
    (hy+dhy*sin(omega * t))*op("Sy", s_i) +
    (hx+dhx*cos(omega * t))*op("Sx", s_i) +
    hz*j*op("Sz", s_i) 
    Gj = exp(-im * tau/2 * hj_inner)
    push!(evo_gates, Gj)
    
    
  end
  append!(evo_gates, reverse(evo_gates))
  return evo_gates
 end

############### Define time dependent part of Hamiltonian in TDVP simulation ##############################
function Ham_time_dependent_TDVP(N::Int,sites,dhx::Float64,dhy::Float64)
  global os_Ham= OpSum()
  for i = 1:N
    global os_Ham += dhy,"Sy",i
    global os_Ham += dhx,"Sx",i
  end
  global H=MPO(os_Ham,sites)
  return H
end

############### Define total Hamiltonian in TDVP simulation ##############################
function Ham_tot_TDVP(N::Int,sites,H_evo::MPO,dhx::Float64,dhy::Float64,omega::Float64,t0::Float64)
  f0 = map(ω -> (t -> 1), 0)
  f1cos = map(ω -> (t -> cos(ω * (t+t0))), omega)
  f1sin = map(ω -> (t -> sin(ω * (t+t0))), omega)
  
  f = (f0,f1cos,f1sin)
  H1cos = Ham_time_dependent_TDVP(N,sites,dhx,0.0)
  H1sin = Ham_time_dependent_TDVP(N,sites,0.0,dhy)
  H = (H_evo , H1cos, H1sin)
  Ht = TimeDependentSum(f, H)

  return Ht
end

function Ham_BC_TDVP(N::Int,sites,BC_width::Int,t_total::Float64,BC_lambda::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)
  d= = map(v -> (t -> t/t_total*(N-BC_width)+1), velocity)
  
  Jx = -J
  Jy = -J + Ky
  Jz = -J - Kz
  hzsites = [hz.*(-atan.((n.-BC_width+1.5-d)/(N/BC_lambda))./pi.-atan.((n.-0.5-d)/(N/BC_lambda))./pi) for n=1:N];
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
  global Ht=MPO(os_Ham,sites)


  return Ht
end

############### Define auto detection of length of psi_evo ##############################
function detect_psi_evo_length(file)
  count = 0
  for key in keys(file)
      if occursin(r"^psi_evo_\d+$", key)  
          count += 1
      end
  end
  return count
end