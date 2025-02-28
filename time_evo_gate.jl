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

function Ham_BC_TDVP(N::Int,sites,BC_width::Int,BC_length::Float64,t_total::Float64,BC_lambda::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)  
  Jx = -J
  Jy = -J + Ky
  Jz = -J - Kz
 
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
  end
   H0=MPO(os_Ham,sites)
  
  Hztime=MPO[]
  for i = 1:N
      global os_Hamt= OpSum()
      global os_Hamt += 1,"Sz",i
      push!(Hztime,MPO(os_Hamt,sites))
  end
  pushfirst!(Hztime, H0)
  hzsites =  Function[t -> hz.*(-atan.((n.-BC_width+1.5-(t/t_total*BC_length+1))/(N/BC_lambda))./pi.-atan.((n.-0.5-(t/t_total*BC_length+1))/(N/BC_lambda))./pi) for n=1:N];
  pushfirst!(hzsites, t -> 1)
  
  Ht = TimeDependentSum(hzsites, Hztime)


  return Ht
end


########## Spin Hamiltonian 2D gate ############
function Heisenberg_Ham2D_TDVP(Nx::Int,Ny::Int,sites,J_inter::Float64,t_total::Float64,Jin_sigma::Float64,J_movinglength::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
  N = Nx * Ny

  J_inter = -J_inter
  Jx = -J
  Jy = -J + Ky
  Jz = -J - Kz
  hzsites = [0.0 for n=1:N]
  hzsites[1]=hz
  hzsites[N-1]=-hz
  if BC == 10
     hzsites[2]= 0
     hzsites[N]= 0
  elseif BC == 11
      hzsites[2]= hz
      hzsites[N]= -hz
  elseif BC == 12
      hzsites[2]= -hz
      hzsites[N]= hz
  elseif BC == 13
      hzsites[2]= hz
      hzsites[N]= hz
  end
  ##### Hamiltonian ######


lattice = square_lattice(Nx, Ny; yperiodic=false)



global os_Ham= OpSum()
for i = 1:(Nx-1)
  site_1 =lattice[i].s1
  site_2 =lattice[i].s2

  site_3 =lattice[i+Nx-1].s1
  site_4 =lattice[i+Nx-1].s2
  os_Ham .+= Jx, "Sx", site_1, "Sx", site_2
  os_Ham .+= Jy, "Sy", site_1, "Sy", site_2
  os_Ham .+= Jz, "Sz", site_1, "Sz", site_2

  os_Ham .+= Jx, "Sx", site_3, "Sx", site_4
  os_Ham .+= Jy, "Sy", site_3, "Sy", site_4
  os_Ham .+= Jz, "Sz", site_3, "Sz", site_4
end

for i = 1:N
  os_Ham .+= hy,"Sy",i
  os_Ham .+= hx,"Sx",i
  os_Ham .+= hzsites[i],"Sz",i
end
global    H0=MPO(os_Ham,sites)
HJintime=MPO[]
for i = 1:Nx
    global os_Hamt= OpSum()
    site_1 =lattice[i+2*Nx-2].s1
    site_2 =lattice[i+2*Nx-2].s2

    os_Ham .+= 1, "Sx", site_1, "Sx", site_2
    os_Ham .+= 1, "Sy", site_1, "Sy", site_2
    os_Ham .+= 1, "Sz", site_1, "Sz", site_2
    push!(HJintime,MPO(os_Hamt,sites))
end

  push!(HJintime,MPO(os_Hamt,sites))


  pushfirst!(Hztime, H0)
  Jinsites =  Function[t -> J_inter.*exp(-((n-(-10+t/t_total*J_movinglength)/Nx/Jin_sigma)^2)) for n=1:Nx];
  pushfirst!(Jinsites, t -> 1)
    
  Ht = TimeDependentSum(Jinsites, Hztime)
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