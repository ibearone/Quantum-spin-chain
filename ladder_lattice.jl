#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the spin    #
# ladder Hanmiltonian with ladder order.    #
############################################# 

########## Spin Hamiltonian  ############
struct LatticeBond
    s1::Int
    s2::Int
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    type::String
end

function LatticeBond(s1::Int, s2::Int)
    return LatticeBond(s1, s2, 0.0, 0.0, 0.0, 0.0, "")
end
  
function LatticeBond(
    s1::Int, s2::Int, x1::Real, y1::Real, x2::Real, y2::Real, bondtype::String=""
  )
    cf(x) = convert(Float64, x)
    return LatticeBond(s1, s2, cf(x1), cf(y1), cf(x2), cf(y2), bondtype)
end
  
const Lattice = Vector{LatticeBond}


function ladder_lattice(Nx::Int, Ny::Int)::Lattice
   
    N = Nx * Ny
    Nbond = 2N - Ny - Nx
    latt = Lattice(undef, Nbond)
    b = 0
    for n in 1:N
      x = mod(n - 1, Nx) + 1
      y = div(n - 1, Nx) + 1
      if x < Nx
        latt[b += 1] = LatticeBond(n, n + 1, x, y, x + 1, y)
      end

    end
    
    for n in 1:N
      x = mod(n - 1, Nx) + 1
      y = div(n - 1, Nx) + 1
      if y >1
       latt[b += 1] = LatticeBond(n-Nx,n , x, y, x , y - 1)
      end
    end
    return latt
end



function Heisenberg_Ham_ladder(Nx::Int,Ny::Int,J_inter::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
    N = Nx * Ny
  
    J_inter = -J_inter
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0.0 for n=1:N]
    hzsites[1]=hz
    hzsites[Nx]=-hz
    if BC == 10
       hzsites[Nx+1]= 0
       hzsites[N]= 0
    elseif BC == 11
        hzsites[Nx+1]= hz
        hzsites[N]= -hz
    elseif BC == 12
        hzsites[Nx+1]= -hz
        hzsites[N]= hz
    elseif BC == 13
        hzsites[Nx+1]= hz
        hzsites[N]= hz
    end
    sites = siteinds("S=1/2",N)
    ##### Hamiltonian ######


 lattice = ladder_lattice(Nx, Ny)



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

 for i = 1:Nx
    site_1 =lattice[i+2*Nx-2].s1
    site_2 =lattice[i+2*Nx-2].s2

    os_Ham .+= J_inter, "Sx", site_1, "Sx", site_2
    os_Ham .+= J_inter, "Sy", site_1, "Sy", site_2
    os_Ham .+= J_inter, "Sz", site_1, "Sz", site_2

 end



 for i = 1:N
    os_Ham .+= hy,"Sy",i
    os_Ham .+= hx,"Sx",i
    os_Ham .+= hzsites[i],"Sz",i
 end
 global H=MPO(os_Ham,sites)
return H,sites
end


function interleave_arrays(arr1, arr2)
    len = min(length(arr1), length(arr2))
    result = []

    for i in 1:len
        push!(result, arr1[i])
        push!(result, arr2[i])
    end

    # 处理剩余的元素
    append!(result, arr1[len+1:end])
    append!(result, arr2[len+1:end])

    return result
end


function Heisenberg_Ham_ladder_single(Nx::Int,Ny::Int,Nc::Int,J_inter::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
    N = Nx * Ny
  
    J_inter = -J_inter
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0.0 for n=1:N]
    hzsites[1]=hz
    hzsites[Nx]=-hz
    if BC == 10
       hzsites[Nx+1]= 0
       hzsites[N]= 0
    elseif BC == 11
        hzsites[Nx+1]= hz
        hzsites[N]= -hz
    elseif BC == 12
        hzsites[Nx+1]= -hz
        hzsites[N]= hz
    elseif BC == 13
        hzsites[Nx+1]= hz
        hzsites[N]= hz
    end
    sites = siteinds("S=1/2",N)
    ##### Hamiltonian ######


 lattice = ladder_lattice(Nx, Ny)



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

 for i = Nc
    site_1 =lattice[i+2*Nx-2].s1
    site_2 =lattice[i+2*Nx-2].s2

    os_Ham .+= J_inter, "Sx", site_1, "Sx", site_2
    os_Ham .+= J_inter, "Sy", site_1, "Sy", site_2
    os_Ham .+= J_inter, "Sz", site_1, "Sz", site_2

 end



 for i = 1:N
    os_Ham .+= hy,"Sy",i
    os_Ham .+= hx,"Sx",i
    os_Ham .+= hzsites[i],"Sz",i
 end
 global H=MPO(os_Ham,sites)
return H,sites
end

function Heisenberg_Ham_ladder_single_2(sites,Nx::Int,Ny::Int,Nc::Int,J_inter::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
    N = Nx * Ny
  
    J_inter = -J_inter
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0.0 for n=1:N]
    hzsites[1]=hz
    hzsites[Nx]=-hz
    if BC == 10
       hzsites[Nx+1]= 0
       hzsites[N]= 0
    elseif BC == 11
        hzsites[Nx+1]= hz
        hzsites[N]= -hz
    elseif BC == 12
        hzsites[Nx+1]= -hz
        hzsites[N]= hz
    elseif BC == 13
        hzsites[Nx+1]= hz
        hzsites[N]= hz
    end
    ##### Hamiltonian ######


 lattice = ladder_lattice(Nx, Ny)



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

 for i = Nc
    site_1 =lattice[i+2*Nx-2].s1
    site_2 =lattice[i+2*Nx-2].s2

    os_Ham .+= J_inter, "Sx", site_1, "Sx", site_2
    os_Ham .+= J_inter, "Sy", site_1, "Sy", site_2
    os_Ham .+= J_inter, "Sz", site_1, "Sz", site_2

 end



 for i = 1:N
    os_Ham .+= hy,"Sy",i
    os_Ham .+= hx,"Sx",i
    os_Ham .+= hzsites[i],"Sz",i
 end
 global H=MPO(os_Ham,sites)
return H
end



function Heisenberg_Ham_ladder_single_mobile(Nx::Int,Ny::Int,Nc::Int,NBC1C1::Int,NBC2C1::Int,NBC1C2::Int,NBC2C2::Int,J_inter::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
    N = Nx * Ny
  
    J_inter = -J_inter
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0.0 for n=1:N]
    hzsites[1:NBC1C1].=hz
    hzsites[NBC2C1:Nx].=-hz
    if BC == 10
       hzsites[Nx+1]= 0
       hzsites[N]= 0
    elseif BC == 11
        hzsites[Nx+1:Nx+NBC1C2].= hz
        hzsites[Nx+NBC2C2:N].= -hz
    elseif BC == 12
        hzsites[Nx+1:Nx+NBC1C2].= -hz
        hzsites[Nx+NBC2C2:N].= hz
    elseif BC == 13
        hzsites[Nx+1]= hz
        hzsites[N]= hz
    end
    sites = siteinds("S=1/2",N)
    ##### Hamiltonian ######


 lattice = ladder_lattice(Nx, Ny)



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

 for i = Nc
    site_1 =lattice[i+2*Nx-2].s1
    site_2 =lattice[i+2*Nx-2].s2

    os_Ham .+= J_inter, "Sx", site_1, "Sx", site_2
    os_Ham .+= J_inter, "Sy", site_1, "Sy", site_2
    os_Ham .+= J_inter, "Sz", site_1, "Sz", site_2

 end



 for i = 1:N
    os_Ham .+= hy,"Sy",i
    os_Ham .+= hx,"Sx",i
    os_Ham .+= hzsites[i],"Sz",i
 end
 global H=MPO(os_Ham,sites)
return H,sites
end