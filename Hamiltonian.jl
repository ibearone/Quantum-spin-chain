#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the spin    #
# ladder Hanmiltonian wirh normal order.    #
############################################# 

########## Spin Hamiltonian  ############
function Heisenberg_Ham(N::Int,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)
    sites = siteinds("S=1/2",N)
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0 for n=1:N]
    hzsites[1]=hz
    hzsites[N]=-hz
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
return H,sites
end

########## Spin Hamiltonian  ############
function Heisenberg_Ham2D(Nx::Int,Ny::Int,J_inter::Float64,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64,BC)
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
    sites = siteinds("S=1/2",N)
    ##### Hamiltonian ######


 lattice = square_lattice(Nx, Ny; yperiodic=false)



 global os_Ham= OpSum()
 for i = 1:(Nx-1)
    site_1 =lattice[3*i-2].s1
    site_2 =lattice[3*i-2].s2

    site_3 =lattice[3*i].s1
    site_4 =lattice[3*i].s2
    os_Ham .+= Jx, "Sx", site_1, "Sx", site_2
    os_Ham .+= Jy, "Sy", site_1, "Sy", site_2
    os_Ham .+= Jz, "Sz", site_1, "Sz", site_2

    os_Ham .+= Jx, "Sx", site_3, "Sx", site_4
    os_Ham .+= Jy, "Sy", site_3, "Sy", site_4
    os_Ham .+= Jz, "Sz", site_3, "Sz", site_4
 end

 for i = 1:(Nx-1)
    site_1 =lattice[3*i-1].s1
    site_2 =lattice[3*i-1].s2

    os_Ham .+= J_inter, "Sx", site_1, "Sx", site_2
    os_Ham .+= J_inter, "Sy", site_1, "Sy", site_2
    os_Ham .+= J_inter, "Sz", site_1, "Sz", site_2

 end

 os_Ham .+= J_inter, "Sx", lattice[end].s1, "Sx", lattice[end].s2
 os_Ham .+= J_inter, "Sy", lattice[end].s1, "Sy", lattice[end].s2
 os_Ham .+= J_inter, "Sz", lattice[end].s1, "Sz", lattice[end].s2

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


########## Spin Hamiltonian mobile ############
function Heisenberg_Ham_mobile(N::Int,NBC1::Int,NBC2::Int,J::Float64,Kz::Float64,Ky::Float64,hx::Float64,hy::Float64,hz::Float64)
    sites = siteinds("S=1/2",N)
    Jx = -J
    Jy = -J + Ky
    Jz = -J - Kz
    hzsites = [0 for n=1:N]
    hzsites[1:NBC1].=hz
    hzsites[NBC2:N].=-hz
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
return H,sites
end