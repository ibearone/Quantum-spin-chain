#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the main code for spin chain      #
# simulation using ITensor library.         #
#                                           #
############################################# 
########## Initialization ############
using ITensors
using ITensorMPS
using HDF5
using LaTeXStrings
#using Plots
using JLD2
using DelimitedFiles
using Dates

########### Modules ##############
include("read_input.jl")

include("Hamiltonian.jl")
include("ladder_lattice.jl")
include("Spin_Operator.jl")
include("DW_chirality_operators.jl")

include("DMRGloop.jl")

########## Reading Inputs (Hamiltonian) ############

file_in = open("Input", "r")

Lattice_type = read_input(file_in,"Lattice",Int,0)
Mobile_DW = read_input(file_in,"Mobile_DW",Int,0)
if Lattice_type ==1 || Lattice_type ==5
    N  = read_input(file_in,"N",Int,0)
    J  = read_input(file_in,"J",Float64,0)
    if Mobile_DW == 1
        NBC  = read_input(file_in,"NBC",Int,1)
    end

 elseif Lattice_type == 2 || Lattice_type ==3
    Nx  = read_input(file_in,"Nx",Int,0)
    Ny  = read_input(file_in,"Ny",Int,0)
    J_inter= read_input(file_in,"J_in",Float64,0)
    J  = read_input(file_in,"J_ex",Float64,0)
    BC_2D = read_input(file_in,"Bondary Condition",Int,0)
 elseif Lattice_type == 4 
    Nx  = read_input(file_in,"Nx",Int,0)
    Ny  = read_input(file_in,"Ny",Int,0)
    Nc   = read_input(file_in,"Nc ",Int,0)
    J_inter= read_input(file_in,"J_in",Float64,0)
    J  = read_input(file_in,"J_ex",Float64,0)
    BC_2D = read_input(file_in,"Bondary Condition",Int,0)
    if Mobile_DW == 1
        NBC_1  = read_input(file_in,"NBC_Chain1",Int,1)
        NBC_2  = read_input(file_in,"NBC_Chain2",Int,1)
    end
end


Initial_Condition = read_input(file_in,"Initial Condition ",String,2)
Kz = read_input(file_in,"Kz",Float64,0)
Ky = read_input(file_in,"Ky",Float64,0)
hx = read_input(file_in,"hx",Float64,0)
hy = read_input(file_in,"hy",Float64,0)
hz = read_input(file_in,"hz",Float64,0)

############ Reading Inputs (DMRG) ############
work_flow = read_input(file_in,"work",String,2)
nsweeps = read_input(file_in,"nsweeps",Int,0)
cutoff = read_input(file_in,"cutoff",Float64,0)
converg = read_input(file_in,"converg",Float64,0)
noise = read_input(file_in,"noise",Float64,1)
max_sweep = read_input(file_in,"max_sweep",Float64,0)
if work_flow == "Continue"
 band_min = read_input(file_in,"band_min",Int,0)
else
end
band_max = read_input(file_in,"band_max",Int,0)
maxdim = read_input(file_in,"maxdim",Int,1)
mindim = read_input(file_in,"mindim",Int,1)
psi_lambda = read_input(file_in,"psi_lambda",Float64,0)
close(file_in)

########## End of Reading Inputs ############

########## Writing Oututs ############
if work_flow == "Start"
    file_out = open("Output", "w")
else
    file_out = open("Output", "a")
end

write(file_out, "###############################################")
write(file_out, "\rOutput file of Quantum Spin Chain simulation")
write(file_out, "\r")

write(file_out, "\r#### Parameters #####")
write(file_out, "\r")

write(file_out, "\rLattice Model: $Lattice_type ")
write(file_out, "\rMobile_DW: $Mobile_DW ")

if Lattice_type == 1 || Lattice_type ==5
    write(file_out, "\rSite Number 'N': $N")
    if Mobile_DW == 1
        write(file_out, "\rSite of BC 'NBC': $NBC")
        #write(file_out, "\rSite of BC2 'NBC2': $NBC2")
    end

 elseif Lattice_type == 2 || Lattice_type ==3
    write(file_out, "\rSite Number 'Nx': $Nx")
    write(file_out, "\rSite Number 'Ny': $Ny")
    write(file_out, "\rBondary Condition (2D) : $BC_2D")
    write(file_out, "\rInter-Exchange Constant 'J_inter': $J_inter")
    global N=Nx*Ny
 elseif Lattice_type == 4 
    write(file_out, "\rSite Number 'Nx': $Nx")
    write(file_out, "\rSite Number 'Ny': $Ny")
    write(file_out, "\rSite of coupling 'Nc': $Nc")
    write(file_out, "\rBondary Condition (2D) : $BC_2D")
    write(file_out, "\rInter-Exchange Constant 'J_inter': $J_inter")
    global N=Nx*Ny
    if Mobile_DW == 1
        write(file_out, "\rSite of BC_Chian1 'NBC_1': $NBC_1")
        write(file_out, "\rSite of BC_Chian2 'NBC_2': $NBC_2")
    end
end



write(file_out, "\rInitial Condition: $Initial_Condition")


write(file_out, "\rExchange Constant 'J': $J")
write(file_out, "\rUniaxial Anisotropy 'Kz': $Kz")
write(file_out, "\rIn-plane Anisotropy 'Ky': $Ky")
write(file_out, "\rExternal Filed 'hx': $hx")
write(file_out, "\rExternal Filed 'hy': $hy")
write(file_out, "\rExternal Filed 'hz': $hz")

write(file_out, "\r")

write(file_out, "\r#### DMRG Parameters ####")
write(file_out, "\r")

write(file_out, "\rWork Flow: $work_flow")
write(file_out, "\rnsweeps: $nsweeps")
write(file_out, "\rcutoff: $cutoff")
write(file_out, "\rconverg: $converg")
write(file_out, "\rnoise: $noise")
write(file_out, "\rmax_sweep: $max_sweep")
if work_flow == "Continue"
    write(file_out, "\rband_min: $band_min")
   else
end
write(file_out, "\rband_max: $band_max")
write(file_out, "\rmaxdim: $maxdim")
write(file_out, "\rmindim: $mindim")
write(file_out, "\rpsi_lambda: $psi_lambda")

write(file_out, "\r")
#flush(file_out)

############ memory conuter #####

############
####### timer  ##################

write(file_out, "\rReading Inputs finished.")
DATE=[]
push!(DATE,Dates.DateTime(Dates.now()))
rightnow=DATE[end]
write(file_out, "\rDate: $rightnow")

write(file_out, "\r")
flush(file_out)

write(file_out, "\r#### Work flow Information ####")
write(file_out, "\r")

redirect_stdout(file_out)
if work_flow == "Start"
    
    ########## Constructing Spin Chain Hamiltonian ############
    if Lattice_type ==1
        if Mobile_DW == 1
           global H,sites=Heisenberg_Ham_mobile(N,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
        else
           global H,sites=Heisenberg_Ham(N,J,Kz,Ky,hx,hy,hz)
        end
    elseif Lattice_type ==2
     global H,sites=Heisenberg_Ham2D(Nx,Ny,J_inter,J,Kz,Ky,hx,hy,hz,BC_2D)
    elseif Lattice_type ==3
     global H,sites=Heisenberg_Ham_ladder(Nx,Ny,J_inter,J,Kz,Ky,hx,hy,hz,BC_2D)
    elseif Lattice_type ==4
        if Mobile_DW == 1
            global H,sites=Heisenberg_Ham_ladder_single_mobile(Nx,Ny,Nc,NBC_1[1],NBC_1[2],NBC_2[1],NBC_2[2],J_inter,J,Kz,Ky,hx,hy,hz,BC_2D)
        else
            global H,sites=Heisenberg_Ham_ladder_single(Nx,Ny,Nc,J_inter,J,Kz,Ky,hx,hy,hz,BC_2D)
        end
    elseif Lattice_type ==5
        sites = siteinds("S=1/2",N)
        os_Ham = OpSum()
        os_Ham += hz, "Sz", 1
        os_Ham += 2*hz, "Sz", 2  
        os_Ham += J, "Sz", 1, "Sz", 1
        os_Ham += J, "Sz", 2, "Sz", 2
        global H=MPO(os_Ham,sites)
    end

    ####### timer  ##################
    write(file_out, "\rConstructing Hamiltonian finished.")
    push!(DATE,Dates.DateTime(Dates.now()))
    rightnow=DATE[end]
    write(file_out, "\rDate: $rightnow")
    write(file_out, "\r")
elseif work_flow == "Continue"
    pattern = r"psi_\d+\.h5" 
    files_folder = readdir()
    psi_files = filter(file -> occursin(pattern, file), files_folder)
    #band_min=length(psi_files)+1
    psi=MPS[]
    for i=1:(band_min-1)
        local f = h5open(string("psi_",i,".h5"),"r")
        local psi_temp=read(f,"psi",MPS)
        close(f)
        push!(psi,psi_temp)
    end

    local f = h5open("Ham.h5","r")
    H=read(f,"Hamiltonian",MPO)
    close(f)

    energy=load("DMRG_data.jld2","energy")
    sweep_num=load("DMRG_data.jld2","sweep_num")
    E=load("DMRG_data.jld2","E")
    #D_E=load("data.jld2","D_E")
    sites = siteinds(psi[1])
    write(file_out, "\rContinue calculation from State $band_min.")
    push!(DATE,Dates.DateTime(Dates.now()))
    rightnow=DATE[end]
    write(file_out, "\rDate: $rightnow")
    write(file_out, "\r")
else
    write(file_out, "\rWork flow is not assigned.")
    push!(DATE,Dates.DateTime(Dates.now()))
    rightnow=DATE[end]
    write(file_out, "\rDate: $rightnow")
    write(file_out, "\r")
    exit()
end

###### Saving Data #########
file_ham = h5open("Ham.h5","w")
write(file_ham,"Hamiltonian",H)
close(file_ham)


##### Spin Operators ######
global SztotMPO,SytotMPO,SxtotMPO = Tot_spin_Op(sites)

##### Initialization for DMRG ######

if Initial_Condition == "DW" || Initial_Condition == "Variation"
    if Lattice_type ==1
    state=[n< (N+1)/2 ? "Up" : "Dn" for n=1:N]
    elseif Lattice_type == 2 || Lattice_type ==3 || Lattice_type ==4
        if BC_2D ==10
            state_1=[n< (N+1)/2 ? "Up" : "Dn" for n=1:2:N]
            state_2=[ "Up" for n=2:2:N]
            state=interleave_arrays(state_1,state_2)
        elseif BC_2D ==11
            state=[n< (N+1)/2 ? "Up" : "Dn" for n=1:N]
        elseif BC_2D ==12
            state_1=[n< (N+1)/2 ? "Up" : "Dn" for n=1:2:N]
            state_2=[n< (N+1)/2 ? "Dn" : "Up" for n=2:2:N]
            state=interleave_arrays(state_1,state_2)

        elseif BC_2D ==13
            state_1=[n< (N+1)/2 ? "Up" : "Dn" for n=1:2:N]
            state_2=[ "Up" for n=2:2:N]
            state=interleave_arrays(state_1,state_2)
        end
    end
    psi_init=MPS(sites,state);
elseif Initial_Condition == "Random"
    psi_init=randomMPS(sites,3)
else
    write(file_out, "\rInitial ondition is not assigned.")
    exit()
end
#psi_init=randomMPS(sites)


if work_flow == "Start"
    #### ground state ######
    write(file_out, "\r !!! Calculation starts from State 1.")
    println("\r##### Start from Sweep: 0  #######")
    
    ##### Initialization ########
    E=[] # energy loop data
    sweep_num=[] # sweep number
    psi=MPS[] # eigenvector
    energy=[] # eigenvalue  

        ##### pre-calculation ########
    E_loop=[]
    sweep_num_loop=[]
    global energy_temp,psi_temp,sweep_num_loop,E_loop = DMRG_loop(H,[],psi_init,nsweeps,maxdim,mindim,cutoff,noise,sweep_num_loop,E_loop,0,converg)
         
         ##### dmrg loop ########
         
         for j=1:10
            nsweeps_dum=nsweeps*2^j
            global energy_temp,psi_temp,sweep_num_loop,E_loop =DMRG_loop(H,[],psi_temp,nsweeps_dum,maxdim,mindim,cutoff,0,sweep_num_loop,E_loop,0,converg)
            if energy_diff < converg || nsweeps_dum >= max_sweep
                break
            end
            GC.gc()
         end

    ##### save for loops ########
    push!(E,E_loop)
    push!(sweep_num,sweep_num_loop)
    push!(psi,psi_temp)
    push!(energy,energy_temp)
        
    ##### save for psi ###########
    local f = open("Ene_psi_1.txt", "w") 
    writedlm(f, [sweep_num[1] E[1]], ' ')
    close(f)

    i = 1
    file_psi = h5open(string("psi_",i,".h5"),"w")
    write(file_psi,"psi",psi[i])
    close(file_psi)
    
    jldsave("DMRG_data.jld2"; energy,E,sweep_num)

    write(file_out, "\rState 1 finished.")
    local rightnow=DATE[end]
    write(file_out, "\rDate: $rightnow")
    write(file_out, "\r")
    if sweep_num_loop[end] >= max_sweep
        write(file_out, "\rDmrg loops of state 1 exceeds maxmium = $max_sweep")
        push!(DATE,Dates.DateTime(Dates.now()))
        close(file_out)
        exit()
    else
    end
    GC.gc()


    #### excited state (i) ######
    band_min=2
   
    for i=band_min:band_max
        write(file_out, "\r !!! Calculation starts from State $i.")
        println("\r##### Start from Sweep: 0  #######")
         E_loop=[]
         sweep_num_loop=[]
         if  Initial_Condition ==  "Variation"
            global psi_init = normalize!(psi[i-1]+psi_lambda*randomMPS(sites,mindim[1]))
         end
        ##### pre-calculation ########
         global energy_temp,psi_temp,sweep_num_loop,E_loop = DMRG_loop(H,psi[1:i-1],psi_init,nsweeps,maxdim,mindim,cutoff,noise,sweep_num_loop,E_loop,1,converg)
    
        ##### dmrg loop ########
        for j=1:10
            nsweeps_dum=nsweeps*2^j
            global energy_temp,psi_temp,sweep_num_loop,E_loop =DMRG_loop(H,psi[1:i-1],psi_temp,nsweeps_dum,maxdim,mindim,cutoff,0,sweep_num_loop,E_loop,1,converg)
            if energy_diff < converg || sweep_num_loop[end] >= max_sweep
                break
            end
            GC.gc()
         end
        ##### save for loops ########
        if length(psi)< i
            push!(E,E_loop)
            push!(sweep_num,sweep_num_loop)
            push!(psi,psi_temp)
            push!(energy,energy_temp)
            write(file_out, "\rState $i finished.")
            
        else
            E[i]=E_loop
            sweep_num[i]=sweep_num_loop
            psi[i]=psi_temp
            energy[i]=energy_temp
            write(file_out, "\rState $i has been rewritten.")
            
        end

        ###### save for psi ##################
        local f = open("Ene_psi_"*string(i)*".txt", "w") 
        writedlm(f, [sweep_num[i] E[i]], ' ')
        close(f)

        local file_psi = h5open(string("psi_",i,".h5"),"w")
        write(file_psi,"psi",psi[i])
        close(file_psi)

        jldsave("DMRG_data.jld2"; energy,E,sweep_num)

        ####### timer  ##################
        push!(DATE,Dates.DateTime(Dates.now()))
        local rightnow=DATE[end]
        write(file_out, "\rDate: $rightnow")
        write(file_out, "\r")

        if sweep_num_loop[end] >= max_sweep
            write(file_out, "\rDmrg loops of state $i exceeds maxmium = $max_sweep")
            break
        else
        end

        GC.gc()
    end

elseif work_flow == "Continue"

    for i=band_min:band_max
        write(file_out, "\rCalculation starts from State $i.")
        write(file_out, "\r")
         E_loop=[]
         sweep_num_loop=[]
         if  Initial_Condition ==  "Variation"
            global psi_init = normalize!(psi[i-1]+psi_lambda*randomMPS(sites,mindim[1]))
         end
        ##### pre-calculation ########
        #global energy_temp,psi_temp =@time dmrg(H,psi_init; nsweeps, maxdim,cutoff,noise,mindim, observer=obs,outputlevel=1,eigsolve_krylovdim=5);
        #push!(sweep_num_loop,nsweeps)
        #push!(E_loop,energy_temp)
         global energy_temp,psi_temp,sweep_num_loop,E_loop = DMRG_loop(H,psi[1:i-1],psi_init,nsweeps,maxdim,mindim,cutoff,noise,sweep_num_loop,E_loop,1,0)
   
       ##### dmrg loop ########
       for j=1:10
         nsweeps_dum=nsweeps*2^j
         global energy_temp,psi_temp,sweep_num_loop,E_loop =DMRG_loop(H,psi[1:i-1],psi_temp,nsweeps_dum,maxdim,mindim,cutoff,0,sweep_num_loop,E_loop,1,converg)
         if energy_diff < converg || sweep_num_loop[end] >= max_sweep
            break
         end
        end
        push!(psi,psi_temp)
        if length(E)< i
            push!(E,E_loop)
            push!(sweep_num,sweep_num_loop)
            push!(energy,energy_temp)
            write(file_out, "\rState $i finished.")
            
        else
            E[i]=E_loop
            sweep_num[i]=sweep_num_loop
            energy[i]=energy_temp
            write(file_out, "\rState $i has been rewritten.")
            
        end

        ###### save for psi ##################
        local f = open("Ene_psi_"*string(i)*".txt", "w") 
        writedlm(f, [sweep_num[i] E[i]], ' ')
        close(f)
        local file_psi = h5open(string("psi_",i,".h5"),"w")
        write(file_psi,"psi",psi[i])
        close(file_psi)

        jldsave("DMRG_data.jld2"; energy,E,sweep_num)

        ####### timer  ##################
        push!(DATE,Dates.DateTime(Dates.now()))
        local rightnow=DATE[end]
        write(file_out, "\rDate: $rightnow")
        write(file_out, "\r")

        if sweep_num_loop[end] >= max_sweep
            write(file_out, "\rDmrg loops of state $i exceeds maxmium = $max_sweep")
            write(file_out, "\r")
            break
        else
        end

        GC.gc()
    end

end





######## Spin profile Plot #############
header = ["Sz" "Sx" "Sy" ]
#global plot_Si_array_1=[]
#global plot_Si_array_2=[]
if     Lattice_type == 1
    for i in eachindex(psi)
        
        Szlocal=[]
        Sylocal=[]
        Sxlocal=[]

        Szlocal=expect(psi[i],"Sz");
        Sylocal=expect(psi[i],"Sy");
        Sxlocal=expect(psi[i],"Sx");

        site_num=Int(length(Szlocal));

        open("Spindata_hy"*string(hy)*"_hx"*string(hx)*"state"*string(i)*".txt"; write=true) do f
            writedlm(f, header)
            writedlm(f, [real(Szlocal[1:site_num]) real(Sxlocal[1:site_num]) real(Sylocal[1:site_num]) ])
        end
    end

elseif Lattice_type == 2
    for i in eachindex(psi)

        Szlocal=[]
        Sylocal=[]
        Sxlocal=[]

        Szlocal=expect(psi[i],"Sz");
        Sylocal=expect(psi[i],"Sy");
        Sxlocal=expect(psi[i],"Sx");

        site_num=Int(length(Szlocal));

        open("Spindata_hy"*string(hy)*"_hx"*string(hx)*"state"*string(i)*"_chain1.txt"; write=true) do f
            writedlm(f, header)
            writedlm(f, [real(Szlocal[1:2:end]) real(Sxlocal[1:2:end]) real(Sylocal[1:2:end]) ])
        end

        open("Spindata_hy"*string(hy)*"_hx"*string(hx)*"state"*string(i)*"_chain2.txt"; write=true) do f
            writedlm(f, header)
            writedlm(f, [real(Szlocal[2:2:end]) real(Sxlocal[2:2:end]) real(Sylocal[2:2:end]) ])
        end
    end


elseif  Lattice_type == 3 || Lattice_type ==4
    for i in eachindex(psi)

        Szlocal=[]
        Sylocal=[]
        Sxlocal=[]

        Szlocal=expect(psi[i],"Sz");
        Sylocal=expect(psi[i],"Sy");
        Sxlocal=expect(psi[i],"Sx");
        site_num=Int(length(Szlocal)/2)

        open("Spindata_hy"*string(hy)*"_hx"*string(hx)*"state"*string(i)*"_chain1.txt"; write=true) do f
            writedlm(f, header)
            writedlm(f, [real(Szlocal[1:site_num]) real(Sxlocal[1:site_num]) real(Sylocal[1:site_num]) ])
        end

        open("Spindata_hy"*string(hy)*"_hx"*string(hx)*"state"*string(i)*"_chain2.txt"; write=true) do f
            writedlm(f, header)
            writedlm(f, [real(Szlocal[site_num+1:site_num*2]) real(Sxlocal[site_num+1:site_num*2]) real(Sylocal[site_num+1:site_num*2]) ])
        end
    end

end




################# Domain Wall Chirality #############
if Lattice_type ==1
    global DWxMPO,DWyMPO,DWzMPO=DWC_operator_1D(N::Int,sites)


    global  DWxData=[]
    global  DWyData=[]
    global  DWzData=[]
    for i in eachindex(psi)

        DWxtemp=inner(psi[i]',DWxMPO,psi[i])
        DWytemp=inner(psi[i]',DWyMPO,psi[i])
        DWztemp=inner(psi[i]',DWzMPO,psi[i])

        push!(DWxData,DWxtemp)
        push!(DWyData,DWytemp)
        push!(DWzData,DWztemp)

    end

    header = ["Cx_Re" "Cx_Im" "Cy_Re" "Cy_Im" "Cz_Re" "Cz_Im"]
    open("DWdata_hy"*string(hy)*"_hx"*string(hx)*".txt"; write=true) do f
        writedlm(f, header)
        writedlm(f, [real(DWxData) imag(DWxData) real(DWyData) imag(DWyData) real(DWzData) imag(DWzData) ])
    end
elseif Lattice_type ==3 || Lattice_type ==4
    global DWxMPO1,DWyMPO1,DWzMPO1,DWxMPO2,DWyMPO2,DWzMPO2=DWC_operator_2D(Nx::Int,Ny::Int,sites)


    global DWxData=[]
    global DWyData=[]
    global DWzData=[]
    for i in eachindex(psi)

        DWxtemp=inner(psi[i]',DWxMPO1,psi[i])
        DWytemp=inner(psi[i]',DWyMPO1,psi[i])
        DWztemp=inner(psi[i]',DWzMPO1,psi[i])

        push!(DWxData,DWxtemp)
        push!(DWyData,DWytemp)
        push!(DWzData,DWztemp)

    end

    global DWxData2=[]
    global DWyData2=[]
    global DWzData2=[]
    for i in eachindex(psi)

        DWxtemp=inner(psi[i]',DWxMPO2,psi[i])
        DWytemp=inner(psi[i]',DWyMPO2,psi[i])
        DWztemp=inner(psi[i]',DWzMPO2,psi[i])

        push!(DWxData2,DWxtemp)
        push!(DWyData2,DWytemp)
        push!(DWzData2,DWztemp)

    end


    using DelimitedFiles
    header = ["Cx1_Re" "Cx1_Im" "Cy1_Re" "Cy1_Im" "Cz1_Re" "Cz1_Im" "Cx2_Re" "Cx2_Im" "Cy2_Re" "Cy2_Im" "Cz2_Re" "Cz2_Im"]
    open("DWdata_hy"*string(hy)*"_hx"*string(hx)*".txt"; write=true) do f
        writedlm(f, header)
        writedlm(f, [real(DWxData) imag(DWxData) real(DWyData) imag(DWyData) real(DWzData) imag(DWzData) real(DWxData2) imag(DWxData2) real(DWyData2) imag(DWyData2) real(DWzData2) imag(DWzData2) ])
    end

end

###### Saving Data #########
   
jldsave("DMRG_data.jld2"; energy,E,sweep_num)

file_ene = open("Ene_data_hy"*string(hy)*"_hx"*string(hx)*".txt", "w") 
writedlm(file_ene, energy)
close(file_ene)

####### timer  ##################
write(file_out, "\rSimulation Finished.")
push!(DATE,Dates.Time(Dates.now()))
rightnow=DATE[end]
runtime=round(Dates.value(DATE[end]-DATE[1])/1E3/60;digits = 2)
write(file_out, "\rDate: $rightnow")
write(file_out, "\rtotal runtime: $runtime mins")
write(file_out, "\r###############################################")
close(file_out)
#GC.gc()

########## End of Writing Outputs ###########