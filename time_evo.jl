#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the  code for time evolution on   #
# spin chain simulation using ITensor       #
# library.                                  #
#                                           #
#                                           #
############################################# 

########## Initialization ############
using ITensors
using ITensorMPS
using HDF5
#using ITensors.HDF5
#using LaTeXStrings
#using Plots
using JLD2
#using DelimitedFiles
using Dates

########### Modules ##############
include("read_input.jl")

include("Hamiltonian.jl")
include("time_evo_gate.jl")
include("ladder_lattice.jl")
include("Spin_Operator.jl")
include("DW_chirality_operators.jl")

########## Reading Inputs (Hamiltonian) ############

file_in = open("Input", "r")

Lattice_type = read_input(file_in,"Lattice",Int,0)
Mobile_DW = read_input(file_in,"Mobile_DW",Int,0)
if Lattice_type ==1
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

band_max = read_input(file_in,"band_max",Int,0)
maxdim = read_input(file_in,"maxdim",Int,1)
mindim = read_input(file_in,"mindim",Int,1)
psi_lambda = read_input(file_in,"psi_lambda",Float64,0)
if work_flow == "time_evo"
    write_psi_evo = read_input(file_in,"write_psi_evo",Int,0)
    tau = read_input(file_in,"tau_dis",Float64,0)
    t_total = read_input(file_in,"t_total",Float64,0)
    band_evo = read_input(file_in,"band_evo",Int,0)
   else
end
close(file_in)

########## End of Reading Inputs ############

########## Writing Oututs ############
#if work_flow == "Start"
file_out = open("Output_evo", "w")
#else
    #file_out = open("Output", "a")
#end

write(file_out, "###############################################")
write(file_out, "\rOutput file of Quantum Spin Chain simulation (time evolution)")
write(file_out, "\r")

write(file_out, "\r#### Parameters #####")
write(file_out, "\r")

write(file_out, "\rLattice Model: $Lattice_type ")
write(file_out, "\rMobile_DW: $Mobile_DW ")

if Lattice_type == 1
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
# if work_flow == "Continue"
#     write(file_out, "\rband_min: $band_min")
#    else
# end
write(file_out, "\rband_max: $band_max")
write(file_out, "\rmaxdim: $maxdim")
write(file_out, "\rmindim: $mindim")
write(file_out, "\rpsi_lambda: $psi_lambda")

write(file_out, "\r#### Time evolution Parameters ####")
write(file_out, "\r")

if work_flow == "time_evo"
    write(file_out, "\rwrite psi_evo: $write_psi_evo")
    write(file_out, "\rtau_dis: $tau")
    write(file_out, "\rt_total: $t_total")
    write(file_out, "\rband_evo: $band_evo")
   else
end
write(file_out, "\r")
#flush(file_out)

############ memory conuter #####

############
####### timer  ##################

write(file_out, "\rReading Inputs finished.")

write(file_out, "\r")
flush(file_out)


###### Reading Data #########
psi=[]
for i=1:band_max
    f = h5open(string("psi_",i,".h5"),"r")
    psi_temp=read(f,"psi",MPS)
    close(f)
    push!(psi,psi_temp)
end

sites = siteinds(psi[1])

write(file_out, "\rReading Psi_* finished.")

write(file_out, "\r")
flush(file_out)

f = h5open("Ham.h5","r")
H=read(f,"Hamiltonian",MPO)
close(f)

write(file_out, "\rReading Hamiltonian finished.")

write(file_out, "\r")
flush(file_out)

energy=load("DMRG_data.jld2","energy")
E=load("DMRG_data.jld2","E")
sweep_num=load("DMRG_data.jld2","sweep_num")

write(file_out, "\rReading jld2 finished.")

write(file_out, "\r")
flush(file_out)

###### Reading Data #########


DATE=[]
push!(DATE,Dates.Time(Dates.now()))
rightnow=DATE[end]
write(file_out, "\rDate: $rightnow")

write(file_out, "\r#### Define gates for time evolution ####")
write(file_out, "\r")


####### define time evolution operator #######
if Mobile_DW == 1
  global H_time,gates=Ham_mobile_gates(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)

else
  NBC=[1,N]
  global H_time,gates=Ham_mobile_gates(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
end
global DWxMPO,DWyMPO,DWzMPO=DWC_operator_1D(N::Int,sites)
write(file_out, "\r#### Work flow Information ####")
write(file_out, "\r")




####### define time evolution operator #######

write(file_out, "\r#### Start calculation time evolution of state $band_evo ####")
write(file_out, "\r")

psi_evo=MPS[]

Ene_H0=[]
Ene_H_time=[]

Sz_i=[]
Sy_i=[]
Sx_i=[]

Cx=[]
Cy=[]
Cz=[]

global psi_temp = apply(gates, psi[band_evo]; cutoff)
push!(psi_evo,psi_temp)
for (i, t) in enumerate(0.0:tau:t_total)
  psi_t = psi_evo[i]

  ene_temp=real(inner(psi_t', H,psi_t))
  ene_temp2=real(inner(psi_t', H_time,psi_t))

  Sz = expect(psi_t, "Sz")
  Sy = expect(psi_t, "Sy")
  Sx = expect(psi_t, "Sx")
  DWxtemp=inner(psi_t',DWxMPO,psi_t)
  DWytemp=inner(psi_t',DWyMPO,psi_t)
  DWztemp=inner(psi_t',DWzMPO,psi_t)

  
  push!(Ene_H0,ene_temp)
  push!(Ene_H_time,ene_temp2)

  push!(Sz_i,Sz)
  push!(Sy_i,Sy)
  push!(Sx_i,Sx)
  
  push!(Cx,DWxtemp)
  push!(Cy,DWytemp)
  push!(Cz,DWztemp)

  t≈t_total && break
    
  global psi_temp = apply(gates, psi_temp; cutoff)
  normalize!(psi_temp)
  push!(psi_evo,psi_temp)

    if i+1 % 10 == 0 && i != 0
      write(file_out, "\rtime step: $t")
      push!(DATE,Dates.Time(Dates.now()))
      local rightnow=DATE[end]
      write(file_out, "\rDate: $rightnow")
      write(file_out, "\r")
      flush(file_out)
  end
end



if write_psi_evo == 1
 # Open an HDF5 file for writing with a name based on band_evo
 file_psi = h5open(string("psi_evo_", band_evo, ".h5"), "w")

 # Iterate over each MPS in the psi_evo array
 for i in 1:length(psi_evo)
   # Create a unique name for each MPS entry
   evo_name = "psi_evo_$(i)"
  
   # Write the MPS to the HDF5 file under the unique name
   write(file_psi, evo_name, psi_evo[i])
 end
   psi_evo_length=length(psi_evo)
   # Close the HDF5 file
   close(file_psi)
   write(file_out, "\r#### Psi_evo_$band_evo Saved ####")
   write(file_out, "\r")
   write(file_out, "\rPsi_evo_$band_evo length : $psi_evo_length")
   write(file_out, "\r")
else
  
end

###### Saving Data #########
   
jldsave("time_evo_data.jld2"; Ene_H0,Ene_H_time,Sz_i,Sy_i,Sx_i,Cz,Cy,Cx)

####### timer  ##################
write(file_out, "\rSimulation Finished.")
push!(DATE,Dates.Time(Dates.now()))
rightnow=DATE[end]
runtime=round(Dates.value(DATE[end]-DATE[1])/1E9/60;digits = 2)
write(file_out, "\rDate: $rightnow")
write(file_out, "\rtotal runtime: $runtime mins")
write(file_out, "\r###############################################")
close(file_out)
