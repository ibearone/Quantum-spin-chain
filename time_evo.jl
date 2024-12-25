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
using ITensorTDVP
using Observers
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
    time_evo_method = read_input(file_in,"time_evo_method",String,2)
    write_psi_evo = read_input(file_in,"write_psi_evo",Int,0)
    tau = read_input(file_in,"tau_dis",Float64,0)
    t_total = read_input(file_in,"t_total",Float64,0)
    band_evo = read_input(file_in,"band_evo",Int,0)
   else
end
close(file_in)

########## End of Reading Inputs ############

########## Writing Oututs ############
if work_flow == "time_evo"
 file_out = open("Output_evo", "w")
 else
  write(file_out, "\rWork flow is not assigned.")
  push!(DATE,Dates.Time(Dates.now()))
  rightnow=DATE[end]
  write(file_out, "\rDate: $rightnow")
  write(file_out, "\r")
  exit()
end

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

write(file_out, "\rband_max: $band_max")
write(file_out, "\rmaxdim: $maxdim")
write(file_out, "\rmindim: $mindim")
write(file_out, "\rpsi_lambda: $psi_lambda")

write(file_out, "\r")
write(file_out, "\r#### Time evolution Parameters ####")
write(file_out, "\r")

if work_flow == "time_evo"
    write(file_out, "\rtime_evo_method: $time_evo_method")
    write(file_out, "\rwrite psi_evo: $write_psi_evo")
    write(file_out, "\rtau_dis: $tau")
    write(file_out, "\rt_total: $t_total")
    write(file_out, "\rband_evo: $band_evo")
   else
end
write(file_out, "\r")


write(file_out, "\r#### Reading DMRG data ####")
write(file_out, "\r")

###### Reading Data #########
psi=[]
for i=1:band_max
    f = h5open(string("psi_",i,".h5"),"r")
    psi_temp=read(f,"psi",MPS)
    close(f)
    push!(psi,psi_temp)
end

sites = siteinds(psi[1])

write(file_out, "\rRead psi_$band_max finished.")
write(file_out, "\r")


f = h5open("Ham.h5","r")
H=read(f,"Hamiltonian",MPO)
close(f)

write(file_out, "\rRead Hamiltonian finished.")
write(file_out, "\r")

energy=load("DMRG_data.jld2","energy")
E=load("DMRG_data.jld2","E")
sweep_num=load("DMRG_data.jld2","sweep_num")

write(file_out, "\rRead jld2 finished.")
write(file_out, "\r")


write(file_out, "\r#### Define gates for time evolution ####")
write(file_out, "\r")
flush(file_out)

####### define time evolution operator #######
if Mobile_DW == 1
  global H_time,gates=Ham_mobile_gates(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)

else
  NBC=[1,N]
  global H_time,gates=Ham_mobile_gates(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
end
global DWxMPO,DWyMPO,DWzMPO=DWC_operator_1D(N::Int,sites)


####### time evolution loops #######

write(file_out, "\r#### Start calculation time evolution of state $band_evo ####")
write(file_out, "\r")
flush(file_out)

####### Initialization ########
if write_psi_evo == 1
  psi_evo=MPS[]
end


Ene_H0=[]
Ene_H_time=[]
S_site=[]
DW_C=[]
DATE =[]

if time_evo_method == "TEBD"
    global psi_temp = apply(gates, psi[band_evo]; cutoff)
    if write_psi_evo == 1
      push!(psi_evo,psi_temp)
    end
      
    for (i, t) in enumerate(0.0:tau:t_total)
      #psi_t = psi_evo[i]

      ene_temp=real(inner(psi_temp', H,psi_temp))
      ene_temp2=real(inner(psi_temp', H_time,psi_temp))

      Sz = expect(psi_temp, "Sz")
      Sy = expect(psi_temp, "Sy")
      Sx = expect(psi_temp, "Sx")

      DWztemp=inner(psi_temp',DWzMPO,psi_temp)
      DWytemp=inner(psi_temp',DWyMPO,psi_temp)
      DWxtemp=inner(psi_temp',DWxMPO,psi_temp)

      push!(Ene_H0,ene_temp)
      push!(Ene_H_time,ene_temp2)

      #push!(S_site,(Sz,Sy,Sx))
      push!(S_temp_z,Sz)
      push!(S_temp_y,Sy)
      push!(S_temp_x,Sx)
      #push!(DW_C,(DWztemp,DWytemp,DWxtemp))
      push!(DW_temp_z,DWztemp)
      push!(DW_temp_y,DWytemp)
      push!(DW_temp_x,DWxtemp)

      t≈t_total && break
        
      global psi_temp = apply(gates, psi_temp; cutoff)
      normalize!(psi_temp)
      if write_psi_evo == 1
        push!(psi_evo,psi_temp)
      end

        if (i-1) % Int(t_total/tau/100) == 0 && i != 0
          write(file_out, "\rTime step: $t/$t_total")
          push!(DATE,Dates.Time(Dates.now()))
          local rightnow=DATE[end]
          write(file_out, "\rDate: $rightnow")
          write(file_out, "\r")
          flush(file_out)
          GC.gc()
      end
    end
    S_site=(S_temp_z,S_temp_y,S_temp_x)
    DW_C=(DW_temp_z,DW_temp_y,DW_temp_x)

 elseif time_evo_method == "TDVP"
  write(file_out, "\r!!! Start TDVP Calculation !!!")
  write(file_out, "\r")
    step(; sweep) = sweep
    #begin
      #if sweep % Int(t_total/tau/100) == 0 && sweep != 0
      #  println(file_out,"Current step: ", sweep)
      #end
    #  return 
    #end
    current_time(; sweep,current_time) = current_time
    #begin
      # if sweep % Int(t_total/tau/100) == 0 && sweep != 0
      #   println(file_out,"Current time: ", round(real(current_time*im),digits=2),"/",t_total, "(real time)")
      # end
    #  return current_time
    #end
    return_state(; state) = state
    measure_Ene0(; sweep,state)= real(inner(state', H, state))
    #begin
    #   local val = 
      # if sweep % Int(t_total/tau/100) == 0 && sweep != 0
      #   println(file_out,"Ene0 = ", round(val,digits=8))
      # end
    #  return val
    #end
    measure_Ene_time(; sweep,state)= real(inner(state', H_time, state))
    #begin
      local val = 
      # if sweep % Int(t_total/tau/100) == 0 && sweep != 0
      #   println(file_out,"Ene_time = ", round(val,digits=8))
      # end
    #  return val
    #end

    timer(; sweep, current_time, state) = begin
      if sweep % Int(t_total/tau/100) == 0 && sweep != 0
        push!(DATE,Dates.Time(Dates.now()))
        local rightnow=DATE[end]
        println(file_out,"Time step: ", round(real(current_time*im),digits=2),"/",t_total,"   Ene0 = ", round(real(inner(state', H, state)),digits=8),"   Ene_time = ", round(real(inner(state', H_time, state)),digits=8),"   Date: ",rightnow)
        write(file_out, "\r")
        flush(file_out)
      end
      return nothing
    end

    measure_sz(; state) = expect(state, "Sz")
    measure_sy(; state) = expect(state, "Sy")
    measure_sx(; state) = expect(state, "Sx")
    measure_Cz(; state) = inner(state',DWzMPO,state)
    measure_Cy(; state) = inner(state',DWyMPO,state)
    measure_Cx(; state) = inner(state',DWxMPO,state)

    obs = observer(
      "steps" => step, "times" => current_time, "states" => return_state, "Ene0" => measure_Ene0, "Ene_time" => measure_Ene_time,
       "sz" => measure_sz, "sy" => measure_sy, "sx" => measure_sx,
       "Cz" => measure_Cz, "Cy" => measure_Cy, "Cx" => measure_Cx,"timer" => timer
    )

    state = tdvp(H_time, -im*t_total, psi[band_evo]; time_step=-im*tau, cutoff, (step_observer!)=obs, outputlevel=0)
    Ene_H0 = obs.Ene0
    Ene_H_time = obs.Ene_time
    S_site=(obs.sz,obs.sy,obs.sx)
    DW_C=(obs.Cz,obs.Cy,obs.Cx)
    if write_psi_evo == 1
     psi_evo=obs.states
    end
    # if obs.steps % Int(t_total/tau/100) == 0 && obs.steps != 0
    #   local t=obs.times
    #   write(file_out, "\rTime step: $t/$t_total")
    #   push!(DATE,Dates.Time(Dates.now()))
    #   local rightnow=DATE[end]
    #   write(file_out, "\rDate: $rightnow")
    #   write(file_out, "\r")
    #   flush(file_out)
    #   GC.gc()
    # end
else

end

write(file_out, "\rTime step: $t_total/$t_total")
push!(DATE,Dates.Time(Dates.now()))
rightnow=DATE[end]
write(file_out, "\rDate: $rightnow")
write(file_out, "\r")
write(file_out, "\rFinished time evolution of state $band_evo")
flush(file_out)
GC.gc()

###### Saving psi_evo Data #########

if write_psi_evo == 1
 file_psi = h5open(string("psi_evo_", band_evo, ".h5"), "w")

 for i in 1:length(psi_evo)
   evo_name = "psi_evo_$(i)"
  
   write(file_psi, evo_name, psi_evo[i])
 end
   psi_evo_length=length(psi_evo)
   close(file_psi)


   write(file_out, "\r#### Psi_evo_$band_evo Saved ####")
   write(file_out, "\r")
   write(file_out, "\rPsi_evo_$band_evo length : $psi_evo_length")
   write(file_out, "\r")
   flush(file_out)
 else
  
end

###### Saving Data #########
   
jldsave("time_evo_data.jld2"; Ene_H0,Ene_H_time,S_site,DW_C)

####### timer  ##################
write(file_out, "\rSimulation Finished.")
push!(DATE,Dates.Time(Dates.now()))
rightnow=DATE[end]
runtime=round(Dates.value(DATE[end]-DATE[1])/1E9/60;digits = 2)
write(file_out, "\rDate: $rightnow")
write(file_out, "\rtotal runtime: $runtime mins")
write(file_out, "\r###############################################")
close(file_out)
