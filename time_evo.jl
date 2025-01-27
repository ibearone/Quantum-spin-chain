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
using JLD2
using Dates
using ITensorTDVP
using Observers

#using ITensors.HDF5
#using LaTeXStrings
#using Plots
#using DelimitedFiles

########### Modules ##############
include("read_input.jl")
include("Hamiltonian.jl")
include("time_evo_gate.jl")
include("ladder_lattice.jl")
include("Spin_Operator.jl")
include("DW_chirality_operators.jl")
include("updaters.jl")
########## Reading Inputs (Hamiltonian) ############

file_in = open("Input", "r")

Lattice_type = read_input(file_in,"Lattice",Int,0)
Mobile_DW = read_input(file_in,"Mobile_DW",Int,0)
if Lattice_type ==1 || Lattice_type == 5
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
    continue_evo = read_input(file_in,"continue_evo",Int,0)
    if continue_evo == 1
      t_end = read_input(file_in,"t_end",Float64,0) 
    else
    end
    band_evo = read_input(file_in,"band_evo",Int,0)
   else
end

if time_evo_method == "TDVP_Ht" || time_evo_method == "TEBD_Ht"
  dhx = read_input(file_in,"dhx",Float64,0)
  dhy = read_input(file_in,"dhy",Float64,0)
  omega = read_input(file_in,"omega",Float64,0)
  nsite = read_input(file_in,"nsite",Int,0)
  band_tar = read_input(file_in,"band_tar",Int,0)
end
close(file_in)

########## End of Reading Inputs ############

DATE =[]
########## Writing Oututs ############
if work_flow == "time_evo"
  file_out = open("Output_evo", "w")

  else
    file_out = open("Output_evo", "w")
   write(file_out, "\rWork flow (time_evo) is not assigned.")
   push!(DATE,Dates.DateTime(Dates.now()))
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
    write(file_out, "\rcontinue_evo: $continue_evo")
    if continue_evo == 1
      write(file_out, "\rt_end: $t_end")
    else
    end
    write(file_out, "\rband_evo: $band_evo")
   else
end

if time_evo_method == "TDVP_Ht" || time_evo_method == "TEBD_Ht"
  write(file_out, "\rdhx: $dhx")
  write(file_out, "\rdhy: $dhy")
  write(file_out, "\romega: $omega")
  write(file_out, "\rnsite: $nsite")
  write(file_out, "\rband_tar: $band_tar")
end

write(file_out, "\r")
push!(DATE,Dates.DateTime(Dates.now()))
rightnow=DATE[end]
write(file_out, "\rDate: $rightnow")
write(file_out, "\r")

write(file_out, "\r#### Reading DMRG data ####")
write(file_out, "\r")

###### Reading Data #########
psi=[]
for i=1:band_max
 file_psi = h5open(string("psi_",i,".h5"),"r")
 psi_temp=read(file_psi,"psi",MPS)
  close(file_psi)
  push!(psi,psi_temp)
end
write(file_out, "\rRead psi_$band_max finished.")
write(file_out, "\r")

if continue_evo == 1
   file_psi_evo_end = h5open(string("psi_evo_end_",band_evo,".h5"),"r")
  global psi_evo_end=read(file_psi_evo_end,"psi",MPS)
  close(file_psi_evo_end)
  write(file_out, "\rRead psi_evo_end_$band_max finished.")
  write(file_out, "\r")
else
end

sites = siteinds(psi[1])

file_ham = h5open("Ham.h5","r")
H=read(file_ham,"Hamiltonian",MPO)
close(file_ham)

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
if Lattice_type == 1 
  if Mobile_DW == 1
    global H_evo = Heisenberg_Ham_mobile_2(sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
  else
    NBC=[1,N]
    global H_evo = Heisenberg_Ham_mobile_2(sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
  end
  global DWxMPO,DWyMPO,DWzMPO=DWC_operator_1D(N,sites)

elseif Lattice_type ==5

  global H_evo=H

  global DWxMPO,DWyMPO,DWzMPO=DWC_operator_1D(N,sites)
else
  write(file_out, "\rHamiltonian not assigned. (Lattice_type is limited to 1 or 5.)")
  push!(DATE,Dates.DateTime(Dates.now()))
  rightnow=DATE[end]
  write(file_out, "\rDate: $rightnow")
  write(file_out, "\r")
  exit()
end


####### time evolution loops #######

write(file_out, "\r#### Start calculation time evolution of state $band_evo ####")
write(file_out, "\r")
flush(file_out)

####### Initialization ########

if continue_evo == 0
  Ene_H0=[]
  Ene_H_evo=[]
  S_site=[]
  DW_C=[]
  p01=[]
  global psi_init=psi[band_evo]
elseif continue_evo == 1
  ################ reading previous data ##############
  if time_evo_method == "TEBD" || time_evo_method == "TDVP" || time_evo_method == "TDVP_Im_time"
    Ene_H0=load("time_evo_data.jld2","Ene_H0")
    Ene_H_evo=load("time_evo_data.jld2","Ene_H_evo")
    S_site=load("time_evo_data.jld2","S_site")
    DW_C=load("time_evo_data.jld2","DW_C")
  elseif time_evo_method == "TDVP_Ht"
    Ene_H0=load("time_evo_data.jld2","Ene_H0")
    S_site=load("time_evo_data.jld2","S_site")
    DW_C=load("time_evo_data.jld2","DW_C")
    p01=load("time_evo_data.jld2","p01")
  end

  global psi_init=psi_evo_end

  if write_psi_evo == 1
    file_psi_evo = h5open(string("psi_evo_", band_evo, ".h5"), "r")      
    psi_evo_length = detect_psi_evo_length(file_psi_evo)

    for i in 1:psi_evo_length
      evo_name = "psi_evo_$(i)"        
        local mps = read(file_psi_evo, evo_name,MPS)
  
      push!(psi_evo, mps)
    end

    close(file_psi_evo)
  end

else

end

if write_psi_evo == 1
  psi_evo=MPS[]
end

####### End of initialization ########




if time_evo_method == "TEBD" 
    obs_Ene0 = []
    obs_Ene_evo = []
    obs_sz = []
    obs_sy = []
    obs_sx = []
    obs_Cz = []
    obs_Cy = []
    obs_Cx = []
    if Lattice_type == 1 
      if Mobile_DW == 1
        global evo_gates = evo_gates_TEBD_LT_1(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
      else
        NBC=[1,N]
        global evo_gates = evo_gates_TEBD_LT_1(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz)
      end
    elseif Lattice_type == 5
      global evo_gates = evo_gates_TEBD_LT_5(N,sites,hx,hy,hz)
    else
    end
 


    global psi_temp = apply(evo_gates, psi_init; cutoff)
      
    for (i, t) in enumerate(tau:tau:t_total)


      push!( obs_Ene0,real(inner(psi_temp', H,psi_temp)))
      push!( obs_Ene_evo,real(inner(psi_temp', H_evo, psi_temp)))
      push!( obs_sz,expect(psi_temp, "Sz"))
      push!( obs_sy,expect(psi_temp, "Sy"))
      push!( obs_sx,expect(psi_temp, "Sx"))
      push!( obs_Cz,inner(psi_temp',DWzMPO,psi_temp))
      push!( obs_Cy,inner(psi_temp',DWyMPO,psi_temp))
      push!( obs_Cx,inner(psi_temp',DWxMPO,psi_temp))
      if write_psi_evo == 1
      global psi_evo=vcat(psi_evo, psi_temp)
      end
      
      t≈t_total && break
        
      global psi_temp = apply(evo_gates, psi_temp; cutoff)
      normalize!(psi_temp)
      if Int(round(t_total/tau/100)) ==0
        push!(DATE,Dates.DateTime(Dates.now()))
        local rightnow=DATE[end]
        local E0_print=round(obs_Ene0[i],digits=8)
        local Et_print=round( obs_Ene_evo[i],digits=8)
        write(file_out, "\rTime step: $t/$t_total   Ene0 = $E0_print   Ene0 = $Et_print  Date: $rightnow")
        write(file_out, "\r")
        flush(file_out)
        GC.gc()
      else
        if i % Int(t_total/tau/100) == 0 && i != 0
          push!(DATE,Dates.DateTime(Dates.now()))
          local rightnow=DATE[end]
          local E0_print=round(obs_Ene0[i],digits=8)
          local Et_print=round( obs_Ene_evo[i],digits=8)
          write(file_out, "\rTime step: $t/$t_total   Ene0 = $E0_print   Ene0 = $Et_print  Date: $rightnow")
          write(file_out, "\r")
          flush(file_out)
          GC.gc()
        end
      end
    end

    if continue_evo == 0
      Ene_H0 =  obs_Ene0
      Ene_H_evo = obs_Ene_evo
      S_site=( obs_sz, obs_sy, obs_sx)
      DW_C=( obs_Cz, obs_Cy, obs_Cx)


    elseif continue_evo == 1
      Ene_H0=vcat(Ene_H0, obs_Ene0)
      Ene_H_evo=vcat(Ene_H_evo, obs_Ene_evo)
      S_site=(vcat(S_site[1],obs_sz),vcat(S_site[2], obs_sy),vcat(S_site[3], obs_sx))
      DW_C = (vcat(DW_C[1], obs_Cz),vcat(DW_C[2], obs_Cy),vcat(DW_C[3], obs_Cx))

    end
elseif time_evo_method == "TEBD_Ht" 
    obs_Ene0 = []
    obs_p = []
    obs_sz = []
    obs_sy = []
    obs_sx = []
    obs_Cz = []
    obs_Cy = []
    obs_Cx = []
    if continue_evo == 1
      t0=t_end
      write(file_out, "\r!!! Start time evolution from t0=$t_end !!!")
      write(file_out, "\r")
    elseif continue_evo == 0
      t0=0.0
      write(file_out, "\r!!! Start time evolution from t0=0 !!!")
      write(file_out, "\r")
    else
    end

    if Lattice_type == 1 
      if Mobile_DW == 1
        global evo_gates = evo_gates_TEBD_LT_1_Ht(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz,dhx,dhy,omega,t0,tau)
      else
        global NBC=[1,N]
        global evo_gates = evo_gates_TEBD_LT_1_Ht(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz,dhx,dhy,omega,t0,tau)
      end
    elseif Lattice_type == 5
     
      global evo_gates =evo_gates_TEBD_LT_5_Ht(N,sites,hx,hy,hz,dhx,dhy,omega,t0,tau)
    else
    end

    global psi_temp = apply(evo_gates, psi_init; cutoff)
      
    for (i, t) in enumerate(tau:tau:t_total)

      if Lattice_type == 1 
        if Mobile_DW == 1
          global evo_gates = evo_gates_TEBD_LT_1_Ht(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz,dhx,dhy,omega,t0+t,tau)
        else
          NBC=[1,N]
          global evo_gates = evo_gates_TEBD_LT_1_Ht(N,sites,NBC[1],NBC[2],J,Kz,Ky,hx,hy,hz,dhx,dhy,omega,t0+t,tau)
        end
      elseif Lattice_type == 5

        global evo_gates =evo_gates_TEBD_LT_5_Ht(N,sites,hx,hy,hz,dhx,dhy,omega,t0+t,tau)
      else
      end
      push!( obs_Ene0,real(inner(psi_temp', H_evo,psi_temp)))
      push!( obs_p,abs(inner(psi_temp,psi[band_tar]))^2)
      push!( obs_sz,expect(psi_temp, "Sz"))
      push!( obs_sy,expect(psi_temp, "Sy"))
      push!( obs_sx,expect(psi_temp, "Sx"))
      push!( obs_Cz,inner(psi_temp',DWzMPO,psi_temp))
      push!( obs_Cy,inner(psi_temp',DWyMPO,psi_temp))
      push!( obs_Cx,inner(psi_temp',DWxMPO,psi_temp))
      if write_psi_evo == 1
      global psi_evo=vcat(psi_evo, psi_temp)
      end
      
      t≈t_total && break
        
      global psi_temp = apply(evo_gates, psi_temp; cutoff)
      normalize!(psi_temp)

      if Int(round(t_total/tau/100)) ==0
        push!(DATE,Dates.DateTime(Dates.now()))
        local rightnow=DATE[end]
        local E0_print=round(obs_Ene0[i],digits=8)

        write(file_out, "\rTime step: $t/$t_total  sweep:$i  Ene0 = $E0_print   Date: $rightnow")
        write(file_out, "\r")
        flush(file_out)
        GC.gc()
      else
        if i % Int(round(t_total/tau/100)) == 0 
            push!(DATE,Dates.DateTime(Dates.now()))
            local rightnow=DATE[end]
            local E0_print=round(obs_Ene0[i],digits=8)
    
            write(file_out, "\rTime step: $t/$t_total  sweep:$i  Ene0 = $E0_print   Date: $rightnow")
            write(file_out, "\r")
            flush(file_out)
            GC.gc()
        end
      end
    end

    if continue_evo == 0
      Ene_H0 =  obs_Ene0
      p01 = obs_p
      S_site=( obs_sz, obs_sy, obs_sx)
      DW_C=( obs_Cz, obs_Cy, obs_Cx)


    elseif continue_evo == 1
      Ene_H0=vcat(Ene_H0, obs_Ene0)
      p01=vcat(Ene_H_evo, obs_p)
      S_site=(vcat(S_site[1],obs_sz),vcat(S_site[2], obs_sy),vcat(S_site[3], obs_sx))
      DW_C = (vcat(DW_C[1], obs_Cz),vcat(DW_C[2], obs_Cy),vcat(DW_C[3], obs_Cx))

    end

elseif time_evo_method == "TDVP" || time_evo_method == "TDVP_Im_time"
    if time_evo_method == "TDVP"
      write(file_out, "\r!!! Start TDVP Calculation (real time) !!!")
      write(file_out, "\r")
    elseif time_evo_method == "TDVP_Im_time"
      write(file_out, "\r!!! Start TDVP Calculation (imaginary time)!!!")
      write(file_out, "\r")
    else
    end
    ################# measures of TDVP loops ################
    step(; sweep) = sweep
    current_time(; current_time) = current_time
    return_state(; state) = state
    measure_Ene0(; state) = real(inner(state', H, state))
    measure_Ene_evo(; state) = real(inner(state', H_evo, state))
    measure_sz(; state) = expect(state, "Sz")
    measure_sy(; state) = expect(state, "Sy")
    measure_sx(; state) = expect(state, "Sx")
    measure_Cz(; state) = inner(state',DWzMPO,state)
    measure_Cy(; state) = inner(state',DWyMPO,state)
    measure_Cx(; state) = inner(state',DWxMPO,state)
    timer(; sweep, current_time, state) = begin
        if sweep % Int(t_total/tau/100) == 0 && sweep != 0
          push!(DATE,Dates.DateTime(Dates.now()))
          local rightnow=DATE[end]
          println(file_out,"Time step: ", round(current_time,digits=2),"/",t_total,"   Ene0 = ", round(real(inner(state', H, state)),digits=8),"   Ene_time = ", round(real(inner(state', H_evo, state)),digits=8),"   Date: ",rightnow)
          write(file_out, "\r")
          flush(file_out)
        end
      return nothing
    end

    obs = observer(
      "steps" => step, "times" => current_time, "states" => return_state, "Ene0" => measure_Ene0, "Ene_evo" => measure_Ene_evo,
        "sz" => measure_sz, "sy" => measure_sy, "sx" => measure_sx,
        "Cz" => measure_Cz, "Cy" => measure_Cy, "Cx" => measure_Cx,"timer" => timer
    )
    ################# measures of TDVP loops ################

    if time_evo_method == "TDVP"
      psi_temp = tdvp(-im*H_evo, t_total, psi_init; time_step=tau, cutoff, (step_observer!)=obs, outputlevel=0)
    elseif time_evo_method == "TDVP_Im_time"
      psi_temp = tdvp(-H_evo, t_total, psi_init; time_step=tau, cutoff, (step_observer!)=obs, outputlevel=0)
    else
    end
    
    Ene_H0 = obs.Ene0
    Ene_H_evo = obs.Ene_evo
    S_site=(obs.sz,obs.sy,obs.sx)
    DW_C=(obs.Cz,obs.Cy,obs.Cx)
    if write_psi_evo == 1
     psi_evo=obs.states
    end

elseif time_evo_method == "TDVP_Ht"
    write(file_out, "\r!!! Start TDVP Calculation with time dependent Hamiltonian !!!")
    write(file_out, "\r")
    
    if continue_evo == 1
      t0=t_end
      write(file_out, "\r!!! Start time evolution from t0=$t_end !!!")
      write(file_out, "\r")
    elseif continue_evo == 0
      t0=0.0
      write(file_out, "\r!!! Start time evolution from t0=0 !!!")
      write(file_out, "\r")
    else
    end
    
    H_evo_total=Ham_tot_TDVP(N,sites,H_evo,dhx,dhy,omega,t0)

    ################# measures of TDVP loops ################
      step(; sweep) = sweep
      current_time(; current_time) = current_time
      if write_psi_evo == 1
        return_state(; state) = state
      end
      measure_Ene0(; state) = real(inner(state',H_evo, state))
      measure_p(; state) = abs(inner(state,psi[band_tar]))^2
      measure_sz(; state) = expect(state, "Sz")
      measure_sy(; state) = expect(state, "Sy")
      measure_sx(; state) = expect(state, "Sx")
      measure_Cz(; state) = inner(state',DWzMPO,state)
      measure_Cy(; state) = inner(state',DWyMPO,state)
      measure_Cx(; state) = inner(state',DWxMPO,state)
      timer(; sweep, current_time, state) = begin
        if Int(round(t_total/tau/100)) ==0
          push!(DATE,Dates.DateTime(Dates.now()))
          local rightnow=DATE[end]
          println(file_out,"Time step: ", round(current_time,digits=2),"/",t_total,"   Ene0 = ", round(real(inner(state', H_evo, state)),digits=8),"   Date: ",rightnow)
          # 内存使用信息
          current_memory = Sys.total_memory() - Sys.free_memory()
          println(file_out, "   Current memory usage: ", current_memory / (1024^3), " GB")
          # 打印单个变量（如 `state`）的内存占用
          state_memory = Base.summarysize(state)
          println(file_out, "   State size: ", state_memory / (1024^2), " MB")
          write(file_out, "\r")
          flush(file_out)
          #GC.gc()
        else
          if sweep % Int(round(t_total/tau/100)) == 0 && sweep != 0
            push!(DATE,Dates.DateTime(Dates.now()))
            local rightnow=DATE[end]
            println(file_out,"Time step: ", round(current_time,digits=2),"/",t_total,"   Ene0 = ", round(real(inner(state', H_evo, state)),digits=8),"   Date: ",rightnow)
            # 内存使用信息
            current_memory = Sys.total_memory() - Sys.free_memory()
            println(file_out, "Current memory usage: ", round(current_memory / (1024^3),digits=8), " GB")
            # 打印单个变量（如 `state`）的内存占用
            state_memory = Base.summarysize(state)
            println(file_out, "State size: ", round(state_memory / (1024^2),digits=8), " MB")
            write(file_out, "\r")
            flush(file_out)
            #GC.gc()
          end
        end
        return nothing
      end
      if write_psi_evo == 1
        obs = observer(
          "steps" => step, "times" => current_time, "states" => return_state, "Ene0" => measure_Ene0, 
           "sz" => measure_sz, "sy" => measure_sy, "sx" => measure_sx,
           "Cz" => measure_Cz, "Cy" => measure_Cy, "Cx" => measure_Cx,"timer" => timer, "p" => measure_p
        )    
      else
        obs = observer(
          "steps" => step, "times" => current_time, "Ene0" => measure_Ene0, 
           "sz" => measure_sz, "sy" => measure_sy, "sx" => measure_sx,
           "Cz" => measure_Cz, "Cy" => measure_Cy, "Cx" => measure_Cx,"timer" => timer, "p" => measure_p
        )     
      end

    ################# measures of TDVP loops ################

    psi_temp = tdvp( -im*H_evo_total,t_total,psi_init;updater=krylov_updater,updater_kwargs=(; tol=converg, eager=true),time_step=tau,cutoff,nsite, (step_observer!)=obs,outputlevel=0)
      
    if continue_evo == 0
      Ene_H0 = obs.Ene0
      S_site=(obs.sz,obs.sy,obs.sx)
      DW_C=(obs.Cz,obs.Cy,obs.Cx)
      p01=obs.p
      if write_psi_evo == 1
        psi_evo=obs.states
      end
    elseif continue_evo == 1
      Ene_H0=vcat(Ene_H0, obs.Ene0)
      p01=vcat(p01,obs.p)
      S_site=(vcat(S_site[1],obs.sz),vcat(S_site[2],obs.sy),vcat(S_site[3],obs.sx))
      DW_C = (vcat(DW_C[1],obs.Cz),vcat(DW_C[2],obs.Cy),vcat(DW_C[3],obs.Cx))

      if write_psi_evo == 1
        psi_evo=vcat(psi_evo,obs.states)
      end
    end

else
  write(file_out, "\rTime_evo_method is not assigned.")
  push!(DATE,Dates.Time(Dates.now()))
  rightnow=DATE[end]
  write(file_out, "\rDate: $rightnow")
  write(file_out, "\r")
  exit()
end

write(file_out, "\rTime step: $t_total/$t_total")
push!(DATE,Dates.DateTime(Dates.now()))
rightnow=DATE[end]
write(file_out, "\rDate: $rightnow")
write(file_out, "\r")
write(file_out, "\rFinished time evolution of state $band_evo")
flush(file_out)
GC.gc()

###### Saving psi_evo Data #######
file_psi = h5open(string("psi_evo_end_",band_evo,".h5"),"w")
  write(file_psi,"psi", psi_temp)
  close(file_psi)
  file_psi = h5open(string("psi_evo_start_",band_evo,".h5"),"w")
  write(file_psi,"psi",psi_init)
  close(file_psi)

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

if time_evo_method == "TEBD"
  jldsave("time_evo_data.jld2"; Ene_H0,Ene_H_evo,S_site,DW_C)
 elseif time_evo_method == "TDVP" || time_evo_method == "TDVP_Im_time"
  jldsave("time_evo_data.jld2"; Ene_H0,Ene_H_evo,S_site,DW_C)
 elseif time_evo_method == "TDVP_Ht" || time_evo_method == "TEBD_Ht"
  jldsave("time_evo_data.jld2"; Ene_H0,S_site,DW_C,p01)
 end

####### timer  ##################
write(file_out, "\rSimulation Finished.")
push!(DATE,Dates.DateTime(Dates.now()))
rightnow=DATE[end]
runtime=round(Dates.value(DATE[end]-DATE[1])/1E3/60;digits = 2)
write(file_out, "\rDate: $rightnow")
write(file_out, "\rtotal runtime: $runtime mins")
write(file_out, "\r###############################################")
close(file_out)
