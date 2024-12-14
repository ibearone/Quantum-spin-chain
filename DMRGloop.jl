#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code of single DMRG loop. #
############################################# 

########## DMRG loop  ############
mutable struct CheckObserver <: AbstractObserver
    energy_tol::Float64
    last_energy::Float64
    
    
    global CheckObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

function ITensors.checkdone!(o::CheckObserver;kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    global energy_diff = abs(energy-o.last_energy)
    if energy_diff < o.energy_tol
      println("Stopping DMRG after sweep $sw with energy difference $energy_diff.")
      return true
    end
    # Otherwise, update last_energy and keep going
    o.last_energy = energy
    return false
end


function ITensors.measure!(o::CheckObserver; kwargs...)
      energy = kwargs[:energy]
      sweep = kwargs[:sweep]
      bond = kwargs[:bond]
      half_sweep = kwargs[:half_sweep]
      psi = kwargs[:psi]
      projected_operator = kwargs[:projected_operator]
  if bond==1 && half_sweep==2
    psi_size =  Base.format_bytes(Base.summarysize(psi))
    PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
    println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
  end
end



function DMRG_loop(H,psi_previous,psi_init,nsweeps,maxdim,mindim,cutoff,noise,sweep_num_loop,E_loop,Control,etol)

    #eigsolve_verbosity = 1::Int
    obs = CheckObserver(etol)
   

    if Control ==0
        global energy_temp,psi_temp =@time dmrg(H,psi_init; nsweeps, maxdim,cutoff,noise,mindim, observer=obs,outputlevel=1,eigsolve_krylovdim=5);
        push!(sweep_num_loop,nsweeps)
        push!(E_loop,energy_temp)
        #local sweep_show=sweep_num_loop[end]
        
        
        #if  length(E_loop) > 1
            #local D_E=abs(E_loop[end]-E_loop[end-1])
            local E_variance = real(inner(H,psi_temp,H,psi_temp) -  inner(psi_temp',H,psi_temp)^2)
            println("Convergence of energy : D_E= $energy_diff")
            println("Variance of energy : Delta_E^2= $E_variance")
            
        #end
        #println("##### End of Sweep: $sweep_show  #######")
        #push!(Sz_tot_loop,real(inner(psi_temp',SztotMPO,psi_temp)))  
    elseif Control ==1
        global energy_temp,psi_temp =@time dmrg(H,psi_previous,psi_init; nsweeps, maxdim,cutoff,noise,mindim,observer=obs,outputlevel=1,eigsolve_krylovdim=5);
        push!(sweep_num_loop,nsweeps)
        push!(E_loop,energy_temp)
        #local sweep_show=sweep_num_loop[end]
        
        
        #if  length(E_loop) > 1
            #local D_E=abs(E_loop[end]-E_loop[end-1])
            local E_variance = real(inner(H,psi_temp,H,psi_temp) -  inner(psi_temp',H,psi_temp)^2)
            println("Convergence of energy : D_E= $energy_diff")
            println("Variance of energy : Delta_E^2= $E_variance")
        #end
        #println("##### End of Sweep: $sweep_show  #######")
        
        #push!(Sz_tot_loop,real(inner(psi_temp',SztotMPO,psi_temp)))
        #push!(D_E_loop,abs(E_loop[end]-E_loop[end-1]))


    end
return energy_temp,psi_temp,sweep_num_loop,E_loop
end

