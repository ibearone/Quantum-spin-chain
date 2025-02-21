#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the total   #
# spin operator Sz,Sy,Sx                    #
############################################# 

########## Total Spin Operators  ############
function Tot_spin_Op(sites)
   
    global os_Sz= OpSum()
    for i =1:N
        global os_Sz += 1,"Sz",i
    end
    SztotMPO=MPO(os_Sz,sites);

    global os_Sy= OpSum()
    for i =1:N
        global os_Sy += 1,"Sy",i
    end
    SytotMPO=MPO(os_Sy,sites);

    global os_Sx= OpSum()
    for i =1:N
        global os_Sx += 1,"Sx",i
    end
    SxtotMPO=MPO(os_Sx,sites);
return SztotMPO,SytotMPO,SxtotMPO

end