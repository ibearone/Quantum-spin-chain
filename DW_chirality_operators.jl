#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code defining the         #
#chirality operator Cz,Cx,Cy                #
############################################# 

########## DW chirality operators  ############
function DWC_operator_1D(N::Int,sites)
    global DWz= OpSum()
    for i =1:N-1
        DWz += "Sx",i,"Sy",i+1
        DWz += -1,"Sy",i,"Sx",i+1
    end


    global DWzMPO=MPO(DWz,sites)

    global DWy= OpSum()
    for i =1:N-1
        DWy += "Sz",i,"Sx",i+1
        DWy += -1,"Sx",i,"Sz",i+1
    end


    global DWyMPO=MPO(DWy,sites)

    global DWx= OpSum()
    for i =1:N-1
        DWx += "Sy",i,"Sz",i+1
        DWx += -1,"Sz",i,"Sy",i+1
    end

    global DWxMPO=MPO(DWx,sites)

return DWxMPO,DWyMPO,DWzMPO
end


function DWC_operator_2D(Nx::Int,Ny::Int,sites)
    global DWz1= OpSum()
    for i = 1:Nx-1
        DWz1 += "Sx",i,"Sy",i+1
        DWz1 += -1,"Sy",i,"Sx",i+1
    end

    global DWzMPO1=MPO(DWz1,sites)

    global DWz2= OpSum()
    for i = Nx+1:2*Nx -1
        DWz2 += "Sx",i,"Sy",i+1
        DWz2 += -1,"Sy",i,"Sx",i+1
    end

    global DWzMPO2=MPO(DWz2,sites)
    global DWy1= OpSum()
    for i =1:Nx-1
        DWy1 += "Sz",i,"Sx",i+1
        DWy1 += -1,"Sx",i,"Sz",i+1
    end

    global DWy2= OpSum()
    for i =Nx+1:2*Nx-1
        DWy2 += "Sz",i,"Sx",i+1
        DWy2 += -1,"Sx",i,"Sz",i+1
    end



    global DWyMPO1=MPO(DWy1,sites)
    global DWyMPO2=MPO(DWy2,sites)

    global DWx1= OpSum()
    for i =1:Nx-1
        DWx1 += "Sy",i,"Sz",i+1
        DWx1 += -1,"Sz",i,"Sy",i+1
    end

    global DWx2= OpSum()
    for i =Nx+1:2*Nx-1
        DWx2 += "Sy",i,"Sz",i+1
        DWx2 += -1,"Sz",i,"Sy",i+1
    end


    global DWxMPO1=MPO(DWx1,sites)
    global DWxMPO2=MPO(DWx2,sites)

return DWxMPO1,DWyMPO1,DWzMPO1,DWxMPO2,DWyMPO2,DWzMPO2
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