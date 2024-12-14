#############################################
# Quantum Spin Chain Simulation by Julia    #
# Author: Guanxiong Qu                      #
# Date: 2024.12.14                          #
# Version: v2                               #
# Email: quguanxiong@gmail.com              #
# This is the sub code for reading inputs   #
# from file "Input".                        #
############################################# 

########## Module for Reading Inputs  ############

function read_input(file_name::IOStream,Data_name::String,Data_Type,Control)
    data_found = false
# Initialize  variables 
    ### Read Int,Float64 ####
    if Control ==0
        for line in eachline(file_name)
            if contains(line, Data_name)
                global Data = parse(Data_Type, strip(split(line, "=")[2]))
                data_found = true
                break
            end

        end
        if !data_found
            error("'$Data_name' not found in the file.")
            # You can handle this case further, for example, by setting a default value for Data
        end
        seekstart(file_name)

    ### Read array of Int Float64 ####
    elseif Control == 1
        for line in eachline(file_name)
            if contains(line, Data_name)
                 Data_dum=split(split(line,"=")[2],",")
                 global Data = [parse(Data_Type,Data_dum[i] ) for i in eachindex(Data_dum)]
                 data_found = true
                break
            end
        end
        if !data_found
            error("$Data_name not found in the file.")
           
        end
        seekstart(file_name)
    ### Read String  ####
    elseif Control == 2
        for line in eachline(file_name)
            if contains(line, Data_name)
                global Data = string(strip(split(line, "=")[2]))
                data_found = true
                break
            end
        end
        if !data_found
            error("'$Data_name' not found in the file.")
          
        end
        seekstart(file_name)
    ### Read others ####
    else
    error(" Wrong control in read_input of '$Data_name'. ")

        
    end
return Data
end