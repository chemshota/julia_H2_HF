using SpecialFunctions #for function erf
using LinearAlgebra

#def Fgamma
Fgamma(x) =1/2*(pi/x)^0.5*erf(x)
#def kinetic_integral
function kinetic_integral(a,b,Ra,Rb)
    if a == 0 && b == 0
        return 0.0
    else
        R_AB_sqr = norm(Ra .- Rb)^2
        temp1=(a*b  / (a+b)) 
        temp2 = (3-2*a*b / (a + b)*R_AB_sqr)
        temp3 = (pi/(a+b))^1.5
        temp4 = exp(-temp1*R_AB_sqr)
        temp1 * temp2 * temp3 * temp4
    end
end
#def coulomb integral
function el_cor_coulomb_integral(a,b,Ra,Rb,Rc,Zc)
    if a == 0 && b == 0
        return 0.0
    else
        Rp=(a*Ra+b*Rb)/(a+b)
        if norm(Rp .- Rc)==0
            return 0.0
        else
            R_AB_sqr = norm(Ra .- Rb)^2
            temp1 = (-2*pi*Zc/(a+b))
            temp2 = exp( - a*b/(a+b)*R_AB_sqr)
            temp3 = Fgamma((a+b)* norm(Rp .- Rc)^2)
            temp1 * temp2 * temp3
        end
    end
end

function calc_kinetic_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom)
    coef=norm_const .* contr_coef
    nrow=size(orbital_expornent)[1]
    ncol=size(orbital_expornent)[2]
    T_tensor=zeros(3,4,3,4)
    for i in 1:nrow, j in 1:nrow, k in 1:ncol,l in 1:ncol
        T_tensor[i,k,j,l]=kinetic_integral(orbital_expornent[i,k],orbital_expornent[j,l],orbital_senter_geom[k],orbital_senter_geom[l])
    end
    kinetic=zeros(4,4)
    for i in 1:ncol, j in 1:ncol
        for k in 1:nrow , l in 1:nrow
        kinetic[i,j] += T_tensor[k,i,l,j]*coef[k,i]*coef[l,j]
        end
    end 
    kinetic
end

function calc_el_cor_coulomb_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom,xyz,Zcharge)
    coef=norm_const .* contr_coef
    nrow=size(orbital_expornent)[1]
    ncol=size(orbital_expornent)[2]
    V_tensor=zeros(3,4,3,4)
    for i in 1:nrow, j in 1:nrow, k in 1:ncol,l in 1:ncol
        for c in 1:length(Zcharge)
            V_tensor[i,k,j,l]+=el_cor_coulomb_integral(orbital_expornent[i,k],orbital_expornent[j,l],orbital_senter_geom[k],orbital_senter_geom[l],xyz[:,c],Zcharge[c])
        end
    end
    coulomb=zeros(4,4)
    for i in 1:ncol, j in 1:ncol
        for k in 1:nrow , l in 1:nrow
        coulomb[i,j] += V_tensor[k,i,l,j]*coef[k,i]*coef[l,j]
        end
    end 
    coulomb
end




#initialize variables
num_of_atom = 2
num_of_elec = num_of_atom
SCFsuccess = false
Nb = 4
Nb2 = Nb*(Nb+1)/2
Nb4 = Nb2*(Nb2+1)/2
lwork=15
nucleus_repulsion = 1.0/1.322808
Zcharge=[1,1]
length(Zcharge)
overlap=zeros(4,4)
kinetic=zeros(4,4)
nucattr=zeros(4,4)

#set geometry
xyz=[0.0 0.0 ; 0.0 0.0  ; 0.661404  -0.661404]

#set orbital centor
orbital_senter=[1,1,2,2]
orbital_senter_geom=[xyz[:,i] for i in orbital_senter]




#set orbital_expornent
orbital_expornent=zeros(3,4)
H2_valcence_3_expornent = [18.731137 , 2.8253997 , 0.6401217]
H2_valcence_1_expornent = [0.1612778 , 0.0 , 0.0]
orbital_expornent[:,[1,3]] .=  H2_valcence_3_expornent
orbital_expornent[:,[2,4]] .=  H2_valcence_1_expornent


#set contr_coef
contr_coef=zeros(3,4)
H2_valcence_3_conrtr_coef = [0.03349460,0.23472695,0.81375733]
H2_valcence_1_conrtr_coef = [1.0 , 0.0 , 0.0]
contr_coef[:,[1,3]] .=  H2_valcence_3_conrtr_coef
contr_coef[:,[2,4]] .=  H2_valcence_1_conrtr_coef


#calc norm_const
norm_const=(2 .+ orbital_expornent ./ pi).^0.75

#calc over_lap


T_mat=calc_kinetic_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom)
V_mat=calc_el_cor_coulomb_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom,xyz,Zcharge)
H_one = T_mat+V_mat
