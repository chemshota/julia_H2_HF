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
function el_at_coulomb_integral(a,b,Ra,Rb,Rc,Zc)
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

function calc_el_at_coulomb_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom,xyz,Zcharge)
    coef=norm_const .* contr_coef
    nrow=size(orbital_expornent)[1]
    ncol=size(orbital_expornent)[2]
    V_tensor=zeros(3,4,3,4)
    for i in 1:nrow, j in 1:nrow, k in 1:ncol,l in 1:ncol
        for c in 1:length(Zcharge)
            V_tensor[i,k,j,l]+=el_at_coulomb_integral(orbital_expornent[i,k],orbital_expornent[j,l],orbital_senter_geom[k],orbital_senter_geom[l],xyz[:,c],Zcharge[c])
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

function two_electron_integral(a,b,c,d,Ra,Rb,Rc,Rd)
    if (a==0 && b==0) || (c==0 && d==0)
        return 0.0
    else
        Rp=(a*Ra+b*Rb)/(a+b)
        Rq=(c*Rc+d*Rd)/(c+d)
        if Rp == Rq
            return 0.0
        else
            temp1 = 2*pi^2.5/((a+b)*(c+d)*(a+b+c+d)^0.5)
            temp2 = exp(-a*b/(a+b)*norm(Ra-Rb)^2)
            temp3 = exp(-c*d/(c+d)*norm(Rc-Rd)^2)
            temp4 = Fgamma((a+b)*(c+d)/(a+b+c+d)*norm(Rp-Rq)^2)
            temp1 * temp2 * temp3 * temp4
        end
    end
end

function calc_two_electron_integral(norm_const,contr_coef,Nb,orbital_senter_geom,orbital_expornent)
    two_el_integral = zeros(4,4,4,4)
    coef=norm_const .* contr_coef
    nrow=size(orbital_expornent)[1]
    for i in 1:Nb, j in 1:Nb, k in 1:Nb, l in 1:Nb
        (Ra,Rb,Rc,Rd)=(orbital_senter_geom[i],orbital_senter_geom[j],orbital_senter_geom[k],orbital_senter_geom[l])
        two_el_integral[i,j,k,l] = sum([ coef[ia,i] * coef[ib,j] * coef[ic,k] * coef[id,l] * two_electron_integral(orbital_expornent[ia,i],orbital_expornent[ib,j],orbital_expornent[ic,k],orbital_expornent[id,l],Ra,Rb,Rc,Rd) for ia in 1:nrow  for ib in 1:nrow  for ic in 1:nrow  for id in 1:nrow ])
    end
    two_el_integral
end

function calc_J_K_matrix(Nb,two_el_integral,Pmatrix) 
    Jmatrix=zeros(4,4)
    Kmatrix=zeros(4,4)
    for i in 1:Nb, j in 1:Nb
        Jmatrix[i,j] = sum([two_el_integral[i,j,xci,theta]*Pmatrix[xci,theta] for xci in 1:Nb for theta in 1:Nb ])
        Kmatrix[i,j] = - 1/2 * sum([two_el_integral[i,theta,xci,j]*Pmatrix[xci,theta] for xci in 1:Nb for theta in 1:Nb ])
    end
    (Jmatrix,Kmatrix)
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

#one_eletron_matrix
T_mat=calc_kinetic_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom)
V_mat=calc_el_at_coulomb_mat(norm_const,contr_coef,orbital_expornent,orbital_senter_geom,xyz,Zcharge)
H_core = T_mat+V_mat

#two_eletron_matrix
two_el_integral=calc_two_electron_integral(norm_const,contr_coef,Nb,orbital_senter_geom,orbital_expornent)

#init density
Pmatrix=diagm(0 => [0.5,0.5,0.5,0.5])
(Jmatrix,Kmatrix)=calc_J_K_matrix(Nb,two_el_integral,Pmatrix) 
Fock_matrix=H_core+Jmatrix-Kmatrix
new_energy=dot(Nb*Nb*(H_core+Fock_matrix),Pmatrix)*0.5+nucleus_repulsion
