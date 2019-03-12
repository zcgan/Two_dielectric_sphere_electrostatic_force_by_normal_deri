%This matlab script calculates the force between two charged dielectric bodies from
%the surface normal derivative on grids, governed by Poisson equantion, and solved by BEM.

%%%%%%%%%%about input files%%%%%%%%%%%%%%%%%%%%%%%%
%All we need is the location of each collocation point, the normal derivative on each collocation point,
%and the weights for each point (I am using CC, so the weights are the triangle areas, but for other schemes, the force can be calculated as well).
%these information are saved in data.dat

%further more, now I am consiering the two sphere case, so for simplicity,
%consider two point charges located at each sphere center. Then we also
%need the location of each sphere center, the source charge magnitude, and
%the dielectric constant inside and outside each sphere, and the number of boundary elements for each sphere,
%these information are saved in source.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=load('data.dat');
%load the collocation pts locations, normal derivative value, and weights from the data file
grid_loc=data(:,2:4); %same as tr_xyz(1:3,i) in tabi
normal_dev=data(:,8); %same as xvct(i), the surface phi norm derv
weights=data(:,10); %same as tr_area(i) in tabi

normal_dev=normal_dev/4/pi; %rescale the unit of weights

%load the source/dielectric object information
sourcedata=load('source.dat');
src_loc=sourcedata(:,2:4); %location of the central charges, same as chrpos(1:3,i)
src_char=sourcedata(:,5); %magnitude of the central charges, same as atmchr(i)
eps_i=sourcedata(1,6); %now we assume each sphere has the same dielectric constant inside (same as eps0)
eps_s=sourcedata(1,7); %eps1
N=sourcedata(:,8); %the number of triangles on each sphere %same as nface/N, here N=2

index=zeros(2,2); %record the start and end index for each sphere
index(1,1)=1; index(1,2)=N(1); %nface1
index(2,1)=N(1)+1;%nface1+1 
index(2,2)=N(1)+N(2);%nface


%obtaining the induced surface charge density using Eq.(1.4) in the writeup.
sigma_b=(1-eps_i/eps_s).*normal_dev;

%Calculate the force using Eq.(2.5) in the writeup.
%initialize the force vectors in xyz directions
Fx=zeros(1,2); Fy=zeros(1,2); Fz=zeros(1,2);

%Step1: calculate the first contribution, i.e., the source-source pairwise Colomb
%interaction
for i=1:2 %two bodies
    for j=1:2
        if j~=i
            bem_ind_start=1+(j-1)*N(i);
            bem_ind_end=j*N(j);
            dx=src_loc(i,1)-src_loc(j,1);
            dy=src_loc(i,2)-src_loc(j,2);
            dz=src_loc(i,3)-src_loc(j,3);
            dr=sqrt(dx*dx+dy*dy+dz*dz);
            
            Fx(i)=Fx(i)+src_char(i)*src_char(j)*dx/eps_i/eps_i/dr/dr/dr;
            Fy(i)=Fy(i)+src_char(i)*src_char(j)*dy/eps_i/eps_i/dr/dr/dr;
            Fz(i)=Fz(i)+src_char(i)*src_char(j)*dz/eps_i/eps_i/dr/dr/dr;
        end
    end
end
            

%Step2: calculate the second contribution, the force on the source charges
%due to the induced surface charges.
for i=1:2 %two bodies
    for j=1:2
        if j~=i
            for k=index(j,1):index(j,2)
            dx=src_loc(i,1)-grid_loc(k,1);
            dy=src_loc(i,2)-grid_loc(k,2);
            dz=src_loc(i,3)-grid_loc(k,3);
            dr=sqrt(dx*dx+dy*dy+dz*dz);
            
            Fx(i)=Fx(i)+src_char(i)*sigma_b(k)*weights(k)*dx/eps_i/dr/dr/dr;
            Fy(i)=Fy(i)+src_char(i)*sigma_b(k)*weights(k)*dy/eps_i/dr/dr/dr;
            Fz(i)=Fz(i)+src_char(i)*sigma_b(k)*weights(k)*dz/eps_i/dr/dr/dr;
            end
        end
    end
end


%Step3: calculate the third contribution, the force on the induced charges
%due to the source charge (maybe symmetric with step2?)
for i=1:2 %two bodies
    for j=1:2
        if j~=i
            for k=index(i,1):index(i,2)
            dx=grid_loc(k,1)-src_loc(j,1);
            dy=grid_loc(k,2)-src_loc(j,2);
            dz=grid_loc(k,3)-src_loc(j,3);
            dr=sqrt(dx*dx+dy*dy+dz*dz);
            
            Fx(i)=Fx(i)+src_char(j)*sigma_b(k)*weights(k)*dx/eps_i/dr/dr/dr;
            Fy(i)=Fy(i)+src_char(j)*sigma_b(k)*weights(k)*dy/eps_i/dr/dr/dr;
            Fz(i)=Fz(i)+src_char(j)*sigma_b(k)*weights(k)*dz/eps_i/dr/dr/dr;
            end
        end
    end
end

%Step4: calculate the fourth contribution, the force on the induced charges
%due to the other induced charges
for i=1:2 %two bodies
    for j=1:2
        if j~=i
            for k=index(i,1):index(i,2)
                for l=index(j,1):index(j,2)
                    
            dx=grid_loc(k,1)-grid_loc(l,1);
            dy=grid_loc(k,2)-grid_loc(l,2);
            dz=grid_loc(k,3)-grid_loc(l,3);
            dr=sqrt(dx*dx+dy*dy+dz*dz);
            
            Fx(i)=Fx(i)+sigma_b(l)*weights(l)*sigma_b(k)*weights(k)*dx/dr/dr/dr;
            Fy(i)=Fy(i)+sigma_b(l)*weights(l)*sigma_b(k)*weights(k)*dy/dr/dr/dr;
            Fz(i)=Fz(i)+sigma_b(l)*weights(l)*sigma_b(k)*weights(k)*dz/dr/dr/dr;
                end
            end
        end
    end
end
Fx=Fx*eps_s; 
Fy=Fy*eps_s; 
Fz=Fz*eps_s; 

disp('The result is:')
disp('sphere#1   sphere#2')
disp('Fx='); disp(Fx)
disp('Fy='); disp(Fy)
disp('Fz='); disp(Fz)