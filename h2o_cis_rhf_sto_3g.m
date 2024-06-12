%%%
%
function [] = h2o_cis_rhf_sto_3g
clc; 
%
format short
%
[En_total ,f, tei_mol] = h2o_rhf_sto3g;

[En_total] % -74.965900481678958 (atomic unit, au), ground state energy
%
[f]; % = [                               
    %-20.251556277131314    % orbital energy
    % -1.257624813925194
    % -0.593933613384772
    % -0.459749506519034
    % -0.392624438360491
    %  0.581943686802284
    %  0.692851608959985];
%
hole = [1, 2, 3, 4, 5];   % hole state
virt = [6, 7];            % virtial (particle) state   
%
%
dim1B = length(f);
A_singlet = zeros(dim1B,dim1B);
A_triplet = zeros(dim1B,dim1B);
%B = zeros(dim1B,dim1B);
II = 0;
for i = hole 
    for a = virt
        II = II + 1;
        JJ = 0;
        for j = hole 
            for b = virt 
                JJ = JJ + 1;
                A_singlet(II,JJ) = ((f(a) - f(i)) * (i==j) * (a==b)) + (2*tei_mol(i,a,b,j) - tei_mol(i,j,a,b)); % singlet singly excited states
                A_triplet(II,JJ) = ((f(a) - f(i)) * (i==j) * (a==b)) + (0*tei_mol(i,a,j,b) - tei_mol(i,j,a,b)); % triplet singly excited states
%
%                B(II,JJ) = tei_data(a,b,i,j);
            end
        end
    end
end
%%%
one_au = 27.211324570273;% eV
%
[eig(A_singlet), eig(A_singlet) * one_au]; % CIS singlet excitaion energies in (au), (eV)
%
1.0e+02 *[  
   0.200866949787128   5.465855766098269
   0.201306908906967   5.477827636505857
   0.010183403834917   0.277103906982081
   0.014425005604119   0.392523509421686
   0.006049171710963   0.164605974808331
   0.014582104106088   0.396798367748276
   0.004590036449483   0.124900971616276
   0.005152251084047   0.140199576515549
   0.007665618424153   0.208591630971489
   0.006722278401474   0.182922099434239]; %  
%
[eig(A_triplet), eig(A_triplet) * one_au]; % CIS triplet excitaion energies in (au), (eV)

1.0e+02 * [
   0.200260539145265   5.449354529299660
   0.200896615895569   5.466663020203799
   0.012438451133452   0.338466730943852
   0.006840949516749   0.186151297669107
   0.013437335146605   0.365647688033810
   0.004694586860706   0.127745926790006
   0.003838438423941   0.104448993796863
   0.004668813546587   0.127044600774279
   0.006269423429664   0.170599315813060
   0.005207567681806   0.141704814411292];

%%%
return
end

%%%
function [En_total ,f, tei_mol] = h2o_rhf_sto3g
%
format long
%
S_ov = [ 
  1.00000143e+00  2.50986878e-01  5.00165930e-02  4.54010988e-01 0.00000000e+00  2.92760956e-01 -2.45538547e-01;
  2.50986878e-01  1.00000143e+00  5.00165930e-02  4.54010988e-01 0.00000000e+00 -2.92760956e-01 -2.45538547e-01;
  5.00165930e-02  5.00165930e-02  1.00000143e+00  2.36703475e-01 0.00000000e+00  0.00000000e+00  2.51146260e-17;
  4.54010988e-01  4.54010988e-01  2.36703475e-01  9.99999786e-01 0.00000000e+00  0.00000000e+00  5.14887930e-17;
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 1.00000019e+00  0.00000000e+00  0.00000000e+00;
  2.92760956e-01 -2.92760956e-01  0.00000000e+00  0.00000000e+00 0.00000000e+00  1.00000019e+00  0.00000000e+00;
 -2.45538547e-01 -2.45538547e-01  2.51146260e-17  5.14887930e-17 0.00000000e+00  0.00000000e+00  1.00000019e+00
];
%
H_core = [
 -4.95657708e+00 -1.56021450e+00 -1.61969970e+00 -3.54335553e+00 0.00000000e+00 -1.89101453e+00  1.65871752e+00;
 -1.56021450e+00 -4.95657708e+00 -1.61969970e+00 -3.54335553e+00 0.00000000e+00  1.89101453e+00  1.65871752e+00;
 -1.61969970e+00 -1.61969970e+00 -3.26851555e+01 -7.60432563e+00 0.00000000e+00  0.00000000e+00  1.86791769e-02;
 -3.54335553e+00 -3.54335553e+00 -7.60432563e+00 -9.30208434e+00 0.00000000e+00  0.00000000e+00  2.22152591e-01;
  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 -7.43085558e+00  0.00000000e+00  0.00000000e+00;
 -1.89101453e+00  1.89101453e+00  0.00000000e+00  0.00000000e+00 0.00000000e+00 -7.56706412e+00  0.00000000e+00;
  1.65871752e+00  1.65871752e+00  1.86791769e-02  2.22152591e-01 0.00000000e+00  0.00000000e+00 -7.52666703e+00];

%%%
%
dim = 7; % size of basis sets & 1s x 2 + 2s1p (2x1 + 1x3 = 5) = 2 + 5 = 7.;
%
%
N_el = 10.;               % number of electron 
%
itermax = 100; tol = 1e-8;
%
tei_n = 2401;             % = 7^4, .i.e., all values of TEI
%
read_tei_data = fopen('tei_h2o_sto_3g.txt', 'r');               % data of two-electron integral in atomic basis set
tei_data_n5 = textscan(read_tei_data, '%d %d %d %d %f');
%
p = zeros(tei_n,1); q = zeros(tei_n,1); r = zeros(tei_n,1); s = zeros(tei_n,1); vals = zeros(tei_n,1);
p(1:tei_n) = tei_data_n5{1};
q(1:tei_n) = tei_data_n5{2};
r(1:tei_n) = tei_data_n5{3};
s(1:tei_n) = tei_data_n5{4};
vals(1:tei_n) = tei_data_n5{5};
for i = 1:tei_n
    tei(p(i),q(i),r(i),s(i)) = vals(i);
%    tei(q(i),p(i),r(i),s(i)) = vals(i);    
%    tei(p(i),q(i),s(i),r(i)) = vals(i);    
%    tei(q(i),p(i),s(i),r(i)) = vals(i);   
    %
%    tei(r(i),s(i),p(i),q(i)) = vals(i);    
%    tei(s(i),r(i),p(i),q(i)) = vals(i);        
%    tei(r(i),s(i),q(i),p(i)) = vals(i);        
%    tei(s(i),r(i),q(i),p(i)) = vals(i);            
end
tei;
%
P_old = 0.5 * ones(dim,dim); % initial charge population
for iter = 1:itermax
    iter;
    P = P_old;
    %
    F = H_core;
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F(p,q) = F(p,q) + P(r,s) * (tei(p,q,r,s) - 0.5.*tei(p,r,q,s));
                end
    
            end
    
        end
    end
    Ham_fock = F ;     % Fock matrix
    S_mat_fock = S_ov;

    [Vec,En] = eig(Ham_fock,S_mat_fock);                                     % Eigenvalue problem: F*c = En*S*c - Roothaan equation
    En = diag(En);
    [foo, ij] = sort(En);
    En = En(ij);
    [En(1), En(2)];  % orbital energies
    %
    Vec = Vec(:,ij);                       % expansion coefficients 
    %
    for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec(p,i) * Vec(q,i) * S_ov(p,q);
            end
        end
        Vec(:,i) = Vec(:,i)/sqrt(norm);
    end
    %
    P_new = zeros(dim,dim);
    for i = 1:N_el/2
        for pp = 1:dim
            for qq = 1:dim
                P_new(pp,qq) = P_new(pp,qq) + 2*Vec(pp,i)*Vec(qq,i);
            end
        end
    end
    %
     if (abs(sum(sum(P_new-P_old))) < tol)
            break 
     end
    %        
    P_old = P_new;

end
%%%
En_0 = (sum(0.5*diag(P(:,:)*(H_core(:,:) + F(:,:))))); % -83.873808940969425 

%%%
h1_xyz = [0.,  1.43233673, -0.96104039];  % geometry for H2O molecule I used
h2_xyz = [0., -1.43233673, -0.96104039];
o_xyz  = [0.,  0.,          0.24026010];
%%%
z1 = 1.; % nuclear charge for H atom
z2 = 1.; % nuclear charge for H atom
z3 = 8.; % nuclear charge for O atom
%
V_nuc_rep = z1*z2/sqrt(sum((h1_xyz - h2_xyz).^2)) + ...
            z1*z3/sqrt(sum((h1_xyz - o_xyz).^2)) + ...
            z2*z3/sqrt(sum((h2_xyz - o_xyz).^2)); % 8.907908459290468
%

En_total = En_0 + V_nuc_rep;  % -74.965900481678958, total energy for ground state in atomic unit (au)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
P = Vec'; 
%
tei_mol = zeros(dim,dim,dim,dim); % two-electron integral in molecular basis set
for ii = 1:dim
    for jj = 1:dim
        for kk = 1:dim
            for ll = 1:dim
                for mm = 1:dim
                    for nn = 1:dim
                        for oo = 1:dim
                            for pp = 1:dim
                                tei_mol(ii,jj,kk,ll) =  tei_mol(ii,jj,kk,ll) + P(ii,mm)*P(jj,nn)*P(kk,oo)*P(ll,pp)*tei(mm,nn,oo,pp); 
                            end
                        end
                    end
                end
            end
        end
    end
end
%%%%%

f = En;

%%%
return
end
