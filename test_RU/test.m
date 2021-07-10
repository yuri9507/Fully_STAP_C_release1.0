clear all;
load ('Ru.mat');

fileID = fopen("Ru_for_matlab.bin","rb");
Ru_from_c = fread(fileID,'double');
fclose(fileID);

Ru_from_c_test = reshape(Ru_from_c,2,409600);
Real_OF_RU_FROM_c = Ru_from_c_test(1,:);
Image_OF_RU_FROM_c = Ru_from_c_test(2,:);
Complex_OF_RU_FROM_c = complex(Real_OF_RU_FROM_c,Image_OF_RU_FROM_c);
RU_TEST = zeros(640,640);

for i = 0:1:639
    for j = 1:1:640
        %RU_TEST(1,1) = Complex_OF_RU_FROM_c(1,0*640+1)
        %RU_TEST(1,2) = Complex_OF_RU_FROM_c(1,0*640+2)...RU_TEST(1,640) =Complex_OF_RU_FROM_c(1,0*640+640)
        %RU_TEST(2,1) = Complex_OF_RU_FROM_c(1,1*640+1) ..RU_TEST(2,640) = Complex_OF_RU_FROM_c(1,1*640+640)
        RU_TEST(i+1,j) = Complex_OF_RU_FROM_c(1,i*640+j);
        
    end
    
end


fileID_1 = fopen("Ru_inverse_for_matlab.bin","rb");
Ru_inverse_from_c = fread(fileID_1,'double');
fclose(fileID_1);


Ru_inverse_from_c_test = reshape(Ru_inverse_from_c,2,409600);
Real_OF_RU_INVERSE_FROM_c = Ru_inverse_from_c_test(1,:);
Image_OF_RU_INVERSE_FROM_c = Ru_inverse_from_c_test(2,:);
Complex_OF_RU_INVERSE_FROM_c = complex(Real_OF_RU_INVERSE_FROM_c,Image_OF_RU_INVERSE_FROM_c);
RU_INVERSE_TEST = zeros(640,640);

for i = 0:1:639
    for j = 1:1:640
        
        RU_INVERSE_TEST(i+1,j) = Complex_OF_RU_INVERSE_FROM_c(1,i*640+j);
        
    end
    
end

identity_matrix= eye(640,640);
A = RU_TEST*RU_INVERSE_TEST - identity_matrix;

Ru_inverse = inv(Ru);

B = Ru * Ru_inverse - identity_matrix; % 1.5972 e-05
%B = Ru_inverse * Ru - identity_matrix;   %1.5210 e-08
%B = Ru \ Ru - identity_matrix;   %1.4079e-08
%B = Ru_inverse \ Ru_inverse - identity_matrix;   %2.9323e-08



n_from_C_Ru_C_Ru_INV = norm(A,'fro');
n_from_Matalb_Ru_Matalb_Ru_inv = norm(B,'fro');

C = RU_TEST * Ru_inverse - identity_matrix;
n_from_C_Ru_Matalb_Ru_inv = norm(C,'fro');

D = Ru * RU_INVERSE_TEST - identity_matrix;
n_from_Matlab_Ru_C_Ru_inv = norm(D,'fro');

E = RU_INVERSE_TEST*RU_TEST - identity_matrix;
n_from_C_Ru_inv_C_Ru = norm(E,'fro');

F = Ru_inverse * Ru - identity_matrix;
n_from_Matlab_Ru_inv_Matlab_Ru = norm(F,'fro');


