function Param = ay_creat_model()
Param = ay_create_state_space(9,4,0,9,eye(9,9),[],[],[1 2 3 4 5 6 7 8 9],[0 0 0 0 0 0 0 0 0]);
RATE_A = 0.1;
RATE_B = 0.1;
RATE_C = 0.1;
% set initial value for Wv,s1 to RATE_A
Param.Ak(1,1)=  1-RATE_A; 
Param.Ak(2,1)=  1;          Param.Ak(2,2)= -RATE_B;
% set initial value for Wv,s2 to RATE_A
Param.Ak(3,3)=  1-RATE_A;
Param.Ak(4,3)=  1;          Param.Ak(4,4) = -RATE_B;
% set initial value for Wv,s3 to RATE_A
Param.Ak(5,5)= 1-RATE_A;
Param.Ak(6,5)= 1;           Param.Ak(6,6)= - RATE_B;
% set initial value for Wv,s4 to RATE_A
Param.Ak(7,7)= 1-RATE_A;
Param.Ak(8,7)= 1;           Param.Ak(8,8)=  - RATE_B;
% set initial value for Wq,a1 to 1
Param.Ak(9,9)= 1-RATE_C;
% set initial value of B
Param.Bk(1,1)   = RATE_A;
Param.Bk(2,1)   = RATE_B;
Param.Bk(3,2)   = RATE_A;
Param.Bk(4,2)   = RATE_B;
Param.Bk(5,3)   = RATE_A;
Param.Bk(6,3)   = RATE_B;
Param.Bk(7,4)   = RATE_A;
Param.Bk(8,4)   = RATE_A;
Param.Bk(9,1:4) = RATE_C;
% set Param.Wk
Param.Wk = 0.01 * eye(9,9);
% Set Param.Ek
Param.Ek(1) = 0;Param.Ek(3) = 0;Param.Ek(5) = 0;Param.Ek(7) = 0;
end