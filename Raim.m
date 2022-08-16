clear all

Satellite_S = importdata('SatellitePostion.txt');
Receiver = importdata('ReceiverPostion.txt');
%���ջ����ݣ�1,2,3-XYZ,4-λ�ø��£�5,6,7-PreciseXYZ��8-λ�ø���
Satellite = Satellite_S(:,[2,3,17,18,19,25,27,28,29,34]);
%�������ݣ�1-TOW,2- PRN,3,4,5-XYZ,6-RANGE,7,8,9- PreciseXYZ,10PreciseRANGE

BugSate1 = find(Satellite(:, 2) == 3); %����3����Ϊ������
len1 = length(BugSate1);
Satellite(BugSate1,7) = Satellite(BugSate1,7) + randi([20,30],len1,1);
Satellite(BugSate1,8) = Satellite(BugSate1,8) + randi([20,30],len1,1);
Satellite(BugSate1,9) = Satellite(BugSate1,9) + randi([20,30],len1,1);

% BugSate2 = find(Satellite(:, 2) == 13); %����13����Ϊ������
% len2 = length(BugSate2);
% Satellite(BugSate2,7) = Satellite(BugSate2,7) + randi([-20,-10],len2,1);
% Satellite(BugSate2,8) = Satellite(BugSate2,8) + randi([-20,-10],len2,1);
% Satellite(BugSate2,9) = Satellite(BugSate2,9) + randi([-20,-10],len2,1);

TimeSpan = unique(Satellite(:, 1));     %%SOW
N = [];
N_Satellite = [];

td = [4.7989, 4.8352, 4.8658, 4.8922, 4.9153, 4.9359, 4.9545, 4.9713, 4.9868, 5.0011, 5.0144, 5.0268, 5.0384, 5.0494, 5.0597, 5.0695];
sigma_T = [4.4172, 4.7985, 5.0894, 5.3360, 5.5548, 5.7539, 5.9379, 6.1100, 6.2722, 6.4262, 6.5731, 6.7139, 6.8492, 6.9797, 7.1058, 7.2281];
for ii = 1:size(TimeSpan)
    this_TOW = TimeSpan(ii);
    index = find(Satellite(:, 1) == this_TOW);
    CurrData = Satellite(index, :);  %ȡ��Ӧʱ�������й۲�����
    nsat = size(CurrData,1);
    sat = CurrData(:,2); %�ɹ۲����Ǻ� 
    Pse = CurrData(:,10); %α��
    if nsat < 5 
        continue        
    end
            
    H = [];        
    R = [];
    for jj = 1:nsat
        R(jj) = sqrt( (Receiver(ii,5) - CurrData(jj,7))^2 + (Receiver(ii,6) - CurrData(jj,8))^2 + (Receiver(ii,7) - CurrData(jj,9))^2);
        H(jj,1:3) = (Receiver(ii,5:7) - CurrData(jj,7:9))/R(jj);
        H(jj,4) = 1;
    end    
    H'
    %��ȡ5����
        PRN = sat(1:5);
        R1 = R(1:5);
        H1 = H(1:5,:);
        S1 = eye(5) - H1 * inv(H1'*H1) * H1';
%         if (nsat-rank(S)-4) == 0   %ֻȡ������Ҫ�������
%             V = R' - CurrData(:,10);
%             err = S*V;
%             del = sqrt(err'*err/(nsat-4));  %���ͳ����
%         else
%             continue
%         end   
          V1 = R1' - Pse(1:5); %α�����
          err = S1*V1;
          Sigma = sqrt((err'* err)/(5-4));
%         del = sqrt(err'*err/(nsat-4));  %���ͳ����
%         N = [N;[nsat,nsat-rank(S)]];  %S�����Ƿ�Ϊnsat-4
          %������Ϊ0.5%�����λ����T^2=7.8794
          Sigma0 = 2; 
          
          %% add
          SigmaT = Sigma0 * sigma_T(1);
           
          
          %T = sqrt(7.8794);
          %SigmaT = Sigma0*T;%ֻȡ5���ǣ�n-4Ϊ1
          %if Sigma < SigmaT  %ѡȡ��5���ǿ���
          if Sigma < SigmaT
              R2 = R(1:5);
              H2 = H(1:5,:);
              Pse2 = Pse(1:5);
              PRN = sat(1:5);
              for k0 = 1:nsat-5
                    R2 = [R2,R(5+k0)];
                    H2 = [H2;H(5+k0,:)];
                    Pse2 = [Pse2;Pse(5+k0)];
                    S = eye(size(H2,1)) - H2 * inv(H2'*H2) * H2';
                    V = R2' - Pse2;
                    err = S*V;
                    Sigma = sqrt(err'*err/(length(PRN)-3));
                    if Sigma < SigmaT
                        PRN = [PRN;sat(5+k0)];  %ֻ��ȡ����Ӧ�����Ǻ�
                    else
                        R2(end) = [];
                        H2(end,:) = [];
                        Pse2(end) = [];
                    end
              end     
          else      %5�����д���һ�Ź�����
              d = abs(V1)./(Sigma0 * sqrt(diag(S1)));  %���߼���б�ʣ�ѡб������
              [Max,passPRN] = max(d);
%               A = inv(H1'*H1)*H1'; 
%               Slope = sqrt((A(2,:).^2 + A(3,:).^2)'./diag(S1));
%               [Max,passPRN] = max(Slope);
    

              
              %%lambda = (S1 * V1 * V1') / (Sigma0 * Sigma0); 
              R(passPRN) = [];   %�ѹ����Ƕ�Ӧ���������
              H(passPRN,:) = [];
              Pse(passPRN) = [];
              sat(passPRN) = [];
              R3 = R(1:4);
              H3 = H(1:4,:);
              Pse3 = Pse(1:4);
              PRN = sat(1:4);
              for t0 = 1:nsat-5
                  R3 = [R3,R(4+t0)];
                  H3 = [H3;H(4+t0,:)];
                  Pse3 = [Pse3;Pse(4+t0)];
                  S = eye(size(H3,1)) - H3 * inv(H3'*H3) * H3';
                  V = R3' - Pse3;%%
                  err = S*V;
                  Sigma = sqrt(err'*err/(length(PRN)-3));
                  if Sigma < SigmaT
                        PRN = [PRN;sat(4+t0)];  %ֻ��ȡ����Ӧ�����Ǻ�
                  else
                      R3(end) = [];      
                      H3(end,:) = [];
                      Pse3(end) = [];
                  end
               end
          end
      N_Satellite = [N_Satellite;PRN];    
end

for tt = 1:nsat-5
  D(tt) = abs(err(tt)) / (Sigma0 * sqrt(S1(tt,tt)));
  if D(tt) > td(5-4)
      tt
  end
end




