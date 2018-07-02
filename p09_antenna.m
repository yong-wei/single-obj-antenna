function [f,F,dx] = p09_antenna(x,mpopt)
% Antenna objectives
% Parameter.Mode ------------ 'A'       %Amplitude 
%                             'A_phi'   %Amplitude + Phase
%                             'A_phi_p' %Amplitude + Phase + Position
% Parameter.position_mode --- 'asymmetry'                 
%                             'symmetry'
%                             'coordinate_coded'
% Example 
%           load('data');
%           F = fitness_antenna(x);
% ---------------------------------------------------------
% Pan Anqi panqi516@gmail.com 29-6-2018
% ---------------------------------------------------------
% Modified by Yong-Wei Zhang 02-07-2018, yongwzhang@gmail.com
% x: [0,1];
% Amplitude: [0,1]
% Phase: [-pi/4,pi/4];
% distance: [0.5*lambda,lambda]
% f: single objective
% F: multi objective for each indivadual
% dx: design parameter decoded from x
[nPop,nVar] = size(x);
% if nargin < 2
    numda = 3e8/1.3e9; % wavelength    
    theta = 0:1/10:90; % angle resolution=1/1.5
    Parameter = struct('Mode', 'A_phi_p'...
                      ,'position_mode', 'symmetry'...
                      ,'N', 28 ...
                      ,'NullType', 'deep'...
                      ,'numda', numda ...    % wavelength    
                      ,'d', 0.5*numda ...
                      ,'k', (2*pi)/numda ... % number of waves
                      ,'theta',[fliplr(-theta(2:end)),theta] ...% angle resolution=1/1.5
                       );                         
    Parameter.ObjIdx = [1 2 5]; % chose objectives
% end
Parameter.draw = 0;
%--------------------------------------------------------------------------
% This part is only for initialize parameters of algorithms
if exist('mpopt','var') %
    if strcmp(mpopt,'initial')
        f.ID = 'p09';
        f.FunctionName = 'Antenna Design';
        f.FeasibleBounds = {zeros(1,nVar),ones(1,nVar)};
        f.FeasibleDimension = 6:3:42;
        f.GlobalMinima =  -101.3973;
        f.BestSoFarFitness =  -101.3973;
        return
    elseif strcmp(mpopt,'ITE')
        f = [para.LB',para.UB'];
        return
    elseif strcmp(mpopt,'draw')
        Parameter.draw = 1;
    end
end
%--------------------------------------------------------------------------

F = zeros(nPop,8);
% M = length(Parameter.ObjIdx); 
Ns = length(Parameter.theta);          % Sample points
broadnull = [50,60];
dx = zeros(nPop,Parameter.N*3);
for n = 1:nPop
Pos = zeros(1,Parameter.N);
switch Parameter.Mode
    case {'A'}
        switch Parameter.position_mode
            case 'asymmetry'
                Amp = x(n,1:Parameter.N);% Asymmetry amp
            case 'symmetry'
                Amp = [flip(x(n,1:Parameter.N/2)),x(n,1:Parameter.N/2)];% symmetry amp
        end 
        i = 1:1:Parameter.N;
        Pos = (i-0.5)*Parameter.d;% Unify distribued position
        Phi = zeros(1,Parameter.N);% Same phase
    case {'A_phi'}
        switch Parameter.position_mode
            case 'asymmetry'
                Amp = x(n,1:Parameter.N);% Asymmetry amp
                Phi = x(n,Parameter.N+1:2*Parameter.N);% Asymmetry phase               
            case 'symmetry'
                Amp = [flip(x(n,1:Parameter.N/2)),x(n,1:Parameter.N/2)];% symmetry amp
                Phi = [flip(x(n,Parameter.N/2+1:Parameter.N)),x(n,Parameter.N/2+1:Parameter.N)];
            case  'a_coordinate_coded' % Asymmetry
                Amp = x(n,1:Parameter.N);% Asymmetry amp
                Phi = Parameter.numda*x(n,Parameter.N+1:2*Parameter.N);
        end
        i = 1:1:Parameter.N;
        Pos = (i-0.5)*Parameter.d;% Unify distribued position
    case {'A_phi_p'}        
        switch Parameter.position_mode
            case 'asymmetry' % asymmetry position
                Amp = x(n,1:Parameter.N);% Asymmetry amp
                Phi = x(n,Parameter.N+1:2*Parameter.N);% Asymmetry phase
                distance = x(n,2*Parameter.N+1:3*Parameter.N);
                Pos = cumsum(distance);
            case 'symmetry' % symmetry position
                Amp = [flip(x(n,1:Parameter.N/2)),x(n,1:Parameter.N/2)];% symmetry amp
                Phi = [flip(x(n,Parameter.N/2+1:Parameter.N)),x(n,Parameter.N/2+1:Parameter.N)];
                distance = x(n,Parameter.N+1:3*Parameter.N/2);% symmetry coordinate
                distance(2:end) = mapX(distance(2:end),[0.5*numda,numda]);
                distance(1) = mapX(distance(1),[0.25*numda,0.5*numda]);
                Pos = [-flip(cumsum(distance)),cumsum(distance)];
            case   'a_coordinate_coded' % coordinate encode
                Amp = x(n,1:Parameter.N);% asymmetry amp
                Phi = x(n,Parameter.N+1:2*Parameter.N);% asymmetry phase
                Pos = x(n,2*Parameter.N+1:3*Parameter.N); 
        end  
end
Phi = mapX(Phi,[-pi/4,pi/4]);
dx(n,:) = [Amp,Phi,Pos];
f1 = sum(repmat(Amp',1,Ns).*cos(Parameter.k*sin(repmat(Parameter.theta,Parameter.N,1).*pi/180).*...
    repmat(Pos',1,Ns)+repmat(Phi',1,Ns)),1);
f_out = 20*log10(abs(f1)/max(abs(f1)));% transferom to Db
if isfield(Parameter,'draw')
    if Parameter.draw == 1
        plot(Parameter.theta,f_out);
        ylim([-60,0]);
    end
end
%% FNBW--first null beam width
fMax = find(f_out==max(f_out));
for l = fMax+1:Ns
    if(f_out(l) > f_out(l-1))
        Parameter.Theta_main = Parameter.theta(l-1);
        break;
    end
end
F(n,5) = Parameter.Theta_main*2;%abs(Parameter.Theta_main*2-dFNBW);
%% SLL
SideLobe = [f_out(1:find(Parameter.theta>=-Parameter.Theta_main,1)),f_out(find(Parameter.theta>=Parameter.Theta_main,1):Ns)];
F(n,1) = max(SideLobe);% Only points other than the main lobe, i.e. the side lobes.
%% NULL
if ~isfield(Parameter,'NullType')
    Parameter.NullType = 'deep';
end
switch Parameter.NullType
    case 'deep'
        % deep nulls
        deepNoo = [find(Parameter.theta>=30,1),find(Parameter.theta>=32.5,1),find(Parameter.theta>=35,1),find(Parameter.theta>=-30,1),find(Parameter.theta>=-32.5,1),find(Parameter.theta>=-35,1)];
        deepNoo = f_out(floor(deepNoo));
        F(n,2) = mean(deepNoo);%mean(abs(deepNoo-deep));
    case 'board'
        % broad nulls
        deepNzoo = [find(Parameter.theta>=-broadnull(2),1),find(Parameter.theta>=-broadnull(1),1); ...
            find(Parameter.theta>=broadnull(1),1),find(Parameter.theta>=broadnull(2),1)];
        % deep zero area, avg. power = zero
        broad1 = f_out(floor(deepNzoo(1,1)):floor(deepNzoo(1,2)));
        broad2 = f_out(floor(deepNzoo(2,1)):floor(deepNzoo(2,2)));
        F(n,2) = 0.5*(mean(broad1)+mean(broad2));
end
%% Average mainlobe
MainLobe = f_out(find(Parameter.theta>=-Parameter.Theta_main,1):find(Parameter.theta>=Parameter.Theta_main,1));% only the main lobe points
F(n,3) = -mean(MainLobe);
%% HPBW--half power beam width
HPBW = find(f_out>=-3);
if numel(HPBW) == 0
    HPBW = 0;
end
F(n,4) = -Parameter.theta(HPBW(1))+Parameter.theta(numel(HPBW));
%% Directivity ????
Pall = sum((abs(f1)).^2.*cos(Parameter.theta.*pi/180)*pi/Ns);
F(n,6) = -10*log10(2*max(abs(f1))^2/Pall);
%% Aperture
F(n,7) = (Pos(Parameter.N)-Pos(1))/Parameter.numda;
%% Shaped array
SideLeft = f_out(1:find(Parameter.theta>=-Parameter.Theta_main,1));
SideRight = f_out(find(Parameter.theta>=Parameter.Theta_main,1):Ns);
F(n,8) =  (sum(abs(SideLeft+20)) + sum(abs(SideRight+30)))/Ns;
%% Combine
end
F = F(:,Parameter.ObjIdx);
f = sum(F,2);
f(isnan(f)) = inf;

function x = mapX(x,b)
x = x*(b(2)-b(1))+b(1);