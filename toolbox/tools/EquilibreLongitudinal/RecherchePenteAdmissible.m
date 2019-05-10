function [valid,PenteEq] = RecherchePenteAdmissible(PenteCom,ThetaCible,param,data)
% Recherche de la pente admissible la plus proche de PenteCom
% par balayage sur tous les points de la courbe

data.AmplitudeTheta = 2*ThetaCible ;
data.n = max(ceil(2*ThetaCible/0.1*max([0.01;abs(2*data.c2*data.AmplitudeTheta)+abs(data.c1);abs(2*data.d2*data.AmplitudeTheta)+abs(data.d1)])),3) ;

%ThetaTest = data.ThetaMinAppx + ((1:(data.n))-1)/((data.n)-1)*data.AmplitudeTheta ;
ThetaTest = -ThetaCible + ((1:(data.n))-1)/((data.n)-1)*2*ThetaCible ;
Beta1Test = data.c2*ThetaTest.*ThetaTest + data.c1*ThetaTest + data.cpi*PenteCom + data.c0t ;
Beta2Test = data.d2*ThetaTest.*ThetaTest + data.d1*ThetaTest + data.dpi*PenteCom + data.d0t ;


% Recueil du besoin de variation de pente pour tous les points de la
% trajectoire et selection du point avec variation de pente minimale
dp_mini = 1e10 ;
for ii=1:(data.n)
    [intersection,delta_pente] = IncrementPente(Beta1Test(ii),Beta2Test(ii),param,data) ;
    if ((intersection)&&(abs(delta_pente) < abs(dp_mini)))
        dp_mini = delta_pente ;
    end
end

if (1e10 == dp_mini)
    % aucun moyen d'équilibrer le sous-marin
    valid = 0 ;
    PenteEq = 0 ;
else
    valid = 1 ;
    PenteEq = PenteCom + dp_mini ;
end


Beta1TestCorr= data.c2*ThetaTest.*ThetaTest + data.c1*ThetaTest + data.cpi*PenteEq + data.c0t ;
Beta2TestCorr = data.d2*ThetaTest.*ThetaTest + data.d1*ThetaTest + data.dpi*PenteEq + data.d0t ;

% figure(1)
% plot(Beta1Test*57.3,Beta2Test*57.3,'r-',Beta1TestCorr*57.3,Beta2TestCorr*57.3,'b-',[param.Beta1Min param.Beta1Max,param.Beta1Max,param.Beta1Min,param.Beta1Min]*57.3,[param.Beta2Min,param.Beta2Min,param.Beta2Max,param.Beta2Max,param.Beta2Min]*57.3,'k--')
% grid on
