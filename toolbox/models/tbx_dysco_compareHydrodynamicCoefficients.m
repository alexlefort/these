function tbx_dysco_compareHydrodynamicCoefficients(...
                                           navData,scenarioResults,options)
% TBX_DYSCO_COMPAREHYDRODYNAMICCOEFFICIENTS takes several
% navData inputs and plots the hydrodynamic coefficients wrt simulation
% data.
%
% Inputs:
%   - navData           : Structure containing manoeuvrability model
%   - scenarioResults   : Cell array with scenario results, represented
%                         with filename or loaded data structure
%   - options           : [Optional] Structure with optional fields
%                          - tstart         : Start time. Default is zero
%                          - tend           : End time. Default is inf
%
% Example
%
% tbx_dysco_compareHydrodynamicCoefficients(...
%   {'ScorpeneChili1.nav','ScorpeneChili2.nav'},...
%   {'DP8C-1-Simulations.Res','DP8C-2-Simulations.Res'});
%
%
% SIREHNA
% GJ
%==========================================================================
% SVN info
% SVN $Id$
% SVN $HeadURL$
%==========================================================================
optionsDef.tstart           = 0;
optionsDef.tend             = Inf;
optionsDef.color            = [];
if nargin < 3
    options = optionsDef;
    if nargin < 2
        scenarioResults = 'uuv.Res';
        if nargin == 0
            navData = 'Fx212vA.nav';
        end
    end
end
options = tbx_struct_addMissingFields(optionsDef,options);
if ischar(navData)
    C = tbx_dysco_readNavData(navData);
elseif iscell(navData)
    for i=1:numel(navData)
        if ischar(navData{i})
            c = tbx_dysco_readNavData(navData{i});
            if i==1
                C = c;
            else
                C = tbx_struct_concatenate(C, c);
            end
        elseif isstruct(navData{i})
            if i==1
                C = navData{i};
            else
                C = tbx_struct_concatenate(C, navData{i});
            end
        end
    end
elseif isstruct(navData)
    C = navData;
end
if ischar(scenarioResults)
    scenarioResults = tbx_dysco_getDyscoResults(scenarioResults,true);
elseif iscell(scenarioResults)
    for i=1:numel(scenarioResults)
        if ischar(scenarioResults{i})
            scenarioResults{i} = tbx_dysco_getDyscoResults(scenarioResults{i}, true);
        end
    end
end
if ~iscell(scenarioResults)
    scenarioResults = {scenarioResults};
end
nScenario = numel(scenarioResults);
if numel(options.tstart)==1 && nScenario>1
    options.tstart = repmat(options.tstart,1,nScenario);
end
if numel(options.tend)==1 && nScenario>1
    options.tend = repmat(options.tend,1,nScenario);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = tbx_dysco_checkNavData(C);
for j=1:nScenario
    id = scenarioResults{j}.t>=options.tstart(j) & scenarioResults{j}.t<=options.tend(j);
    if ~all(id)
        f = fieldnames(scenarioResults{j});
        for i=1:numel(f)
            scenarioResults{j}.(f{i}) = scenarioResults{j}.(f{i})(id);
        end
    elseif ~any(id)
        error('No time elements to consider');
    end
end
H=zeros(1,6);
coeffsNames = {'Cx','Cy','Cz','Cl','Cm','Cn'};
for i=1:6
    H(i) = figure('Name', coeffsNames{i}); hold on; grid on; box on;
    ylabel(coeffsNames{i});
end
cc = 'krgb';
for i=1:numel(C)
    [resConcatenate,offset,t] = tbx_dysco_concatenateHydrodynamicCoefficients(C(i),scenarioResults);

    figure(H(1));
    plot(t,sum(resConcatenate.Cx.values,2),'color',cc(i));
    figure(H(2));
    plot(t,sum(resConcatenate.Cy.values,2),'color',cc(i));
    figure(H(3));
    plot(t,sum(resConcatenate.Cz.values,2),'color',cc(i));
    figure(H(4));
    plot(t,sum(resConcatenate.Cl.values,2),'color',cc(i));
    figure(H(5));
    plot(t,sum(resConcatenate.Cm.values,2),'color',cc(i));
    figure(H(6));
    plot(t,sum(resConcatenate.Cn.values,2),'color',cc(i));
end
scenariiNames = cell(1,nScenario);
for j = 1:nScenario
    if isfield(scenarioResults{j}, 'scenarioName')
        scenariiNames{j} = scenarioResults{j}.scenarioName;
    else
        scenariiNames{j} = '';
    end
end
for i=1:numel(H)
    figure(H(i));
    hold on;
    rr = axis;
    dx = rr(2)-rr(1);
    dy = rr(4)-rr(3);
    axis([rr(1) rr(2)+0.1*dx rr(3)-0.2*dy rr(4)]);
    rr = axis;
    dy = rr(4)-rr(3);
    v = [offset(:,3);offset(end,3)+scenarioResults{end}.t(end)];
    for j=1:nScenario+1
        plot(v(j)*[1 1],rr(3:4),'k-','LineWidth',2);
    end
    for j=1:nScenario
        text(0.5*(v(j)+v(j+1)),rr(3)+0.1*dy,...
            {char(64+j);scenariiNames{j}},...
            'HorizontalAlignment','center',...
            'Interpreter','none',...
            'Rotation',90,...
            'BackgroundColor',[1.0 1.0 1.0]);
    end
    set(gca,'XGrid','off');
end



%   C = tbx_dysco_readNavData('ScorpeneChili.nav', true);
%   scenarioResults1 = tbx_dysco_getDyscoResults('DP8C-1-Simulations.Res', true);
%   scenarioResults2 = tbx_dysco_getDyscoResults('DP8C-2-Simulations.Res', true);
%   [res,offset,t] = tbx_dysco_concatenateHydrodynamicCoefficients(C, {scenarioResults1, scenarioResults2});
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,h] = plotResult(uvwpqr, res, options)
optionsDef.plotOnlyNonZeroResults = true;
optionsDef.titleHeader = '';
optionsDef.color = [];
if nargin < 3
    options = optionsDef;
end
options = tbx_struct_addMissingFields(optionsDef,options);
if isstruct(uvwpqr)
    t = uvwpqr.t;
else
    t = uvwpqr;
end
f = fieldnames(res);
nf = numel(f);
h = cell(1,nf);
H = zeros(1,nf);
for i=1:nf
    nC = size(res.(f{i}).values,2);
    if options.plotOnlyNonZeroResults
        plotId = max(abs(res.(f{i}).values))>0;
    else
        plotId = true(1,nC);
    end
    plotInd = [1 1+find(plotId)];
    if any(plotId)
        H(i) = figure('Name',[options.titleHeader ' ' f{i}]);
        box on;
        h{i} = plot(t,[sum(res.(f{i}).values,2),res.(f{i}).values(:,plotId)]);
        if isempty(options.color)
            % color = tbx_figure_color(plotInd,'lines');
            v = 1:(nC+1);
            v = v(plotInd);
            color = [[1,0,0];tbx_doe_sobol_i4_sobol_generate(3,nC,666)'];
            color = color(v,:);
        else
            color = options.color;
        end
        for j=1:nC
            set(h{i}(j),'color', color(j,:));
        end
        legend(h{i},[f(i) res.(f{i}).coeff(plotId)]);
    end
    grid on;
    xlabel('Time (s)');
    title(f{i});
end
