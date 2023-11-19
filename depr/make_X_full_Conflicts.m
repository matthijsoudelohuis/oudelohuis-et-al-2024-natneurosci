function [X_full,params] = make_X_full_Conflicts(params,tempsessionData, temptrialData, tempvideoData)

%get relevant dimensions:
params.nTrials                     = length(temptrialData.session_ID);
params.nTotalTimebins              = params.nTrials * params.nTimebins;

%Zscore the center position of the pupil in the video:
tempvideoData.center_x{1}           = zscore(tempvideoData.center_x{1});
tempvideoData.center_y{1}           = zscore(tempvideoData.center_y{1});

tempvideoData.zarea{1}(isnan(tempvideoData.zarea{1})) = 0;

if params.smoothSVD
    for iSVD = 1:params.nSVDs
        tempvideoData.motSVD{1}(:,iSVD) = smooth(tempvideoData.motSVD{1}(:,iSVD),params.smoothSVD);
    end
end
tempvideoData.motSVD{1}             = zscore(tempvideoData.motSVD{1},[],1);

%Create predictor matrix w/ kernels:
switch params.modelVersion
    case 1
        [dm,output.x_label]                              = MOL_molGLM_createModel(params,tempsessionData,temptrialData,tempvideoData); %see function for details
    case 2
        [dm,output.x_label]                              = MOL_molGLM_createModel_v2(params,tempsessionData,temptrialData,tempvideoData); %see function for details
    case 3
        %Create predictor matrix w/ kernels:
        [dspec,dm]                              = MOL_neuroGLM_createModel(params,tempsessionData,temptrialData,tempvideoData); %see function for details
        if sum([dspec.covar.edim]) ~= size(dm.X,2); error('Dimensions dont match'); end
        output.x_label = {}; %Store the label that belongs to each predictor (because of temp basis functions, multiple predictors belong to single task vars)
        for iC = 1:length(dspec.covar)
            output.x_label = [output.x_label; repmat({dspec.covar(iC).desc},dspec.covar(iC).edim,1)];
        end
end


X_full = dm.X;

end