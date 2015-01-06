rad_models_ml = zeros(4,1);
rad_models_ml(1) = computeMarginal('posteriorModel_1_RAD.mat');
rad_models_ml(2) = computeMarginal('posteriorModel_2_RAD.mat');
rad_models_ml(3) = computeMarginal('posteriorModel_3_RAD.mat');
rad_models_ml(4) = computeMarginal('posteriorModel_4_RAD.mat');

atad_models_ml = zeros(4,1);
atad_models_ml(1) = computeMarginal('posteriorModel_1_ATAD.mat');
atad_models_ml(2) = computeMarginal('posteriorModel_2_ATAD.mat');
atad_models_ml(3) = computeMarginal('posteriorModel_3_ATAD.mat');
atad_models_ml(4) = computeMarginal('posteriorModel_4_ATAD.mat');

rad_models_ml_p = exp(rad_models_ml)./sum(exp(rad_models_ml));
rad_models_ml_ratios = bsxfun(@(x,y) x./y, repmat(exp(rad_models_ml),1,4), exp(rad_models_ml)');

atad_models_ml_p = exp(atad_models_ml)./sum(exp(atad_models_ml));
atad_models_ml_ratios = bsxfun(@(x,y) x./y, repmat(exp(atad_models_ml),1,4), exp(atad_models_ml)');


