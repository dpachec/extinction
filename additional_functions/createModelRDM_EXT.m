function[modelRDM] = createModelRDM_EXT(neuralRDM, allCSPIts)
    
    modelRDM = zeros(size(neuralRDM, 1));
    modelRDM(1:size(allCSPIts, 1), 1:size(allCSPIts, 1)) = 1; 

    modelRDM1(1, :, :) = modelRDM ;
    modelRDM = modelRDM1; 
end