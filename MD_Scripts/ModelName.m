% Function which generates default structure parameters
function Model_Scaled = ModelName(Model,Damp,TF_Param,S_D,S_R,S_MMD,S_MXD,S_XXD,S_E,S_S,S_A)

    % Parameter set for TF model
    if TF_Param ~= 0 && strcmp(Model,'TF')
        Model_Scaled = [Model num2str(TF_Param)];
    elseif TF_Param ~= 0 % Skip JC models
        Model_Scaled = '';
        return
    else
        Model_Scaled = Model;
    end
    
    % Damping function
    if Damp ~= 0
        Model_Scaled = [Model_Scaled 'd' num2str(Damp)];
    end

    if ~ismembertol(1.0,S_D,1e-5)
        mtxt = strrep(num2str(S_D,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_D' mtxt];
    end
    if ~ismembertol(1.0,S_R,1e-5)
        mtxt = strrep(num2str(S_R,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_R' mtxt];
    end
    
    if ~ismembertol(1.0,S_E,1e-5) && ~strcmp(Model,'TF')
        mtxt = strrep(num2str(S_E,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_E' mtxt];
    elseif ~ismembertol(1.0,S_E,1e-5) && strcmp(Model,'TF')
        Model_Scaled = ''; % epsilon scaling not defined for TF model
        return
    end
    
    if ~ismembertol(1.0,S_S,1e-5) && ~strcmp(Model,'TF')
        mtxt = strrep(num2str(S_S,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_S' mtxt];
    elseif ~ismembertol(1.0,S_S,1e-5) && strcmp(Model,'TF')
        Model_Scaled = ''; % sigma scaling not defined for TF model
        return
    end
    
    if ~ismembertol(1.0,S_A,1e-5) && strcmp(Model,'TF')
        mtxt = strrep(num2str(S_A,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_A' mtxt];
    elseif ~ismembertol(1.0,S_A,1e-5) && ~strcmp(Model,'TF')
        Model_Scaled = ''; % Alpha scaling not defined for JC model
        return
    end
    
    if ~ismembertol(1.0,S_MMD,1e-5)
        mtxt = strrep(num2str(S_MMD,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_MMD' mtxt];
    end
    if ~ismembertol(1.0,S_XXD,1e-5)
        mtxt = strrep(num2str(S_MMD,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_XXD' mtxt];
    end
    if ~ismembertol(1.0,S_MXD,1e-5)
        mtxt = strrep(num2str(S_MMD,'%10.5f'),'-','N');
        Model_Scaled = [Model_Scaled '_MXD' mtxt];
    end
    
end