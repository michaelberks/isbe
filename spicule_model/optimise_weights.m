function [weights_out, output] =...
    optimise_weights(mass_files, weights_in)    
      
    options = optimset('LargeScale', 'off',...
        'MaxIter', 100, 'MaxFunEvals', 100);
    [weights_out, fval, exitflag, output] =...
        fmincon(@myfun,weights_in,[],[],[],[], [1e-6 1e-6], [inf inf],...
            [], options);
    %weights_out = fminunc(@myfun,weights_in);
    
    function er = myfun(x)
%         [e_shape e_tex e_scale er_ss] = model_errors(mass_model, x);
%         er = mean(sum([error_weights(1)*er_ss' error_weights(2)*e_tex']));
       [e_c] = model_errors3(mass_files, x);
       er = mean(e_c.combined);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Old stuff, see model_errors.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     W_shape = weights(1);
%     W_tex   = weights(2);
%     W_scale = 1;%weights(3);
%     
%     combined_data = [W_shape*mass_model.B_shape; W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';
% 
%     [mean_c, P_c, B_c] = pca(combined_data, 0.98);
%     k_shape     = size(mass_model.B_shape, 1);
%     k_tex       = size(mass_model.B_tex, 1);
%     
%     Q_shape = P_c(1:k_shape,:); 
%     Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
%     Q_scale = P_c(end, :);
%     
%     B_shape_c = Q_shape*B_c / W_shape;
%     B_tex_c = Q_tex*B_c / W_tex;
%     B_scale_c = Q_scale*B_c / W_scale;
%     
%     er_shape = mean(sum((B_shape_c(:) - mass_model.B_shape(:)).^2))/...
%         var((B_shape_c(:) - mass_model.B_shape(:)).^2);
%     er_tex = mean(sum((B_tex_c(:) - mass_model.B_tex(:)).^2))/...
%         var((B_tex_c(:) - mass_model.B_tex(:)).^2);
%     er_scale = mean(sum((B_scale_c(:) - mass_model.B_scale(:)).^2))/...
%         var((B_scale_c(:) - mass_model.B_scale(:)).^2);
% 
%     er_shape = mean(sum((B_shape_c(:) - mass_model.B_shape(:)).^2));
%     er_tex = mean(sum((B_tex_c(:) - mass_model.B_tex(:)).^2));
%     er_scale = mean(sum((B_scale_c(:) - mass_model.B_scale(:)).^2));
    
%     er_shape = mean(sum((B_shape_c(:) - mass_model.B_shape(:)).^2))/...
%         mean(mass_model.B_shape(:).^2);
%     er_tex = mean(sum((B_tex_c(:) - mass_model.B_tex(:)).^2))/...
%         mean(mass_model.B_tex(:).^2);
%     er_scale = mean(sum((B_scale_c(:) - mass_model.B_scale(:)).^2))/...
%         mean(mass_model.B_scale(:).^2);
%     
%     er_shape = mean(sum((B_shape_c(:) - mass_model.B_shape(:)).^2))/...
%         (k_shape*mean(mass_model.B_shape(:).^2));
%     er_tex = mean(sum((B_tex_c(:) - mass_model.B_tex(:)).^2))/...
%         (k_tex*mean(mass_model.B_tex(:).^2));
%     er_scale = mean(sum((B_scale_c(:) - mass_model.B_scale(:)).^2))/...
%         mean(mass_model.B_scale(:).^2);
% 
%     
%     error = er_shape + er_tex + er_scale;
