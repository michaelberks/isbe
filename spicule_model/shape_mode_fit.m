function shape_mode_fit(mass_model, shape_idx)

    if ~isfield(mass_model, 'X_pro');
        display('Compute X_pro in mass_model first');
        return
    end
    
    for jj = 1:length(shape_idx)
        ii = shape_idx(jj);
        x_shape_all = mass_model.P_shape*mass_model.B_shape(:,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
        x_shape_1 = mass_model.P_shape(:,1)*mass_model.B_shape(1,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
        x_shape_5 = mass_model.P_shape(:,1:5)*mass_model.B_shape(1:5,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
        x_shape_10 = mass_model.P_shape(:,1:10)*mass_model.B_shape(1:10,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
        x_shape_15 = mass_model.P_shape(:,1:15)*mass_model.B_shape(1:15,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
        x_shape_20 = mass_model.P_shape(:,1:20)*mass_model.B_shape(1:20,ii)...
            + mass_model.mean_shape'/mass_model.scale_factor;
    
        figure('name', ['Mass ', num2str(ii)]); hold on;
        plot(mass_model.X_pro(ii, 1:end/2), mass_model.X_pro(ii,end/2+1:end));
        plot(x_shape_all(1:500), x_shape_all(501:end), 'r');
        plot(x_shape_1(1:500), x_shape_1(501:end), 'g:');
        plot(x_shape_5(1:500), x_shape_5(501:end), 'm:');
        plot(x_shape_10(1:500), x_shape_10(501:end), 'y:');
        plot(x_shape_20(1:500), x_shape_20(501:end), 'k:');
        plot(x_shape_15(1:500), x_shape_15(501:end), 'c:');
    end
end