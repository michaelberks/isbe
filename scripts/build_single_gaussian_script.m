%build single gaussian model - makes use of increased memory on cluster

args.LoadIdx = 1;
args.Idx = [mberksroot, 'background/idx/normal_2/cluster_idx001'];
args.ImageDir = [mberksroot, 'background/pyramid/normal_2/'];
args.CutOffLevel = 4;
args.MaxMemory = 512;
args.SaveFile = [mberksroot, 'background/results/normal_2_single_gaussian'];
try
    [model] = mb_build_gaussian_pyramid_model(args);
catch
   %report last error
   rethrow(lasterror);
end