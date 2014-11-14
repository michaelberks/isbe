args.PathToTextureGMM = 'C:/isbe/dev/background/models/dual_tree/normal512_k10_w1_9_w2_0_3_6.mat';
args.TargetImage = dual_tree{3}(:,:,6);
args.FilledImage = ones(size(args.TargetImage));
args.FilledImage(9:56,9:56) = 0;
args.ComplexMode = 1;
syn_image = mb_gmm_tex_synthesis(args);