[~,~,raw] = xlsread('P:\isbe\ASSURE\CANCERS_PROCAS_VOLPARA_COMBO_copy.xls');
%%
case_ids = raw(2:end, 6);
vbd_scores = cell2mat(raw(2:end, 23));
view_labels = raw(2:end, 11);
breast_labels = raw(2:end,10);

l_breast_idx = strcmpi(breast_labels, 'left');
r_breast_idx = strcmpi(breast_labels, 'right');
view_labels(l_breast_idx) = strcat('L', view_labels(l_breast_idx));
view_labels(r_breast_idx) = strcat('R', view_labels(r_breast_idx));

volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels, 'bob', 0);