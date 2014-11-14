function [out_cell] = conjstr(in_cell)
    out_cell = in_cell;
    for ii = 1:numel(in_cell)
        out_cell{ii} = [in_cell{ii} '*'];
    end
end