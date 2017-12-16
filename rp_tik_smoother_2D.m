function [mat_out] = rp_tik_smoother_2D(mat_in,reg_param_x,reg_param_y,deriv_flag)
% function [tr_out] = rp_tik_smoother(tr_in,t,reg_param,deriv_flag)
%
% Smooth a signal using tikhonov regularization
%
% Inputs:
%       - tr_in - input signal
%       - t - input time vector
%       - reg_param - regularization paramter
%       - deriv_flag - 1st or 2nd derivative operator choice (default is 2)
%
% Outputs:
%       - tr_out - smoothed signal
%
% Made in Germany, while tired, so check it for errors
%

assert(reg_param_x>0,'please set the regularization parameter (x) to a positive number')
assert(reg_param_y>0,'please set the regularization parameter (y) to a positive number')

if nargin==3
    fprintf('no flag chosen, deriv_flag set to 2\n');
    deriv_flag=2;
end

if deriv_flag~=1 && deriv_flag~=2
    fprintf('your choice of deriv_flag makes no sense, setting to 2\n')
    deriv_flag=2;
end


[M1a,M1b] = rp_deriv_op1(mat_in,deriv_flag);

mat_out =reshape(( speye(length(mat_in(:))) + ( reg_param_y*(M1a'*M1a))  + ( reg_param_x*(M1b'*M1b)) )\mat_in(:),size(mat_in,1),size(mat_in,2));