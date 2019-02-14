%   This file is part of acados.
%
%   acados is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 3 of the License, or (at your option) any later version.
%
%   acados is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with acados; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%

function generate_c_code_nonlinear_least_squares( model, opts )

%% import casadi
import casadi.*

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

%% load model
% x
x = model.sym_x;
nx = length(x);
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
% u
u = model.sym_u;
nu = length(u);
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
	np = length(p);
else
    if isSX
        p = SX.sym('p',0, 0);
    else
        p = MX.sym('p',0, 0);
    end
	np = 0;
end

model_name = model.name;

if isfield(model, 'expr_y')
	y = model.expr_y;
	% generate jacobians
	jac_x       = jacobian(y, x);
	jac_u       = jacobian(y, u);
	% Set up functions
	if (strcmp(model.param_y, 'true'))
		y_fun_jac_ut_xt = Function([model_name,'_y_fun_jac_ut_xt'], {x, u, p}, {y, [jac_u'; jac_x']});
	else
		y_fun_jac_ut_xt = Function([model_name,'_y_fun_jac_ut_xt'], {x, u}, {y, [jac_u'; jac_x']});
	end
	% generate C code
	y_fun_jac_ut_xt.generate([model_name,'_y_fun_jac_ut_xt'], casadi_opts);
end

if isfield(model, 'expr_y_e')
	y_e = model.expr_y_e;
	% generate jacobians
	jac_x_e     = jacobian(y_e, x);
	% Set up functions
	if (strcmp(model.param_y_e, 'true'))
		y_e_fun_jac_ut_xt = Function([model_name,'_y_e_fun_jac_ut_xt'], {x, p}, {y_e, jac_x_e'});
	else
		y_e_fun_jac_ut_xt = Function([model_name,'_y_e_fun_jac_ut_xt'], {x}, {y_e, jac_x_e'});
	end
	% generate C code
	y_e_fun_jac_ut_xt.generate([model_name,'_y_e_fun_jac_ut_xt'], casadi_opts);
end

