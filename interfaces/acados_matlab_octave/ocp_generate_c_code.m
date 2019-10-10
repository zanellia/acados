function ocp_generate_c_code(obj)
        % generate C code for CasADi functions
        if (strcmp(obj.model_struct.dyn_type, 'explicit'))
            acados_template_mex.generate_c_code_explicit_ode(obj.acados_ocp_nlp_json.model);
        elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
            if (strcmp(obj.opts_struct.sim_method, 'irk'))
                opts.generate_hess = 1;
                acados_template_mex.generate_c_code_implicit_ode(obj.acados_ocp_nlp_json.model, opts);
            end
        end

        % add checks for 
        % 1) nonlinear least-squares (jeez)
        % 2) Vz (already implemented?)
        % 3) exact Hessian
        % 4) external cost
        % 5) GNSF
        % 6) discrete dynamics
                
        % set include and lib path
        acados_folder = getenv('ACADOS_INSTALL_DIR');
        obj.acados_ocp_nlp_json.acados_include_path = [acados_folder, '/include'];
        obj.acados_ocp_nlp_json.acados_lib_path = [acados_folder, '/lib'];
        % strip non-numerical data
        
        % model
        model.name = obj.acados_ocp_nlp_json.model.name;
        
        obj.acados_ocp_nlp_json.model = [];
        obj.acados_ocp_nlp_json.model = model;
        
        con_h.name = obj.acados_ocp_nlp_json.con_h.name;
        
        obj.acados_ocp_nlp_json.con_h = [];
        obj.acados_ocp_nlp_json.con_h = con_h;
        
        con_h_e.name = obj.acados_ocp_nlp_json.con_h_e.name;
        
        obj.acados_ocp_nlp_json.con_h_e = [];
        obj.acados_ocp_nlp_json.con_h_e = con_h_e;
        
        con_p.name = obj.acados_ocp_nlp_json.con_p.name;
        
        obj.acados_ocp_nlp_json.con_p = [];
        obj.acados_ocp_nlp_json.con_p = con_p;
        
        con_p_e.name = obj.acados_ocp_nlp_json.con_p_e.name;
        
        obj.acados_ocp_nlp_json.con_p_e = [];
        obj.acados_ocp_nlp_json.con_p_e = con_p_e;
        
        % post process numerical data (mostly cast scalars to 1-dimensional cells)
        constr = obj.acados_ocp_nlp_json.constraints;
        %props = properties(constr);
        props = fieldnames(constr);
        for iprop = 1:length(props)
            thisprop = props{iprop};
            %%%Add logic here if you want to work with select properties
            thisprop_value = constr.(thisprop);
            %%%Add logic here if you want to do something based on the property's value
            if size(thisprop_value) == [1 1]
                constr.(thisprop) = num2cell(constr.(thisprop));
            end
        end
        obj.acados_ocp_nlp_json.constraints = constr;
        
        cost = obj.acados_ocp_nlp_json.cost;
        %props = properties(cost);
        props = fieldnames(cost);
        for iprop = 1:length(props)
            thisprop = props{iprop};
            %%%Add logic here if you want to work with select properties
            thisprop_value = cost.(thisprop);
            %%%Add logic here if you want to do something based on the property's value
            if norm(size(thisprop_value) - [1, 1]) == 0
                cost.(thisprop) = num2cell(cost.(thisprop));
            end
        end
        obj.acados_ocp_nlp_json.cost = cost;
        
        % load JSON layout
        acados_folder = getenv('ACADOS_INSTALL_DIR');

        acados_layout = jsondecode(fileread([acados_folder,...
            '/interfaces/acados_template/acados_template/acados_layout.json']));

        dims = obj.acados_ocp_nlp_json.dims;
        % reshape constraints
        constr = obj.acados_ocp_nlp_json.constraints;
        constr_l = acados_layout.constraints;
        fields = fieldnames(constr_l);
        for i = 1:numel(fields)
            if strcmp(constr_l.(fields{i}){1}, 'ndarray')
                if length(constr_l.(fields{i}){2}) == 1
                    this_dims = [dims.(constr_l.(fields{i}){2}{1}), 1];
                else
                    this_dims = [dims.(constr_l.(fields{i}){2}{1}), dims.(constr_l.(fields{i}){2}{1})];
                end
                constr.(fields{i}) = reshape(constr.(fields{i}), this_dims);
            end
        end
        obj.acados_ocp_nlp_json.constraints = constr;
        
        % reshape cost
        cost = obj.acados_ocp_nlp_json.cost;
        cost_l = acados_layout.cost;
        fields = fieldnames(cost_l);
        for i = 1:numel(fields)
            if strcmp(cost_l.(fields{i}){1}, 'ndarray')
                if length(cost_l.(fields{i}){2}) == 1
                    this_dims = [dims.(cost_l.(fields{i}){2}{1}), 1];
                else
                    this_dims = [dims.(cost_l.(fields{i}){2}{1}), dims.(cost_l.(fields{i}){2}{2})];
                end
                cost.(fields{i}) = reshape(cost.(fields{i}), this_dims);
                % convert 1-dimensional arrays to cells
                if length(cost_l.(fields{i}){2}) == 2 && (this_dims(1) == 1 || this_dims(2) == 1)
                    field_as_cell = {};
                    for j = 1:max(this_dims(1), this_dims(2))
                        field_as_cell{end+1} = num2cell(cost.(fields{i})(j));
                    end
                    cost.(fields{i}) = field_as_cell;
                end
            end
        end
        obj.acados_ocp_nlp_json.cost = cost;
        
        % dump JSON file
        json_string = jsonencode(obj.acados_ocp_nlp_json);
        fid = fopen('acados_ocp_nlp.json', 'w');
        if fid == -1, error('Cannot create JSON file'); end
        fwrite(fid, json_string, 'char');
        fclose(fid);
        % render templated C code
        % old call (Python + Jinja)
        % acados_template_mex.generate_solver('acados_ocp_nlp.json', '/home/andrea/.acados_t/bin/python3')
        acados_template_mex.generate_solver_matlab('acados_ocp_nlp.json')
end