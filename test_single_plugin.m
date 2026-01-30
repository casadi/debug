function test_single_plugin(plugin_str)
    import casadi.*

    parts = strsplit(plugin_str,'::');
    type = parts{1};
    name = parts{2};

    disp(['Testing plugin: ' plugin_str])

    % First load the plugin
    if strcmp(type,'XmlFile')
        % do nothing
    elseif strcmp(type,'Importer') || strcmp(type,'Linsol')
        eval([type '.load_plugin(''' name ''')'])
    else
        eval(['load_' lower(type) '(''' name ''')'])
    end

    % Then test if applicable
    if strcmp(type,'Nlpsol')
        if ismember(name, {'scpgen', 'ampl', 'qrsqp', 'fatrop'})
            disp(['Skipping nlpsol test for ' name])
            return
        end

        x = MX.sym('x');
        nlp = struct;
        nlp.x = x;
        nlp.f = x;
        nlp.g = 2*x+1;

        solver = nlpsol('solver', name, nlp);
        solver('lbg',0,'ubg',2);
        disp(['SUCCESS: nlpsol ' name])
    end

    if strcmp(type,'Conic')
        if ismember(name, {'hpipm', 'nlpsol', 'qrqp', 'ipqp', 'fatrop', 'daqp'})
            disp(['Skipping conic test for ' name])
            return
        end

        x = MX.sym('x');
        nlp = struct;
        nlp.x = x;
        nlp.f = x;
        nlp.g = 2*x+1;

        solver = qpsol('solver', name, nlp);
        solver('lbg',0,'ubg',2);
        disp(['SUCCESS: qpsol ' name])
    end

    disp(['DONE: ' plugin_str])
end
