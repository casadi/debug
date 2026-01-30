import casadi.*
plugins = strsplit(CasadiMeta.plugins(),';');
fid = fopen('plugins.txt','w');
for i=1:length(plugins)
    fprintf(fid,'%s\n',plugins{i});
end
fclose(fid);
disp(['Wrote ' num2str(length(plugins)) ' plugins to plugins.txt']);
