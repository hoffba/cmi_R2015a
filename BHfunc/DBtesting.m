
% Find SQLite drivers
driverpath = dir(fullfile(fileparts(which('cmi.m')),'@DBclass','sqlite-jdbc*.jar'));
driverpath = fullfile(driverpath(end).folder,driverpath(end).name);

% Find SQLite database on server:
netname = 'maize.umhsnas.med.umich.edu\Radiology_CMI';
dbpath = 'pipeline_data\Clinical.db';
[~,r] = dos('net use');
r = strsplit(r,'\n');
ind = find(cellfun(@(x)contains(x,netname),r),1);
if isempty(ind)
    error('Cound not find database file. Please connect to Radiology_CMI network drive.');
else
    r = r{ind};
    drivebase = regexp(r,'([A-Z]*:)','tokens');
    dbpath = fullfile(drivebase{1}{1},netname,dbpath);
end

% Set up db connection:
vendor = "Other";
opts = databaseConnectionOptions("jdbc",vendor);
opts.setoptions('DataSourceName',"SQLite", ...
    'JDBCDriverLocation',driverpath, ...
    'Driver',"org.sqlite.JDBC", ...
    'URL',['jdbc:sqlite:' dbpath]);
status = testConnection(opts,'','');