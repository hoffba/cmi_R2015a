function con = cmiDB
% Access CMI image database

dbpath = '/mnt/cmi/projects/CMI/Animal Documentation (tumor monitoring, etc)/DataCatalog.accdb';
url = [['jdbc:odbc:Driver={Microsoft Access Driver (*.mdb, *.accdb)};DSN='';DBQ='] dbpath];
con = database('','','','sun.jdbc.odbc.JdbcOdbcDriver', url); 

