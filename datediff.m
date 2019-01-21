% Calculate time difference between two date strings
function dt = datediff(str1,str2)

dt = (datenum(str2,'yyyymmdd') - datenum(str1,'yyyymmdd')) / 365;